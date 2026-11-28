# -*- coding: utf-8 -*-
"""
modelo.py - o preditor de alvos (PP4CKD): carga, fingerprint e rotulos.

Uma MLP que recebe o fingerprint de uma molecula e devolve uma probabilidade
para cada um dos 7.676 alvos do ChEMBL 36:

    SMILES -> ECFP4 (Morgan r=2, 4096 bits) -> Linear(4096,1000) -> ReLU
           -> Linear(1000,500) -> ReLU -> Linear(500,7676) -> sigmoide

Pesos pre-treinados do PP4CKD (Zenodo 10.5281/zenodo.23147411; pesos
CC-BY-4.0, codigo Apache-2.0). Dos 7.676 alvos, 4.140 sao de *Homo sapiens* e
3.536 nao sao -- **alvo nao humano e de primeira classe** aqui, e e disso que
o canal negativo da fase 3 e deduzido.

TRES COISAS QUE NAO SAO DETALHE
-------------------------------
`weights_only=True` na carga e OBRIGATORIO. Um `.pt` e um pickle, e carregar
sem a trava executaria codigo arbitrario contido no arquivo baixado.

Os INDICES do `Sequential` sao preservados (0, 3, 6), e o Dropout vira
`Identity` em vez de desaparecer: e por esses indices que o `state_dict`
publicado encaixa sem renomear nada.

**O numero de alvos vem do proprio checkpoint**, nao de uma constante. Quem
confere se checkpoint e rotulos descrevem o MESMO espaco e a etapa 1, no ponto
em que pontua -- a coluna *i* de um nao e o alvo *i* do outro.

O tarball do Zenodo esta TRUNCADO na origem: MD5 identico ao publicado, mas o
tar termina antes do fim, e so 22 dos ~77 modelos saem integros --
`ecfp4_full_model.pt` nao e um deles. O contorno e um checkpoint de fold, que e
um modelo treinado de verdade (9/10 dos dados) e **preferivel** na validacao,
porque permite medir fora do que o modelo viu.

Este modulo nunca chama `sys.exit`: ele e biblioteca. Dependencia ausente
levanta ImportError, para que a etapa 1 possa rodar em modo degradado (so com
evidencia medida) em vez de o processo hospedeiro morrer.
"""
import csv
from pathlib import Path

try:
    import torch
except ImportError as e:  # pragma: no cover - ambiente sem torch
    raise ImportError(
        "O preditor de alvos precisa do PyTorch: pip install torch. "
        "Sem ele a etapa 1 roda em modo degradado, so com evidencia medida."
    ) from e

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
except ImportError as e:  # pragma: no cover - ambiente sem rdkit
    raise ImportError(
        "O preditor de alvos precisa do RDKit: pip install rdkit"
    ) from e

# Dimensoes e nomes de arquivo conforme o PP4CKD/PPB3.
FINGERPRINTS = {
    "ECFP4": {"arquivo": "ecfp4_full_model.pt", "prefixo": "ecfp4",
              "raio": 2, "bits": 4096},
    "ECFP6": {"arquivo": "ecfp6_full_model.pt", "prefixo": "ecfp6",
              "raio": 3, "bits": 4096},
}
N_ALVOS = 7676

# Nomes comuns para os organismos mais pedidos, para que a lista de alvos possa
# dizer "rato" em vez de "Rattus norvegicus". Nao e uma taxonomia: e um atalho
# de digitacao. Qualquer trecho do nome cientifico funciona sem estar aqui.
APELIDOS_ORGANISMO = {
    "humano": ["Homo sapiens"],
    "human": ["Homo sapiens"],
    "rato": ["Rattus"],
    "camundongo": ["Mus musculus"],
    "roedor": ["Rattus", "Mus musculus", "Cavia porcellus", "Cricetulus"],
    "coelho": ["Oryctolagus"],
    "porco": ["Sus scrofa"],
    "boi": ["Bos taurus"],
    "bovino": ["Bos taurus"],
    "cao": ["Canis lupus"],
    "cachorro": ["Canis lupus"],
    "galinha": ["Gallus gallus"],
    "levedura": ["Saccharomyces"],
    "tuberculose": ["Mycobacterium tuberculosis"],
    "malaria": ["Plasmodium"],
    "chagas": ["Trypanosoma cruzi"],
    "leishmania": ["Leishmania"],
    "esquistossomose": ["Schistosoma"],
}


# =========================================================================
# Modelo
# =========================================================================

def construir_rede(bits: int, n_alvos: int) -> "torch.nn.Module":
    """A rede PPB3: duas camadas ocultas e uma saida multirrotulo.
    O Dropout nao e recriado porque so atua em treino e este modulo e de
    inferencia; os indices do Sequential sao mantidos para que o state_dict
    original encaixe sem renomear nada."""
    from torch import nn
    return nn.Sequential(
        nn.Linear(bits, 1000),      # 0
        nn.ReLU(),                  # 1
        nn.Identity(),              # 2  (era Dropout(0.2))
        nn.Linear(1000, 500),       # 3
        nn.ReLU(),                  # 4
        nn.Identity(),              # 5  (era Dropout(0.2))
        nn.Linear(500, n_alvos),    # 6
    )


def _extrair_state_dict(bruto) -> dict:
    """Aceita tanto um state_dict puro quanto um checkpoint que o embrulhe."""
    if isinstance(bruto, dict):
        for chave in ("state_dict", "model_state_dict", "model"):
            if chave in bruto and isinstance(bruto[chave], dict):
                bruto = bruto[chave]
                break
    if not isinstance(bruto, dict):
        raise RuntimeError(
            f"Formato de checkpoint nao reconhecido: {type(bruto)}")
    # Remove prefixos deixados por nn.Sequential aninhado ou DataParallel.
    limpo = {}
    for k, v in bruto.items():
        for prefixo in ("module.", "net."):
            if k.startswith(prefixo):
                k = k[len(prefixo):]
        limpo[k] = v
    return limpo


def escolher_checkpoint(diretorio: Path, cfg: dict) -> Path:
    """Decide qual .pt usar dentro da pasta de modelos.

    Prefere o modelo completo (`<fp>_full_model.pt`, treinado em 100% dos
    dados). Se ele nao estiver ali, cai para um checkpoint de fold
    (`<fp>_fold_N.pt`), que e um modelo treinado de verdade, so que em 9/10
    dos dados -- serve para previsao e e o que se deve usar na validacao, para
    medir fora do que o modelo viu.
    """
    completo = diretorio / cfg["arquivo"]
    if completo.is_file():
        return completo
    folds = sorted(diretorio.glob(f"{cfg['prefixo']}_fold_*.pt"))
    if folds:
        print(f"[AVISO] '{cfg['arquivo']}' nao esta em {diretorio}; usando o "
              f"checkpoint de fold '{folds[0].name}'. E um modelo treinado em "
              f"9/10 dos dados -- valido, e preferivel na validacao.")
        return folds[0]
    raise FileNotFoundError(
        f"Nenhum peso para este fingerprint em {diretorio}.\n"
        f"Esperado '{cfg['arquivo']}' ou '{cfg['prefixo']}_fold_N.pt'.\n"
        f"Baixe pp4ckd_models.tar.gz de https://zenodo.org/records/23147411 "
        f"e extraia. O tarball esta truncado na origem: use um fold.")


def carregar_modelo(diretorio: str | None, fp: str = "ECFP4",
                    dispositivo: str | None = None,
                    checkpoint: str | None = None):
    """Carrega os pesos pre-treinados.
    Devolve (rede, config_fp, dispositivo, nome_do_checkpoint).

    SEGURANCA: `weights_only=True` e obrigatorio aqui. Um arquivo .pt e um
    pickle, e carrega-lo sem essa trava executaria codigo arbitrario contido
    no arquivo baixado.
    """
    if fp not in FINGERPRINTS:
        raise ValueError(f"Fingerprint '{fp}' nao suportado por este modulo. "
                         f"Use um de: {', '.join(FINGERPRINTS)}")
    cfg = FINGERPRINTS[fp]
    if checkpoint:
        caminho = Path(checkpoint).expanduser()
        if not caminho.is_file():
            raise FileNotFoundError(f"Checkpoint nao encontrado: {caminho}")
    else:
        if not diretorio:
            raise ValueError(
                "Informe a pasta de modelos ou o arquivo .pt do checkpoint.")
        caminho = escolher_checkpoint(Path(diretorio).expanduser(), cfg)

    if dispositivo is None:
        dispositivo = "cuda" if torch.cuda.is_available() else "cpu"

    bruto = torch.load(caminho, map_location=dispositivo, weights_only=True)
    state = _extrair_state_dict(bruto)

    # O numero de alvos vem do PROPRIO checkpoint, nao de uma constante.
    n_alvos = N_ALVOS
    for chave in ("6.weight", "net.6.weight"):
        if chave in state:
            n_alvos = state[chave].shape[0]
            break

    rede = construir_rede(cfg["bits"], n_alvos)
    faltando, sobrando = rede.load_state_dict(state, strict=False)
    criticas = [k for k in faltando if k.endswith(("weight", "bias"))]
    if criticas:
        raise RuntimeError(f"Pesos ausentes no checkpoint: {criticas}")
    if sobrando:
        print(f"[AVISO] Chaves ignoradas no checkpoint: {list(sobrando)[:5]}")

    rede.to(dispositivo).eval()
    return rede, cfg, dispositivo, caminho.name


# =========================================================================
# Fingerprints
# =========================================================================

def fingerprint_de(smiles: str, cfg: dict):
    """Vetor de bits Morgan, ou None se o RDKit nao ler o SMILES."""
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        return None
    gerador = rdFingerprintGenerator.GetMorganGenerator(
        radius=cfg["raio"], fpSize=cfg["bits"])
    return gerador.GetFingerprintAsNumPy(mol)


# =========================================================================
# Alvos
# =========================================================================

def carregar_rotulos(caminho: str) -> list[dict]:
    """Le target_labels.tsv -> lista indexada por `target_idx`.

    A posicao na lista E a coluna da saida da rede: o indice *i* tem de ser o
    alvo *i*, senao a probabilidade de um alvo leva o nome de outro.
    """
    caminho = Path(caminho).expanduser()
    if not caminho.is_file():
        raise FileNotFoundError(f"Rotulos de alvo nao encontrados: {caminho}")
    por_indice: dict[int, dict] = {}
    with open(caminho, "r", encoding="utf-8-sig", newline="") as f:
        for linha in csv.DictReader(f, delimiter="\t"):
            try:
                idx = int(linha["target_idx"])
            except (KeyError, TypeError, ValueError):
                continue
            try:
                n_ativos = int(float(linha.get("active_molecules") or 0))
            except (TypeError, ValueError):
                n_ativos = 0
            por_indice[idx] = {
                # Normalizado AQUI, no ponto de leitura: o ID do ChEMBL e
                # maiusculo, e quem compara a jusante nao repete a regra.
                "chembl_id": (linha.get("chembl_id") or "").strip().upper(),
                "nome": (linha.get("pref_name") or "").strip(),
                "organismo": (linha.get("organism") or "").strip(),
                "tipo": (linha.get("target_type") or "").strip(),
                # Quantas moleculas ativas treinaram ESTE alvo. Entra na saida
                # porque e a unica forma de o pesquisador ver que uma
                # probabilidade de 0,9 sobre um alvo com 8 ativos nao vale o
                # mesmo que sobre um alvo com 4.000.
                "n_ativos": n_ativos,
            }
    if not por_indice:
        raise RuntimeError(f"Nenhum alvo lido de {caminho}")
    return [por_indice.get(i, {"chembl_id": "", "nome": f"idx_{i}",
                               "organismo": "", "tipo": "", "n_ativos": 0})
            for i in range(max(por_indice) + 1)]


def dividir_organismos(valores) -> list[str]:
    """Normaliza o que veio da opcao de organismo (repetida e/ou com virgulas)
    numa lista achatada. Devolve [] quando nao ha filtro."""
    if not valores:
        return []
    if isinstance(valores, str):
        valores = [valores]
    saida = []
    for v in valores:
        for parte in str(v).split(","):
            parte = parte.strip()
            if parte:
                saida.append(parte)
    return saida


def expandir_organismo(termo: str) -> list[str]:
    """Apelido -> nomes cientificos. Termo que nao e apelido volta como esta,
    porque qualquer trecho do nome cientifico funciona."""
    return APELIDOS_ORGANISMO.get(termo.strip().lower(), [termo.strip()])
