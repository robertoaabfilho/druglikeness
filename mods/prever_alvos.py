#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
prever_alvos.py

Responde à pergunta central deste pipeline:

    dada uma lista de compostos e uma lista de alvos proteicos de interesse,
    quais compostos têm chance de atuar sobre esses alvos?

Usa os pesos pré-treinados do PP4CKD (Zenodo 10.5281/zenodo.23147411,
CC-BY-4.0; código Apache-2.0), uma rede PPB3 que recebe o fingerprint de uma
molécula e devolve uma probabilidade para cada um dos 7.676 alvos do ChEMBL 36.

    SMILES -> ECFP4 (Morgan r=2, 4096 bits) -> Linear(4096,1000) -> ReLU
           -> Linear(1000,500) -> ReLU -> Linear(500,7676) -> sigmoide

Por que reimplementar o forward em vez de importar o pacote `pp4ckd`: ele pede
`numpy<2.0` e `mhfp`, e o modelo é uma MLP de três camadas cujos pesos são um
`state_dict` comum. Carregar o `state_dict` direto evita o conflito de versões
e deixa a dependência em só torch + rdkit.

O que este módulo NÃO faz: decidir a biologia da doença. A lista de alvos é
ENTRADA. Quem a fornece (ou a enriquece a partir da doença) é o pesquisador ou
a etapa de enriquecimento -- aqui ela é respeitada como dada.

ALVOS DE QUALQUER ORGANISMO
---------------------------
Dos 7.676 alvos do modelo, 4.140 são de Homo sapiens e 3.536 não são: 666 de
camundongo, 646 de rato, e centenas de patógeno -- Mycobacterium tuberculosis,
Plasmodium falciparum, Trypanosoma cruzi, Staphylococcus aureus, Escherichia
coli, entre 1.084 organismos. Alvo não humano é caso de primeira classe aqui,
por dois usos distintos:

  afirmativo  doença infecciosa ou parasitária. O alvo é do patógeno, e o
              modelo traz inclusive alvos de ORGANISMO INTEIRO (ensaio
              fenotípico), como CHEMBL360 = M. tuberculosis com 11.049
              moléculas ativas no treino.

  negativo    filtro de atividade indesejada. Medido no conjunto-controle
              deste repositório (6 fármacos que tratam dor neuropática contra
              6 que a causam), a AUC sobe de 0,750 com alvos humanos para
              0,806 ao multiplicar por (1 - p_bacteriano): o ciprofloxacino
              cai do 4º para o 7º lugar, porque atua em DNA girase e não em
              alvo de dor.

DUAS ARMADILHAS, as duas tratadas no código
-------------------------------------------
1. ORTÓLOGO SILENCIOSO. 'noradrenaline transporter' existe em humano
   (CHEMBL222), camundongo (CHEMBL2370) e rato (CHEMBL304, sob outro nome).
   Pedir por nome sem dizer o organismo traz os três. Por isso há
   `--organismo`, o sufixo `| organismo` por linha, e a coluna `organismo`
   em toda linha da saída. Pedir por ID ChEMBL vence o filtro, porque um ID
   é uma escolha explícita.

2. ALVO POUCO TREINADO. Alvo não humano é bem mais esparso: mediana de 26
   moléculas ativas contra 67 dos humanos, 67% deles abaixo de 50 e 798
   abaixo de 10. Uma probabilidade de 0,9 sobre um alvo treinado com 8
   moléculas é ruído com aparência de número. Por isso a saída traz
   `moleculas_treino` e `--minimo-treino` avisa (sem remover) quando o alvo
   pedido está abaixo do limiar.

LIMITE CONHECIDO: PRÓ-FÁRMACO. O modelo é baseado em ligante, isto é, julga a
estrutura que recebe. A etionamida pontua 0,80 em M. tuberculosis
(célula inteira, que é fenótipo) e 0,0017 em InhA (a proteína), porque ela só
inibe a InhA depois de ativada pela EthA -- a molécula-mãe não se liga. Mesma
razão pela qual gabapentina e pregabalina, aminoácidos pequenos e polares,
pontuam baixo nos alvos de dor. Não é defeito de implementação; é o que um
modelo de ligante pode e não pode afirmar.

PROCEDÊNCIA DO RESULTADO
------------------------
Toda linha da saída diz de onde veio a afirmação:

    predito (PP4CKD p=0.93)     previsão de modelo, hipótese a testar
    medido (ChEMBL pchembl=7.2) bioatividade experimental

As duas NUNCA saem na mesma coluna sem distinção: uma previsão de 0,93 não
tem o mesmo peso que um ensaio. Este módulo produz as previsões; a evidência
medida entra por `--evidencias`, quando houver.

USO
---
    # baixar e extrair os pesos (uma vez):
    #   https://zenodo.org/records/23147411 -> pp4ckd_models.tar.gz
    #   tar -xzf pp4ckd_models.tar.gz -C ~/modelos_baixados/pp4ckd/

    python mods/prever_alvos.py \
        --entrada input.exemplo.tsv \
        --modelos ~/modelos_baixados/pp4ckd/models \
        --rotulos ~/modelos_baixados/pp4ckd/target_labels.tsv \
        --alvos alvos.txt \
        --saida results/alvos_preditos.tsv

`--alvos` é um arquivo de texto com um alvo por linha, aceitando ID de alvo do
ChEMBL (`CHEMBL230`) ou trecho do nome preferencial (`cyclooxygenase`).
Linhas vazias e iniciadas por '#' são ignoradas, e '#' também vale como
comentário no fim da linha. Sem `--alvos`, devolve os `--top-k` alvos mais
prováveis de cada composto -- de qualquer organismo, salvo se `--organismo`
restringir.

    # alvos do patógeno (doença infecciosa)
    CHEMBL360                      # M. tuberculosis, célula inteira
    CHEMBL1849                     # InhA de Mtb
    DNA gyrase | Escherichia coli  # nome + organismo fixado nesta linha
    dihydrofolate reductase | tuberculose

    # ou restrinja tudo de uma vez:
    #   --organismo humano
    #   --organismo Plasmodium          (alcança as cepas e isolados)
    #   --organismo rato,camundongo

Para descobrir a grafia exata, sem rodar previsão nenhuma:

    python mods/prever_alvos.py --rotulos target_labels.tsv --listar-organismos
    python mods/prever_alvos.py --rotulos target_labels.tsv \
        --listar-organismos Plasmodium
    python mods/prever_alvos.py --rotulos target_labels.tsv \
        --buscar-alvo gyrase --organismo Staphylococcus

NOTA SOBRE OS PESOS: o pp4ckd_models.tar.gz publicado no Zenodo está
truncado na origem (o MD5 confere com o publicado, mas o tar termina antes
do fim: só 22 dos ~77 modelos saem íntegros, e `ecfp4_full_model.pt` não é
um deles). Por isso `--checkpoint` aceita um .pt específico e a pasta
`--modelos` cai automaticamente para `<fp>_fold_N.pt`, que é um modelo
treinado em 9/10 dos dados -- válido para previsão e preferível na validação.

REQUISITOS
----------
    pip install torch rdkit
"""

import argparse
import csv
import sys
from pathlib import Path

try:
    import torch
except ImportError:
    print("[ERRO] Este módulo precisa do PyTorch: pip install torch", file=sys.stderr)
    sys.exit(1)

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
except ImportError:
    print("[ERRO] Este módulo precisa do RDKit: pip install rdkit", file=sys.stderr)
    sys.exit(1)

# Dimensões e nomes de arquivo conforme o PP4CKD/PPB3.
FINGERPRINTS = {
    "ECFP4": {"arquivo": "ecfp4_full_model.pt", "prefixo": "ecfp4",
              "raio": 2, "bits": 4096},
    "ECFP6": {"arquivo": "ecfp6_full_model.pt", "prefixo": "ecfp6",
              "raio": 3, "bits": 4096},
}
N_ALVOS = 7676

# Nomes comuns para os organismos mais pedidos, para que a lista de alvos e a
# opção --organismo possam dizer "rato" em vez de "Rattus norvegicus". Não é
# uma taxonomia: é um atalho de digitação. Qualquer trecho do nome científico
# funciona sem estar aqui.
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
    """A rede PPB3: duas camadas ocultas e uma saída multirrótulo.
    O Dropout não é recriado porque só atua em treino e este módulo é de
    inferência; os índices do Sequential são mantidos para que o state_dict
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
        raise RuntimeError(f"Formato de checkpoint não reconhecido: {type(bruto)}")
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
    dados). Se ele não estiver ali, cai para um checkpoint de fold
    (`<fp>_fold_N.pt`), que é um modelo treinado de verdade, só que em 9/10
    dos dados -- serve para previsão e é o que se deve usar na validação, para
    medir fora do que o modelo viu.
    """
    completo = diretorio / cfg["arquivo"]
    if completo.is_file():
        return completo
    folds = sorted(diretorio.glob(f"{cfg['prefixo']}_fold_*.pt"))
    if folds:
        print(f"[AVISO] '{cfg['arquivo']}' não está em {diretorio}; usando o "
              f"checkpoint de fold '{folds[0].name}'. É um modelo treinado em "
              f"9/10 dos dados -- válido, e preferível na validação.")
        return folds[0]
    raise FileNotFoundError(
        f"Nenhum peso para este fingerprint em {diretorio}.\n"
        f"Esperado '{cfg['arquivo']}' ou '{cfg['prefixo']}_fold_N.pt'.\n"
        f"Baixe pp4ckd_models.tar.gz de https://zenodo.org/records/23147411 "
        f"e extraia.")


def carregar_modelo(diretorio: str | None, fp: str = "ECFP4",
                    dispositivo: str | None = None,
                    checkpoint: str | None = None):
    """Carrega os pesos pré-treinados.
    Devolve (rede, config_fp, dispositivo, nome_do_checkpoint).

    SEGURANÇA: `weights_only=True` é obrigatório aqui. Um arquivo .pt é um
    pickle, e carregá-lo sem essa trava executaria código arbitrário contido
    no arquivo baixado.
    """
    if fp not in FINGERPRINTS:
        raise ValueError(f"Fingerprint '{fp}' não suportado por este módulo. "
                         f"Use um de: {', '.join(FINGERPRINTS)}")
    cfg = FINGERPRINTS[fp]
    if checkpoint:
        caminho = Path(checkpoint).expanduser()
        if not caminho.is_file():
            raise FileNotFoundError(f"Checkpoint não encontrado: {caminho}")
    else:
        if not diretorio:
            raise ValueError("Informe --modelos (pasta) ou --checkpoint (arquivo .pt).")
        caminho = escolher_checkpoint(Path(diretorio).expanduser(), cfg)

    if dispositivo is None:
        dispositivo = "cuda" if torch.cuda.is_available() else "cpu"

    bruto = torch.load(caminho, map_location=dispositivo, weights_only=True)
    state = _extrair_state_dict(bruto)

    # O número de alvos vem do próprio checkpoint, não de uma constante.
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
    """Vetor de bits Morgan, ou None se o RDKit não ler o SMILES."""
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
    """Lê target_labels.tsv -> lista indexada por target_idx."""
    caminho = Path(caminho).expanduser()
    if not caminho.is_file():
        raise FileNotFoundError(f"Rótulos de alvo não encontrados: {caminho}")
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
                # Quantas moléculas ativas treinaram ESTE alvo. Entra na saída
                # porque é a única forma de o pesquisador ver que uma
                # probabilidade de 0,9 sobre um alvo com 8 ativos não vale o
                # mesmo que sobre um alvo com 4.000.
                "n_ativos": n_ativos,
            }
    if not por_indice:
        raise RuntimeError(f"Nenhum alvo lido de {caminho}")
    return [por_indice.get(i, {"chembl_id": "", "nome": f"idx_{i}",
                               "organismo": "", "tipo": "", "n_ativos": 0})
            for i in range(max(por_indice) + 1)]


def dividir_organismos(valores) -> list[str]:
    """Normaliza o que veio de --organismo (repetida e/ou com vírgulas) numa
    lista achatada. Devolve [] quando não há filtro."""
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
    """Um filtro de organismo -> os trechos que ele representa.
    'rato' -> ['Rattus']; 'Plasmodium' -> ['Plasmodium'] (passa direto)."""
    return APELIDOS_ORGANISMO.get(termo.casefold().strip(), [termo])


def casa_organismo(organismo_do_rotulo: str, filtros: list[str]) -> bool:
    """True quando o organismo do rótulo satisfaz algum filtro pedido.

    A comparação é por TRECHO, de propósito: as bases registram cepa e isolado
    dentro do nome do organismo, então 'Mycobacterium tuberculosis' precisa
    alcançar 'Mycobacterium tuberculosis (strain ATCC 25618 / H37Rv)' e
    'Plasmodium' precisa alcançar 'Plasmodium falciparum (isolate 3D7)'.
    Sem filtro, qualquer organismo serve.
    """
    if not filtros:
        return True
    alvo = organismo_do_rotulo.casefold()
    for f in filtros:
        for trecho in expandir_organismo(f):
            if trecho.casefold() in alvo:
                return True
    return False


def ler_lista_de_alvos(caminho: str) -> list[tuple[str, list[str]]]:
    """Um alvo por linha; devolve pares (termo, organismos_da_linha).

    '#' inicia comentário, na linha inteira ou depois do alvo --
    `CHEMBL1919  # canal de cálcio` vale `CHEMBL1919`, para que a lista possa
    documentar por que cada alvo está ali.

    Um organismo pode ser fixado por linha depois de '|', e vale só para
    aquela linha (sobrepõe --organismo):

        CHEMBL364                        # P. falciparum, ensaio de célula inteira
        enoyl reductase | Staphylococcus aureus
        cyclooxygenase  | humano
        transporter     | rato, camundongo

    O separador é '|' e não ',' nem ';' porque nome preferencial de alvo
    contém vírgula, parêntese e barra com frequência ("Voltage-dependent
    calcium channel subunit alpha-2/delta-1").
    """
    pedidos: list[tuple[str, list[str]]] = []
    with open(Path(caminho).expanduser(), "r", encoding="utf-8-sig") as f:
        for linha in f:
            linha = linha.split("#", 1)[0].strip()
            if not linha:
                continue
            termo, sep, organismo = linha.partition("|")
            termo = termo.strip()
            if not termo:
                continue
            pedidos.append((termo, dividir_organismos(organismo) if sep else []))
    return pedidos


def resolver_alvos(pedidos, rotulos: list[dict],
                   organismos_padrao: list[str] | None = None,
                   minimo_treino: int = 0) -> tuple[list[int], list[str]]:
    """Casa cada alvo pedido com os índices do modelo.

    `pedidos` são os pares (termo, organismos) de `ler_lista_de_alvos`; uma
    lista de strings simples também é aceita. Quando a linha não traz
    organismo, vale `organismos_padrao` (a opção --organismo).

    Duas formas de pedir um alvo:

      ID ChEMBL  casamento exato, e o filtro de organismo NÃO se aplica. Se
                 você escreveu CHEMBL304, você quer o transportador de
                 noradrenalina de rato -- pedir por ID é dizer exatamente qual.

      Nome       casamento por trecho no nome preferencial, aí o filtro de
                 organismo se aplica. É aqui que ele importa: 'noradrenaline
                 transporter' existe em humano, rato e camundongo, e sem
                 filtro os três entram.

    Devolve (índices, não encontrados).
    """
    por_chembl = {r["chembl_id"].upper(): i for i, r in enumerate(rotulos) if r["chembl_id"]}
    organismos_padrao = organismos_padrao or []
    indices: list[int] = []
    nao_encontrados: list[str] = []

    for pedido in pedidos:
        if isinstance(pedido, (tuple, list)):
            termo, organismos_linha = pedido[0], list(pedido[1] or [])
        else:
            termo, organismos_linha = pedido, []
        t = str(termo).strip()
        if not t:
            continue
        filtros = organismos_linha or organismos_padrao

        if t.upper() in por_chembl:
            achados = [por_chembl[t.upper()]]
            r = rotulos[achados[0]]
            # Avisa, sem bloquear: o ID é uma escolha explícita do pesquisador.
            if filtros and not casa_organismo(r["organismo"], filtros):
                print(f"[AVISO] {t} é de '{r['organismo'] or 'organismo não informado'}', "
                      f"fora do filtro {filtros}. Mantido: pedido por ID ChEMBL vence "
                      f"o filtro de organismo.")
        else:
            alvo = t.casefold()
            candidatos = [i for i, r in enumerate(rotulos) if alvo in r["nome"].casefold()]
            achados = [i for i in candidatos
                       if casa_organismo(rotulos[i]["organismo"], filtros)]
            descartados = len(candidatos) - len(achados)
            if descartados:
                print(f"[INFO] '{termo}': {descartados} alvo(s) de outro organismo "
                      f"descartados pelo filtro {filtros}.")
            if len(achados) > 1:
                por_org: dict[str, int] = {}
                for i in achados:
                    org = rotulos[i]["organismo"] or "(sem organismo)"
                    por_org[org] = por_org.get(org, 0) + 1
                resumo = "; ".join(f"{o}: {n}" for o, n in
                                   sorted(por_org.items(), key=lambda kv: -kv[1])[:6])
                print(f"[AVISO] '{termo}' casou com {len(achados)} alvos pelo nome "
                      f"({resumo}). Use o ID ChEMBL, ou '{termo} | <organismo>', "
                      f"para escolher.")

        if achados:
            indices.extend(achados)
        else:
            nao_encontrados.append(t)

    # Preserva a ordem de pedido, sem repetir.
    vistos, unicos = set(), []
    for i in indices:
        if i not in vistos:
            vistos.add(i)
            unicos.append(i)

    if minimo_treino > 0:
        fracos = [i for i in unicos if rotulos[i]["n_ativos"] < minimo_treino]
        if fracos:
            print(f"[AVISO] {len(fracos)} alvo(s) pedidos foram treinados com menos de "
                  f"{minimo_treino} moléculas ativas -- a previsão neles é frágil:")
            for i in sorted(fracos, key=lambda j: rotulos[j]["n_ativos"])[:10]:
                r = rotulos[i]
                print(f"         {r['chembl_id']:<14} {r['nome'][:44]:<44} "
                      f"{r['organismo'][:26]:<26} {r['n_ativos']} ativos")
    return unicos, nao_encontrados


# =========================================================================
# Exploração: sem isto não há como escrever uma lista de alvos não humanos,
# porque o modelo cobre 1.084 organismos e a grafia tem de ser a do ChEMBL.
# =========================================================================

def listar_organismos(rotulos: list[dict], filtro: str = "", limite: int = 60) -> None:
    """Organismos presentes no modelo, com quantos alvos e quanto treino."""
    agregado: dict[str, list[int]] = {}
    for r in rotulos:
        org = r["organismo"] or "(sem organismo)"
        if filtro and filtro.casefold() not in org.casefold():
            continue
        acc = agregado.setdefault(org, [0, 0])
        acc[0] += 1
        acc[1] += r["n_ativos"]
    if not agregado:
        print(f"[INFO] Nenhum organismo contém '{filtro}'.")
        return
    ordenado = sorted(agregado.items(), key=lambda kv: -kv[1][0])
    titulo = f" contendo '{filtro}'" if filtro else ""
    print(f"[INFO] {len(agregado)} organismo(s){titulo} no modelo "
          f"(mostrando até {limite}):\n")
    print(f"{'organismo':<52} {'alvos':>6} {'moléculas ativas':>18}")
    print("-" * 78)
    for org, (n_alvos, n_ativos) in ordenado[:limite]:
        print(f"{org[:52]:<52} {n_alvos:>6} {n_ativos:>18,}")
    if len(ordenado) > limite:
        print(f"... e {len(ordenado) - limite} outro(s). Refine com um trecho do nome.")


def buscar_alvos(rotulos: list[dict], termo: str,
                 organismos: list[str] | None = None, limite: int = 60) -> None:
    """Alvos cujo nome ou ID contém `termo`, respeitando o filtro de organismo."""
    t = termo.casefold()
    achados = [(i, r) for i, r in enumerate(rotulos)
               if (t in r["nome"].casefold() or t == r["chembl_id"].casefold())
               and casa_organismo(r["organismo"], organismos or [])]
    if not achados:
        extra = f" com organismo em {organismos}" if organismos else ""
        print(f"[INFO] Nenhum alvo com '{termo}'{extra}. "
              f"Tente --listar-organismos para conferir a grafia.")
        return
    achados.sort(key=lambda par: -par[1]["n_ativos"])
    print(f"[INFO] {len(achados)} alvo(s) com '{termo}' (mostrando até {limite}), "
          f"do mais para o menos treinado:\n")
    print(f"{'chembl_id':<14} {'alvo':<44} {'organismo':<30} {'tipo':<16} {'ativos':>8}")
    print("-" * 116)
    for _, r in achados[:limite]:
        print(f"{r['chembl_id']:<14} {r['nome'][:44]:<44} {r['organismo'][:30]:<30} "
              f"{r['tipo'][:16]:<16} {r['n_ativos']:>8,}")
    if len(achados) > limite:
        print(f"... e {len(achados) - limite} outro(s).")


# =========================================================================
# Previsão
# =========================================================================

def prever(registros: list[dict], rede, cfg: dict, dispositivo: str,
           lote: int = 256) -> tuple[list[dict], "torch.Tensor"]:
    """Devolve (registros_validos, probabilidades) com uma linha de
    probabilidades por registro válido."""
    validos, vetores = [], []
    for r in registros:
        fp = fingerprint_de(r.get("smiles", ""), cfg)
        if fp is None:
            continue
        validos.append(r)
        vetores.append(fp)

    if not validos:
        return [], torch.empty(0)

    import numpy as np
    saidas = []
    with torch.no_grad():
        for i in range(0, len(vetores), lote):
            bloco = np.asarray(vetores[i:i + lote], dtype="float32")
            x = torch.from_numpy(bloco).to(dispositivo)
            saidas.append(torch.sigmoid(rede(x)).cpu())
    return validos, torch.cat(saidas, dim=0)


def main():
    parser = argparse.ArgumentParser(
        description="Prevê, para cada composto da entrada, a probabilidade de atuar "
                    "sobre os alvos de interesse (pesos do PP4CKD). Os alvos podem ser "
                    "de qualquer organismo coberto pelo modelo, não só humanos.")
    parser.add_argument("--entrada", default=None,
                        help="Arquivo de entrada do pipeline (TSV com cabeçalho '#').")
    parser.add_argument("--modelos", default=None,
                        help="Pasta com os .pt extraídos de pp4ckd_models.tar.gz.")
    parser.add_argument("--checkpoint", default=None,
                        help="Arquivo .pt específico, em vez de deixar a pasta "
                             "--modelos escolher. Use para fixar um fold na validação.")
    parser.add_argument("--rotulos", required=True, help="target_labels.tsv do PP4CKD.")
    parser.add_argument("--alvos", default=None,
                        help="Arquivo com os alvos de interesse (um por linha; "
                             "'termo | organismo' fixa o organismo daquela linha). "
                             "Sem ele, reporta os --top-k alvos de cada composto.")
    parser.add_argument("--organismo", action="append", default=None, metavar="ORG",
                        help="Restringe os alvos casados POR NOME a estes organismos. "
                             "Repita a opção ou separe por vírgula. Aceita apelido "
                             "('humano', 'rato', 'tuberculose', 'malaria') e trecho do "
                             "nome científico ('Plasmodium', que alcança as cepas). "
                             "Um alvo pedido por ID ChEMBL ignora o filtro. SEM a "
                             "opção, alvo de qualquer organismo é aceito.")
    parser.add_argument("--minimo-treino", type=int, default=50, metavar="N",
                        help="Avisa (não remove) quando um alvo pedido foi treinado com "
                             "menos de N moléculas ativas. Padrão 50: dois terços dos "
                             "alvos não humanos estão abaixo disso. Use 0 para calar.")
    parser.add_argument("--listar-organismos", nargs="?", const="", default=None,
                        metavar="TRECHO",
                        help="Só lista os organismos do modelo (filtrando por trecho, "
                             "se dado) e sai. Use para achar a grafia exata do ChEMBL.")
    parser.add_argument("--buscar-alvo", default=None, metavar="TERMO",
                        help="Só lista os alvos cujo nome contém TERMO e sai. "
                             "Respeita --organismo. Use para montar a lista de alvos.")
    parser.add_argument("--saida", default="alvos_preditos.tsv", help="TSV de saída.")
    parser.add_argument("--fp", default="ECFP4", choices=sorted(FINGERPRINTS),
                        help="Fingerprint/modelo a usar (padrão: ECFP4).")
    parser.add_argument("--top-k", type=int, default=10,
                        help="Quantos alvos reportar por composto quando não há --alvos.")
    parser.add_argument("--corte", type=float, default=0.0,
                        help="Descarta previsões com probabilidade abaixo deste valor.")
    args = parser.parse_args()

    rotulos = carregar_rotulos(args.rotulos)
    organismos = dividir_organismos(args.organismo)

    # Modos de exploração: só precisam dos rótulos, não dos pesos nem da entrada.
    if args.listar_organismos is not None:
        listar_organismos(rotulos, args.listar_organismos)
        return
    if args.buscar_alvo:
        buscar_alvos(rotulos, args.buscar_alvo, organismos)
        return

    ausentes = []
    if not args.entrada:
        ausentes.append("--entrada")
    if not args.modelos and not args.checkpoint:
        ausentes.append("--modelos ou --checkpoint")
    if ausentes:
        parser.error("faltam argumentos para prever: " + ", ".join(ausentes))

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from farmaco_completo import ler_entrada

    registros = ler_entrada(args.entrada)
    if not registros:
        print(f"[ERRO] Nenhum composto lido de {args.entrada}", file=sys.stderr)
        sys.exit(1)
    print(f"[INFO] {len(registros)} composto(s) lidos de '{args.entrada}'.")

    n_humanos = sum(1 for r in rotulos if r["organismo"] == "Homo sapiens")
    print(f"[INFO] {len(rotulos)} alvos no modelo "
          f"({n_humanos} de Homo sapiens, {len(rotulos) - n_humanos} de outros "
          f"organismos).")
    if organismos:
        print(f"[INFO] Filtro de organismo para alvos pedidos por nome: {organismos}")

    rede, cfg, dispositivo, nome_ckpt = carregar_modelo(
        args.modelos, args.fp, checkpoint=args.checkpoint)
    print(f"[INFO] Modelo {args.fp} ('{nome_ckpt}') carregado em '{dispositivo}'.")

    indices_pedidos = None
    if args.alvos:
        pedidos = ler_lista_de_alvos(args.alvos)
        indices_pedidos, nao_achados = resolver_alvos(
            pedidos, rotulos, organismos, args.minimo_treino)
        print(f"[INFO] {len(pedidos)} alvo(s) pedidos -> {len(indices_pedidos)} "
              f"encontrados no modelo.")
        if nao_achados:
            print(f"[AVISO] Sem correspondência no modelo ({len(nao_achados)}): "
                  f"{', '.join(nao_achados[:8])}"
                  f"{'...' if len(nao_achados) > 8 else ''}")
            print("        Procure a grafia com --buscar-alvo TERMO "
                  "ou --listar-organismos TRECHO.")
        if not indices_pedidos:
            print("[ERRO] Nenhum dos alvos pedidos existe no modelo. Encerrando.",
                  file=sys.stderr)
            sys.exit(1)
        orgs_resolvidos = sorted({rotulos[i]["organismo"] or "(sem organismo)"
                                  for i in indices_pedidos})
        print(f"[INFO] Organismos nos alvos resolvidos: {', '.join(orgs_resolvidos)}")

    validos, probs = prever(registros, rede, cfg, dispositivo)
    invalidos = len(registros) - len(validos)
    if invalidos:
        print(f"[AVISO] {invalidos} composto(s) com SMILES ilegível, ignorados.")

    colunas = ["nome", "inchikey_parental", "chembl_id_composto", "alvo_chembl_id",
               "alvo", "organismo", "tipo_alvo", "moleculas_treino",
               "probabilidade", "procedencia"]
    linhas = 0
    saida = Path(args.saida)
    if saida.parent != Path(""):
        saida.parent.mkdir(parents=True, exist_ok=True)
    with open(saida, "w", encoding="utf-8", newline="") as f:
        f.write("#" + "\t".join(colunas) + "\n")
        escritor = csv.DictWriter(f, fieldnames=colunas, delimiter="\t",
                                  lineterminator="\n", extrasaction="ignore")
        for posicao, registro in enumerate(validos):
            linha_probs = probs[posicao]
            if indices_pedidos is not None:
                selecao = sorted(indices_pedidos,
                                 key=lambda i: float(linha_probs[i]), reverse=True)
            else:
                selecao = torch.topk(linha_probs, min(args.top_k, len(linha_probs)))\
                    .indices.tolist()
            for i in selecao:
                p = float(linha_probs[i])
                if p < args.corte:
                    continue
                escritor.writerow({
                    "nome": registro.get("nome", ""),
                    "inchikey_parental": registro.get("inchikey_parental", ""),
                    "chembl_id_composto": registro.get("chembl_id", ""),
                    "alvo_chembl_id": rotulos[i]["chembl_id"],
                    "alvo": rotulos[i]["nome"],
                    "organismo": rotulos[i]["organismo"],
                    # ORGANISM/CELL-LINE não é proteína: é ensaio fenotípico de
                    # célula inteira. Afirmação diferente, tem de sair rotulada.
                    "tipo_alvo": rotulos[i]["tipo"],
                    "moleculas_treino": rotulos[i]["n_ativos"],
                    "probabilidade": f"{p:.4f}",
                    "procedencia": f"predito (PP4CKD {args.fp} {nome_ckpt})",
                })
                linhas += 1

    print(f"[INFO] {linhas} previsão(ões) gravadas em '{saida}'.")


if __name__ == "__main__":
    main()
