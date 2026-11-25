#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
importar_prospector.py

Converte a tabela de compostos do pipeline de prospecção de plantas
(`resultados/atividades_compostos_parcial.csv`) no formato de entrada deste
pipeline (TSV com cabeçalho), aplicando o MESMO tratamento químico que a
prospecção usa, para que as duas bases possam ser cruzadas depois.

Os dois projetos são repositórios SEPARADOS: este script não importa nada
do código da prospecção. Ele só conhece o nome das colunas que a prospecção
escreve, documentados abaixo, e as reimplementa de forma independente.

COMO A PROSPECÇÃO OBTÉM E TRATA OS DADOS (contrato replicado aqui)
-------------------------------------------------------------------
Nome       `composto_nome`, com `nome_pubchem` como alternativa. Nomes de
           produto natural costumam conter vírgula ("1,8-cineol",
           "3,3'-diindolilmetano") -- 42,8% deles na base real --, por isso
           a saída é TSV e nunca CSV.

SMILES     A prospecção coleta o SMILES em até quatro bases e escolhe o da
           primeira fonte cujo SMILES o RDKit consegue ler, nesta ordem de
           prioridade:

               NPASS > COCONUT > LOTUS > PubChem

           O escolhido fica em `smiles`, a fonte em `fonte_smiles`, e os
           candidatos de cada base continuam em `smiles_npass`,
           `smiles_coconut`, `smiles_lotus` e `smiles_pubchem`.

Conferência
           Quando duas ou mais bases têm SMILES para o mesmo composto, a
           prospecção compara as ESTRUTURAS canônicas (não o texto) e grava
           o veredito em `conferencia_smiles`:

               "mesma molécula"
               "estereoquímica diferente: NPASS/PubChem x COCONUT"
               "estruturas diferentes: ..."
               "SMILES inválido: <fonte>"

           Compostos cujo veredito aponta divergência estrutural entram na
           saída marcados, não descartados -- a decisão é do pesquisador.

InChIKey   `InChIKey` é a chave de identidade declarada pelas bases. Aqui
           recalculamos duas chaves a partir do SMILES:

             inchikey            da molécula como está (preserva sal e carga)
             inchikey_parental   depois de Cleanup -> maior fragmento
                                 orgânico -> neutralização

           A chave parental é a que reconhece que "cloridrato de X" e "X"
           são o mesmo composto, e é a coluna correta para cruzar as duas
           bases. Divergência entre o InChIKey declarado e o recalculado é
           reportada na coluna `conferencia`, sem descartar a linha.

USO
---
    python mods/importar_prospector.py \
        --origem ~/prospector/pipeline_prospeccao_v4/resultados/atividades_compostos_parcial.csv \
        --saida input.txt

    # amostra para teste rápido, priorizando compostos com evidência ChEMBL
    python mods/importar_prospector.py --origem ... --saida input.exemplo.tsv \
        --limite 200 --preferir-com-chembl

REQUISITOS
----------
    pip install rdkit
"""

import argparse
import csv
import sys
from pathlib import Path

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem.MolStandardize import rdMolStandardize
    RDLogger.DisableLog("rdApp.*")
except ImportError:
    print("[ERRO] Este script precisa do RDKit: pip install rdkit", file=sys.stderr)
    sys.exit(1)

# Ordem de prioridade das fontes de SMILES, igual à da prospecção.
ORDEM_FONTES_SMILES = [
    ("NPASS", "smiles_npass"),
    ("COCONUT", "smiles_coconut"),
    ("LOTUS", "smiles_lotus"),
    ("PubChem", "smiles_pubchem"),
]

COLUNAS_SAIDA = [
    "nome", "smiles", "inchikey", "inchikey_parental",
    "fonte_smiles", "conferencia", "chembl_id",
]

_PADRONIZADORES: dict = {}


def _padronizadores() -> dict:
    if not _PADRONIZADORES:
        _PADRONIZADORES["frag"] = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
        _PADRONIZADORES["uncharger"] = rdMolStandardize.Uncharger()
    return _PADRONIZADORES


def texto(valor) -> str:
    """Normaliza um campo de CSV para string limpa ('' quando ausente)."""
    if valor is None:
        return ""
    s = str(valor).strip()
    return "" if s.lower() in ("", "nan", "none", "-", "null") else s


def molecula(smiles: str):
    if not smiles:
        return None
    try:
        return Chem.MolFromSmiles(smiles)
    except Exception:
        return None


def forma_parental(mol):
    """Cleanup -> maior fragmento orgânico -> neutralização. Mesmo tratamento
    que a prospecção aplica para reconhecer sais e formas carregadas."""
    if mol is None:
        return None
    try:
        p = _padronizadores()
        m = rdMolStandardize.Cleanup(mol)
        m = p["frag"].choose(m)
        return p["uncharger"].uncharge(m)
    except Exception:
        return None


def inchikey_de(mol) -> str:
    if mol is None:
        return ""
    try:
        return Chem.MolToInchiKey(mol) or ""
    except Exception:
        return ""


def escolher_smiles(linha: dict) -> tuple[str, str]:
    """Primeira fonte cujo SMILES o RDKit consegue ler, na ordem da
    prospecção. Devolve (smiles, fonte). Se nenhuma candidata for válida,
    cai para a coluna `smiles` já escolhida pela prospecção."""
    for rotulo, coluna in ORDEM_FONTES_SMILES:
        smi = texto(linha.get(coluna))
        if smi and molecula(smi) is not None:
            return smi, rotulo
    smi = texto(linha.get("smiles"))
    if smi and molecula(smi) is not None:
        return smi, texto(linha.get("fonte_smiles")) or "prospeccao"
    return "", ""


def converter_linha(linha: dict) -> dict | None:
    """Uma linha da prospecção -> uma linha do formato de entrada, ou None
    se o composto não tiver nome ou SMILES utilizável."""
    nome = texto(linha.get("composto_nome")) or texto(linha.get("nome_pubchem"))
    smiles, fonte = escolher_smiles(linha)
    if not nome or not smiles:
        return None

    mol = molecula(smiles)
    pai = forma_parental(mol)
    ik = inchikey_de(mol)
    ik_pai = inchikey_de(pai)

    # O veredito da prospecção é preservado; acrescentamos a conferência da
    # identidade declarada contra a recalculada.
    notas = []
    veredito = texto(linha.get("conferencia_smiles"))
    if veredito:
        notas.append(veredito)
    declarado = texto(linha.get("InChIKey"))
    if declarado and ik and declarado != ik:
        if declarado.split("-")[0] == ik.split("-")[0]:
            notas.append(f"InChIKey declarado difere na estereoquímica/carga: {declarado}")
        else:
            notas.append(f"InChIKey declarado é outro esqueleto: {declarado}")

    return {
        "nome": nome,
        # Grava o SMILES canônico: a prospecção compara estruturas, não texto.
        "smiles": Chem.MolToSmiles(mol),
        "inchikey": ik,
        "inchikey_parental": ik_pai,
        "fonte_smiles": fonte,
        "conferencia": "; ".join(notas),
        "chembl_id": texto(linha.get("chembl_id")),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Converte a tabela de compostos da prospecção de plantas no "
                    "formato de entrada deste pipeline (TSV).")
    parser.add_argument("--origem", required=True,
                        help="atividades_compostos_parcial.csv da prospecção.")
    parser.add_argument("--saida", default="input.txt", help="Arquivo TSV de saída.")
    parser.add_argument("--limite", type=int, default=None,
                        help="Grava no máximo N compostos (amostra para teste).")
    parser.add_argument("--preferir-com-chembl", action="store_true",
                        help="Com --limite, preenche metade da amostra com compostos que "
                             "têm chembl_id (úteis para validação) e metade sem.")
    parser.add_argument("--manter-divergentes", action="store_true",
                        help="Mantém compostos cujo veredito aponta estruturas diferentes "
                             "entre as bases (o padrão já os mantém, marcados).")
    args = parser.parse_args()

    origem = Path(args.origem).expanduser()
    if not origem.is_file():
        print(f"[ERRO] Arquivo não encontrado: {origem}", file=sys.stderr)
        sys.exit(1)

    csv.field_size_limit(10 ** 9)

    vistos: set[str] = set()
    com_chembl: list[dict] = []
    sem_chembl: list[dict] = []
    lidas = invalidas = duplicadas = 0

    # Com --limite, basta encher as cestas: não há motivo para padronizar as
    # 100 mil moléculas da base inteira só para gravar uma amostra.
    if args.limite:
        if args.preferir_com_chembl:
            alvo_com = args.limite // 2
            alvo_sem = args.limite - alvo_com
        else:
            alvo_com = alvo_sem = args.limite
    else:
        alvo_com = alvo_sem = None

    def cestas_cheias() -> bool:
        return (alvo_com is not None
                and len(com_chembl) >= alvo_com and len(sem_chembl) >= alvo_sem)

    with open(origem, "r", encoding="utf-8-sig", newline="") as f:
        for linha in csv.DictReader(f):
            lidas += 1
            if lidas % 20000 == 0:
                print(f"[INFO]   {lidas:,} linhas lidas, {len(vistos):,} compostos únicos...",
                      flush=True)
            conv = converter_linha(linha)
            if conv is None:
                invalidas += 1
                continue
            # Deduplica pela forma parental: sal e base livre são o mesmo composto.
            chave = conv["inchikey_parental"] or conv["inchikey"] or conv["smiles"]
            if chave in vistos:
                duplicadas += 1
                continue
            vistos.add(chave)
            cesta = com_chembl if conv["chembl_id"] else sem_chembl
            limite_cesta = alvo_com if conv["chembl_id"] else alvo_sem
            if limite_cesta is None or len(cesta) < limite_cesta:
                cesta.append(conv)
            if cestas_cheias():
                print(f"[INFO]   amostra completa em {lidas:,} linhas lidas.", flush=True)
                break

    if args.limite and args.preferir_com_chembl:
        metade = args.limite // 2
        selecionados = com_chembl[:metade] + sem_chembl[:args.limite - metade]
    else:
        selecionados = com_chembl + sem_chembl
        if args.limite:
            selecionados = selecionados[:args.limite]

    saida = Path(args.saida)
    with open(saida, "w", encoding="utf-8", newline="") as f:
        # O cabeçalho começa com '#': é legível como comentário e o leitor
        # deste pipeline o usa para identificar as colunas.
        f.write("#" + "\t".join(COLUNAS_SAIDA) + "\n")
        escritor = csv.DictWriter(f, fieldnames=COLUNAS_SAIDA, delimiter="\t",
                                  lineterminator="\n", extrasaction="ignore")
        for linha in selecionados:
            escritor.writerow(linha)

    divergentes = sum(1 for l in selecionados if "estruturas diferentes" in l["conferencia"])
    parcial = " (leitura interrompida ao completar a amostra)" if cestas_cheias() else ""
    print(f"[INFO] Linhas lidas da prospecção      : {lidas:,}{parcial}")
    print(f"[INFO]   sem nome ou SMILES utilizável : {invalidas:,}")
    print(f"[INFO]   duplicadas (mesma parental)   : {duplicadas:,}")
    print(f"[INFO] Compostos únicos vistos         : {len(vistos):,}")
    print(f"[INFO] Gravados em '{saida}'           : {len(selecionados):,} "
          f"({sum(1 for l in selecionados if l['chembl_id'])} com chembl_id)")
    if divergentes:
        print(f"[AVISO] {divergentes} composto(s) com estruturas divergentes entre as bases "
              f"(marcados na coluna 'conferencia' -- revise antes de confiar no resultado).")


if __name__ == "__main__":
    main()
