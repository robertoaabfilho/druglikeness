# -*- coding: utf-8 -*-
"""
propriedades.py

Propriedades que o RDKit calcula localmente: descritores de drogabilidade, QED
e alerta estrutural NOMEADO.

Por que daqui e nao do ADMET-AI
-------------------------------
O ADMET-AI tambem devolve `PAINS_alert`, `BRENK_alert` e `NIH_alert`, mas como
CONTAGEM, sem dizer qual filtro casou. "BRENK = 1" nao se confere; "BRENK:
nitro_group" se confere. As duas colunas convivem: a do ADMET-AI sai com o
nome original dele (e com percentil contra os aprovados do DrugBank), e a daqui
sai nomeada, com `procedencia` dizendo RDKit. Duas colunas respondendo a mesma
pergunta sao aceitaveis quando a procedencia de cada uma esta escrita; uma
coluna sem procedencia nao e.

E, na pratica, esta e a unica fonte de alerta estrutural quando o `admet_ai`
nao esta instalado -- o repositorio roda sem ele (modo degradado).

Medido (ver o desenho da etapa 2)
---------------------------------
Em `input.exemplo.tsv` (200 produtos naturais): PAINS 15%, BRENK 65%, NIH 30%.
No conjunto-controle: o PAINS pega a AMITRIPTILINA, que trata dor neuropatica.
Por isso alerta estrutural MARCA; quem reprova por ele apaga o proprio
positivo de referencia do pipeline.
"""

from mods.importar_prospector import molecula

# Catalogos do RDKit, pelo nome que o proprio RDKit usa. PAINS entra nas tres
# severidades (A, B, C) porque o artigo original as define como um conjunto.
CATALOGOS = {
    "alerta_pains": ("PAINS_A", "PAINS_B", "PAINS_C"),
    "alerta_brenk": ("BRENK",),
    "alerta_nih": ("NIH",),
}

_catalogos_prontos: dict | None = None


def _carregar_catalogos():
    """Monta os catalogos uma vez por processo: cada `FilterCatalog` compila
    centenas de SMARTS, e refaze-lo por molecula dominaria o custo."""
    global _catalogos_prontos
    if _catalogos_prontos is None:
        from rdkit.Chem import FilterCatalog
        prontos = {}
        for coluna, nomes in CATALOGOS.items():
            params = FilterCatalog.FilterCatalogParams()
            for nome in nomes:
                params.AddCatalog(getattr(
                    FilterCatalog.FilterCatalogParams.FilterCatalogs, nome))
            prontos[coluna] = FilterCatalog.FilterCatalog(params)
        _catalogos_prontos = prontos
    return _catalogos_prontos


def alertas_estruturais(mol) -> dict[str, str]:
    """{coluna: nomes dos filtros casados, separados por ';'}. Vazio quando o
    catalogo nao casa nada -- distinto de None, que nao se sabe."""
    saida = {}
    for coluna, catalogo in _carregar_catalogos().items():
        nomes = sorted({m.GetDescription() for m in catalogo.GetMatches(mol)})
        saida[coluna] = ";".join(nomes)
    return saida


def descritores(mol) -> dict:
    """Os mesmos descritores de `farmaco_completo.calcular_descritores`, mais o
    QED. Nao importa aquele modulo: ele carrega a API do Deep-PK, que esta
    etapa substitui, e importa `configuracao` do diretorio de scripts."""
    from rdkit.Chem import Descriptors, QED
    return {
        "peso_molecular": round(Descriptors.MolWt(mol), 2),
        "doadores_hb": Descriptors.NumHDonors(mol),
        "aceptores_hb": Descriptors.NumHAcceptors(mol),
        "tpsa": round(Descriptors.TPSA(mol), 2),
        "ligacoes_rotacionaveis": Descriptors.NumRotatableBonds(mol),
        "refratividade_molar": round(Descriptors.MolMR(mol), 2),
        "atomos_pesados": Descriptors.HeavyAtomCount(mol),
        # QED e coluna, nunca criterio: no conjunto-controle os dois QED mais
        # altos (0,89) sao ciprofloxacino e linezolida, do grupo `causa`.
        "qed": round(QED.qed(mol), 4),
    }


def versao_rdkit() -> str:
    import rdkit
    return getattr(rdkit, "__version__", "desconhecida")


def linha_de_propriedades(ligante: dict) -> tuple[dict, str | None]:
    """(linha de `propriedades.tsv` sem as colunas do ADMET-AI, problema).

    SMILES que o RDKit nao le devolve a linha de identidade com as colunas de
    propriedade VAZIAS e o problema descrito: o ligante nao desaparece da
    saida so porque nao pode ser avaliado."""
    identidade = {
        "ligante": ligante.get("nome", ""),
        "ligante_inchikey_parental": ligante.get("inchikey_parental", ""),
        "ligante_chembl_id": ligante.get("chembl_id", ""),
        "smiles": ligante.get("smiles", ""),
    }
    mol = molecula(ligante.get("smiles", ""))
    if mol is None:
        return identidade, "SMILES que o RDKit nao le"
    linha = {**identidade, **descritores(mol), **alertas_estruturais(mol),
             "procedencia": f"RDKit {versao_rdkit()}"}
    return linha, None
