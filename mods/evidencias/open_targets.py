# -*- coding: utf-8 -*-
"""
open_targets.py

Do sintoma aos alvos: o Open Targets associa doencas (EFO) a alvos com um
escore; o ChEMBL diz qual alvo do modelo corresponde a cada UniProt.

Regras
------
* Sintoma -> EFO nunca e escolhido aqui: `buscar_efo` devolve candidatos e o
  chamador decide. Um EFO errado produz uma lista de alvos errada sem sinal.
* O escore se espalha pela familia do gene (22 subunidades de canal de calcio
  com 0,603-0,605 no topo da dor neuropatica). `colapsar_familia` so dobra
  para dentro de um representante os paralogos AUSENTES do espaco do modelo.
  Todo paralogo presente no modelo continua sendo uma linha propria: a classe
  "SLC06 neurotransmitter transporter" contem NET e SERT, e fundi-los
  apagaria o mecanismo duplo que explica a duloxetina. O que foi dobrado fica
  registrado em `familia_colapsados`.
* Enriquecimento so acrescenta: nenhum alvo do usuario sai, muda de ordem ou
  e sobrescrito.
* Falha por item (`ListaTruncada`, `RespostaInesperada`) vira `observacao` na
  linha; o item nao some e os demais seguem.
* `obter_json(url, corpo)` e injetavel; `corpo` e a consulta GraphQL (POST).
"""

from typing import Callable

from mods.evidencias.atomico import escrever_texto
from mods.evidencias.chembl import (
    ListaTruncada, RespostaInesperada, alvos_por_uniprot, erro_de_item)
from mods.evidencias.esquema import escapar

URL = "https://api.platform.opentargets.org/api/v4/graphql"
ORGANISMO = "Homo sapiens"
PAGINA = 100

COLUNAS_ENRIQUECIDOS = [
    "alvo_chembl_id", "alvo", "organismo", "tipo_alvo", "gene", "ensembl",
    "origem", "efo_id", "sintoma", "escore_associacao", "escore_rna_expression",
    "familia", "familia_colapsados", "via_origem", "alvo_semente",
    "n_moleculas_treino", "no_modelo", "observacao",
]
# `ensembl` e `escore_rna_expression` ja vinham na resposta do Open Targets e
# eram descartados na escrita. Ficam na saida porque a etapa 2 precisa dos
# dois: o Ensembl e a chave da expressao por tecido (`baselineExpression`), e
# `rna_expression` e o escore de expressao DIFERENCIAL da associacao -- a
# contribuicao do Expression Atlas, que indexa os experimentos do GEO. Buscar
# o Ensembl de novo na etapa 2 custaria uma requisicao por alvo para um dado
# que esta etapa ja tinha em maos.

_BUSCA = """query Busca($q: String!) {
  search(queryString: $q, entityNames: ["disease"], page: {index: 0, size: 10}) {
    hits { id name }
  }
}"""

_ASSOCIADOS = """query Associados($efo: String!, $i: Int!, $n: Int!) {
  disease(efoId: $efo) {
    associatedTargets(page: {index: $i, size: $n}) {
      count
      rows {
        score
        datatypeScores { id score }
        target { id approvedSymbol proteinIds { id source } }
      }
    }
  }
}"""


def _post_json_requests(url: str, corpo: dict) -> dict:
    import requests
    resposta = requests.post(url, json=corpo, timeout=60)
    resposta.raise_for_status()
    return resposta.json()


def _dados(resposta: dict, consulta: str) -> dict:
    if resposta.get("errors") or not isinstance(resposta.get("data"), dict):
        raise RespostaInesperada(f"{consulta}: resposta GraphQL com erro ou sem 'data'")
    return resposta["data"]


def buscar_efo(sintoma: str, *,
               obter_json: Callable = _post_json_requests) -> list[dict]:
    """Candidatos de EFO (`id`, `name`) para o sintoma, melhor primeiro. Nunca
    escolhe: quem chama mostra e decide (ou fixa com --efo)."""
    corpo = {"query": _BUSCA, "variables": {"q": sintoma}}
    dados = _dados(obter_json(URL, corpo), f"busca de {sintoma!r}")
    try:
        acertos = dados["search"]["hits"]
    except (KeyError, TypeError):
        raise RespostaInesperada(f"busca de {sintoma!r}: resposta sem search.hits")
    return [{"id": h["id"], "name": h["name"]} for h in acertos]


def alvos_associados(efo_id: str, *, corte: float,
                     obter_json: Callable = _post_json_requests) -> list[dict]:
    """Alvos associados ao EFO com `score >= corte`, do maior para o menor.
    Cada um traz `gene`, `ensembl`, `score`, `datatypes` ({tipo: escore}) e
    `uniprot` (o Swiss-Prot; vazio se o Open Targets nao tem)."""
    achados: list[dict] = []
    indice = 0
    while True:
        corpo = {"query": _ASSOCIADOS,
                 "variables": {"efo": efo_id, "i": indice, "n": PAGINA}}
        dados = _dados(obter_json(URL, corpo), f"alvos de {efo_id}")
        try:
            bloco = dados["disease"]["associatedTargets"]
            linhas, total = bloco["rows"], bloco["count"]
        except (KeyError, TypeError):
            raise RespostaInesperada(f"alvos de {efo_id}: EFO desconhecido ou "
                                     "resposta sem associatedTargets")
        for linha in linhas:
            if linha["score"] < corte:
                return achados  # a lista vem em ordem decrescente de escore
            alvo = linha["target"]
            uniprot = next((p["id"] for p in alvo.get("proteinIds") or []
                            if p.get("source") == "uniprot_swissprot"), "")
            achados.append({
                "gene": alvo["approvedSymbol"], "ensembl": alvo["id"],
                "score": linha["score"], "uniprot": uniprot, "efo_id": efo_id,
                "datatypes": {d["id"]: d["score"]
                              for d in linha.get("datatypeScores") or []}})
        indice += 1
        if not linhas or indice * PAGINA >= total:
            return achados


def alvos_do_chembl(candidatos: list[dict], rotulos: list[dict], *,
                    obter_json_chembl: Callable[[str], dict],
                    organismo: str = ORGANISMO) -> list[dict]:
    """Casa cada candidato do Open Targets com um alvo do ChEMBL (via UniProt)
    e marca se esta no espaco do modelo. Candidato que nao casa fica na saida,
    com `observacao`; nada some."""
    no_modelo = {r["chembl_id"]: r for r in rotulos}
    saida = []
    for c in candidatos:
        linha = {**c, "origem": "open_targets", "alvo_chembl_id": "",
                 "alvo": "", "organismo": organismo, "tipo_alvo": "",
                 "no_modelo": False, "n_moleculas_treino": "", "observacao": "",
                 "escore_rna_expression": (c.get("datatypes") or {}).get(
                     "rna_expression", "")}
        saida.append(linha)
        if not c.get("uniprot"):
            linha["observacao"] = "sem UniProt Swiss-Prot no Open Targets"
            continue
        try:
            alvos = alvos_por_uniprot(c["uniprot"], organismo,
                                      obter_json=obter_json_chembl)
        except (ListaTruncada, OSError, ValueError) as e:
            if not erro_de_item(e):
                raise
            linha["observacao"] = f"ChEMBL: {e}"
            continue
        if not alvos:
            linha["observacao"] = f"{c['uniprot']} sem alvo no ChEMBL"
            continue
        alvo = alvos[0]
        linha.update(alvo_chembl_id=alvo["target_chembl_id"],
                     alvo=alvo.get("pref_name") or "",
                     tipo_alvo=alvo.get("target_type") or "")
        if alvo.get("target_type") != "SINGLE PROTEIN":
            linha["observacao"] = (f"{c['uniprot']} sem alvo SINGLE PROTEIN; "
                                   f"usado {alvo.get('target_type')}")
        modelo = no_modelo.get(linha["alvo_chembl_id"])
        if modelo:
            linha["no_modelo"] = True
            linha["n_moleculas_treino"] = modelo["n_ativos"]
    return saida


def colapsar_familia(alvos: list[dict],
                     relatorio: dict | None = None) -> list[dict]:
    """Dobra, em cada classe de proteina do ChEMBL, os alvos AUSENTES do modelo
    (`no_modelo` False) no de maior escore que fica. Quem esta no modelo, quem
    e do usuario e quem nao tem classe nunca e dobrado. O dobrado sai da lista;
    seu ID vai para `familia_colapsados` do representante e a classe para
    `familia`. A ordem de entrada e preservada."""
    por_classe: dict[str, list[dict]] = {}
    for a in alvos:
        if a.get("classe"):
            por_classe.setdefault(a["classe"], []).append(a)
    rel = relatorio if relatorio is not None else {}
    rel.setdefault("dobrados_por_familia", 0)
    rel.setdefault("dobrados_fora_do_modelo", [])
    dobrados: set[int] = set()
    for classe, membros in por_classe.items():
        ordem = sorted(membros, key=lambda m: -(m.get("score") or 0))
        fica = [m for m in ordem
                if m.get("no_modelo") is not False or m.get("origem") == "usuario"]
        representante = (fica or ordem)[0]
        presentes = {id(m) for m in fica}
        ausentes = [m for m in ordem
                    if m is not representante and id(m) not in presentes]
        if not ausentes:
            continue
        representante["familia"] = classe
        representante["familia_colapsados"] = ";".join(
            m["alvo_chembl_id"] for m in ausentes)
        dobrados.update(id(m) for m in ausentes)
        # O que dobrou tem de ficar CONTADO: ele sai da lista, e sem a
        # contagem desapareceria do teto do modelo -- o numero que justifica
        # retreinar o preditor. Os IDs dos que estavam fora do modelo saem
        # separados porque sao a classe que a evidencia medida ainda alcanca.
        rel["dobrados_por_familia"] += len(ausentes)
        rel["dobrados_fora_do_modelo"].extend(
            m["alvo_chembl_id"] for m in ausentes
            if m.get("no_modelo") is False)
    return [a for a in alvos if id(a) not in dobrados]


def texto_dos_alvos_enriquecidos(alvos: list[dict]) -> str:
    """O arquivo inteiro como texto. A lista de alvos tem centenas de linhas
    (nao milhoes, como a de evidencias), entao monta-la em memoria e o que
    permite prepara-la antes de publicar qualquer saida."""
    linhas = ["#" + "\t".join(COLUNAS_ENRIQUECIDOS)]
    for a in alvos:
        celulas = []
        for c in COLUNAS_ENRIQUECIDOS:
            v = a.get("score") if c == "escore_associacao" and "score" in a else a.get(c)
            if isinstance(v, bool):
                v = "sim" if v else "nao"
            celulas.append(escapar(v))
        linhas.append("\t".join(celulas))
    return "\n".join(linhas) + "\n"


def escrever_alvos_enriquecidos(caminho, alvos: list[dict]) -> None:
    """Segunda saida da etapa. Mesmo formato do `evidencias.tsv`: TSV, cabecalho
    '#', escape reversivel, ponto decimal. Grava tudo que recebe.

    Escrita atomica: este e o arquivo que o README manda revisar a mao, e uma
    falha no meio da escrita nao pode destruir a versao revisada."""
    escrever_texto(caminho, texto_dos_alvos_enriquecidos(alvos))
