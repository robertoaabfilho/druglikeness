# -*- coding: utf-8 -*-
"""
expressao.py

Expressao do ALVO por tecido, para responder "este alvo existe no tecido da
doenca?". Fonte: `target.baselineExpression` do Open Targets (GTEx, RNA-seq de
tecido inteiro, TPM; e Tabula Sapiens, scRNA-seq, CPM). CC0, sem download.

Por que isto MARCA e nunca reprova
----------------------------------
Medido nos 9 alvos curados de dor neuropatica, maior TPM mediano do GTEx em
tecido nervoso:

    GRIN1 255,3 | CACNA2D2 39,2 | TRPV1 22,7 | CACNA2D1 16,7 | SCN1A 12,0
    OPRM1 5,9 | SLC6A4 (SERT) 0,45 | SLC6A2 (NET) 0,08

Um portao de ">= 1 TPM em tecido nervoso" reprovaria SERT e NET: os dois
alvos da duloxetina e da amitriptilina, os mais bem validados da lista. A
causa e conhecida -- RNA-seq de tecido inteiro dilui nucleos pequenos (locus
coeruleus, nucleos da rafe) no cortex amostrado. Expressao baixa aqui nao e
evidencia de ausencia.

GEO
---
O GEO cru nao entra. A via mantida e harmonizada para esse mesmo dado e o
Expression Atlas, que o Open Targets indexa como fonte `expression_atlas` e
resume no escore de tipo de dado `rna_expression` (coluna
`escore_rna_expression`, gravada pela etapa 1). Medido: para `EFO_0005762`
(dor neuropatica) essa fonte devolve ZERO evidencias em SLC6A2, TRPV1 e SCN9A
-- nao existe transcriptomica humana de caso-controle para a doenca, porque
nao se biopsia ganglio da raiz dorsal. A coluna sai vazia e isso fica dito, em
vez de a ausencia parecer um numero.

Casamento de tecido
-------------------
Por identificador UBERON e por ANCESTRALIDADE, nunca por trecho de nome: numa
medicao preliminar, um padrao de nome para "nervoso" casou `cortex of kidney`.
`--tecidos` aceita `UBERON_0000955` ou o nome exato do tecido; um ancestral
(`UBERON_0000955`, cerebro) alcanca `frontal cortex` porque o proprio Open
Targets devolve a lista de ancestrais de cada amostra.
"""

import re
from typing import Callable

from mods.evidencias.open_targets import URL, _dados, _post_json_requests

PAGINA = 2000
# GTEx e a fonte comparavel entre tecidos (TPM de tecido inteiro); o Tabula
# Sapiens usa CPM de pseudo-bulk e nao se compara com ela na mesma coluna.
FONTE_PADRAO = "gtex"

_EXPRESSAO = """query Expressao($e: String!, $i: Int!, $n: Int!) {
  target(ensemblId: $e) {
    approvedSymbol
    baselineExpression(page: {index: $i, size: $n}) {
      count
      rows {
        datasourceId unit median max
        tissueBiosample { biosampleId biosampleName ancestors }
      }
    }
  }
}"""

_BUSCA_ALVO = """query BuscaAlvo($q: String!) {
  search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 5}) {
    hits { id name }
  }
}"""

_RE_UBERON = re.compile(r"(?:UBERON|CL|EFO)[_:]\d+", re.IGNORECASE)


def normalizar_tecido(termo: str) -> str:
    """`uberon:0000955` e `UBERON_0000955` sao o mesmo tecido."""
    termo = termo.strip()
    if _RE_UBERON.fullmatch(termo):
        return termo.replace(":", "_").upper()
    return termo.casefold()


def ler_tecidos(caminho) -> list[str]:
    """Um tecido por linha, '#' comenta. Aceita ID UBERON ou nome exato."""
    pedidos = []
    with open(caminho, encoding="utf-8-sig") as f:
        for linha in f:
            linha = linha.split("#", 1)[0].strip()
            if linha:
                pedidos.append(normalizar_tecido(linha))
    return list(dict.fromkeys(pedidos))


def ensembl_do_gene(gene: str, *, obter_json: Callable = _post_json_requests) -> str:
    """Ensembl de um simbolo de gene. So e usado quando a lista de alvos vem de
    um `alvos_enriquecidos.tsv` antigo, sem a coluna `ensembl`."""
    corpo = {"query": _BUSCA_ALVO, "variables": {"q": gene}}
    dados = _dados(obter_json(URL, corpo), f"busca de alvo {gene!r}")
    for h in (dados.get("search") or {}).get("hits") or []:
        if str(h.get("id", "")).startswith("ENSG"):
            return h["id"]
    return ""


def expressao_basal(ensembl: str, *, obter_json: Callable = _post_json_requests,
                    fonte: str = FONTE_PADRAO) -> list[dict]:
    """Linhas de expressao basal da fonte pedida, uma por tecido.

    Uma pagina de 2000 cobre o maior caso visto (1.428 linhas para SLC6A2);
    quando nao cobrir, pagina. Deduplicada por tecido pela maior mediana: o
    Tabula Sapiens repete o mesmo tecido por tipo celular."""
    linhas: list[dict] = []
    indice = 0
    while True:
        corpo = {"query": _EXPRESSAO,
                 "variables": {"e": ensembl, "i": indice, "n": PAGINA}}
        dados = _dados(obter_json(URL, corpo), f"expressao de {ensembl}")
        alvo = dados.get("target")
        if alvo is None:
            return []
        bloco = alvo["baselineExpression"]
        linhas.extend(bloco["rows"])
        indice += 1
        if not bloco["rows"] or indice * PAGINA >= bloco["count"]:
            break
    por_tecido: dict[str, dict] = {}
    for linha in linhas:
        if fonte and linha.get("datasourceId") != fonte:
            continue
        amostra = linha.get("tissueBiosample") or {}
        chave = amostra.get("biosampleId") or amostra.get("biosampleName") or ""
        if not chave:
            continue
        registro = {
            "tecido_id": chave,
            "tecido": amostra.get("biosampleName") or "",
            "ancestrais": [a.replace(":", "_").upper()
                           for a in (amostra.get("ancestors") or [])],
            "mediana": linha.get("median"),
            "maximo": linha.get("max"),
            "unidade": linha.get("unit") or "",
            "fonte": linha.get("datasourceId") or "",
        }
        anterior = por_tecido.get(chave)
        if anterior is None or (registro["mediana"] or 0) > (anterior["mediana"] or 0):
            por_tecido[chave] = registro
    return sorted(por_tecido.values(), key=lambda r: -(r["mediana"] or 0))


def casa_tecido(linha: dict, pedidos: list[str]) -> bool:
    """O tecido da linha e um dos pedidos, ou descende de um deles.

    Por ID e por ancestralidade. Nome e comparado inteiro, em minusculas --
    nunca por trecho: `cortex of kidney` contem `cortex`."""
    if not pedidos:
        return True
    identificadores = {linha["tecido_id"].replace(":", "_").upper(),
                       *linha["ancestrais"]}
    nome = (linha["tecido"] or "").casefold()
    return any(p in identificadores or p == nome for p in pedidos)


def anotar(linhas: list[dict], pedidos: list[str], *, corte: float):
    """Colunas de expressao de um alvo, a partir das linhas da fonte.

    `corte` nao reprova nada por si: o veredito sai `marca` quando a expressao
    nos tecidos pedidos fica abaixo dele, e a politica decide o resto. Sem
    linha nenhuma o veredito e `indeterminado` -- alvo nao humano nao tem GTEx,
    e isso nao e expressao zero."""
    if not linhas:
        return {"expressao_maxima": "", "tecido_maximo": "",
                "expressao_nos_tecidos_pedidos": "",
                "tecidos_pedidos_casados": "", "unidade_expressao": "",
                "fonte_expressao": "", "expressao_veredito": "indeterminado",
                "expressao_observacao": "sem expressao basal nesta fonte"}
    melhor = linhas[0]
    casados = [l for l in linhas if casa_tecido(l, pedidos)]
    no_pedido = max((l["mediana"] or 0) for l in casados) if casados else None
    if not pedidos:
        veredito, observacao = "indeterminado", "nenhum tecido pedido (--tecidos)"
    elif not casados:
        veredito, observacao = ("indeterminado",
                                "nenhum tecido pedido existe nesta fonte")
    elif no_pedido >= corte:
        veredito, observacao = "passa", ""
    else:
        veredito = "marca"
        observacao = (f"{no_pedido:g} < {corte:g} nos tecidos pedidos; RNA-seq "
                      "de tecido inteiro dilui nucleo pequeno (SERT e NET "
                      "ficam abaixo de 1 TPM e sao alvos validados)")
    return {
        "expressao_maxima": melhor["mediana"],
        "tecido_maximo": melhor["tecido"],
        "expressao_nos_tecidos_pedidos": "" if no_pedido is None else no_pedido,
        "tecidos_pedidos_casados": ";".join(sorted(l["tecido"] for l in casados)),
        "unidade_expressao": melhor["unidade"],
        "fonte_expressao": melhor["fonte"],
        "expressao_veredito": veredito,
        "expressao_observacao": observacao,
    }
