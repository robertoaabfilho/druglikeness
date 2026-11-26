# -*- coding: utf-8 -*-
"""
bula.py

Evidencia NEGATIVA: a bula (openFDA) lista a propria condicao-alvo como efeito
adverso. O ciprofloxacino tem tarja preta da FDA para neuropatia periferica e
aparecia no topo de uma busca por analgesia -- este modulo e a resposta a isso.

A regra anterior errava nos dois sentidos
-----------------------------------------
`mods/checar_alertas_adversos.py` busca com `limit=1` em quatro secoes, sem
descontar a indicacao. Medido no conjunto-controle de 12 farmacos:

* 4 de 6 FALSOS POSITIVOS no grupo que TRATA (pregabalina, duloxetina,
  amitriptilina, carbamazepina): a bula de um farmaco para dor neuropatica
  lista parestesia e neuropatia entre as reacoes adversas;
* FALSO NEGATIVO no caso que originou o script: o ciprofloxacino NAO era
  sinalizado, porque `limit=1` devolve uma bula arbitraria entre centenas e a
  que traz a tarja preta nao e a primeira.

Tres correcoes, e os numeros depois delas
-----------------------------------------
1. `limit=50` bulas por principio ativo, AGREGADAS: qualquer bula com o termo
   conta. (Foi isso que trouxe a tarja preta do ciprofloxacino.)
2. Secao FORTE (`boxed_warning`, `warnings`, `warnings_and_cautions`) conta
   como alerta; `adverse_reactions` e registrada e NAO conta.
3. INDICACAO desconta: termo em `indications_and_usage` significa que o
   farmaco TRATA a condicao, e cancela o alerta.

Resultado no mesmo controle: 0 de 6 no grupo `trata` (pregabalina e duloxetina
batem em secao forte, mas sao indicadas) e 4 de 6 no grupo `causa`
(ciprofloxacino, linezolida, etambutol, metronidazol). Os dois nao pegos sao
honestos: cloranfenicol tem UMA bula no openFDA e etionamida nenhuma -- ambos
saem `indeterminado`, nunca "sem alerta".

O que isto NAO faz
------------------
Nao entende negacao: uma bula que diga "no cases of neuropathy were reported"
ainda casa o termo. Por isso o trecho vai inteiro para `alertas_bula.tsv`, com
a secao e o numero de bulas, para conferencia humana. E por isso tambem a
busca por termo e por PALAVRA inteira: sem isso "neuritis" casaria dentro de
outras palavras e o falso positivo nao seria nem conferivel.
"""

import re
from typing import Callable
from urllib.parse import quote

BASE = "https://api.fda.gov/drug/label.json"

# Secoes da bula onde o termo tem significado de RISCO, por forca.
SECOES_FORTES = ("boxed_warning", "warnings", "warnings_and_cautions")
SECOES_FRACAS = ("adverse_reactions",)
# Nao e secao de risco: e o contrario. Entra para DESCONTAR.
SECOES_INDICACAO = ("indications_and_usage",)

# Ordem de busca do nome: generico, principio ativo, nome comercial.
CAMPOS_NOME = ("generic_name", "substance_name", "brand_name")

LIMITE_BULAS = 50
JANELA_TRECHO = 120


SECOES_POR_FORCA = {"forte": SECOES_FORTES, "fraca": SECOES_FRACAS,
                    "indicacao": SECOES_INDICACAO}

# Ordenacao da janela. Medido: a capsaicina tem 211 bulas e, sem ordenar, as 50
# que a busca devolve sao arbitrarias -- a unica que a indica para dor
# neuropatica ficava de fora. Por data decrescente ela entra (effective_time
# 2026-07-07). O openFDA NAO publica acesso, popularidade nem contagem de uso
# (a bula tem 34 campos e nenhum deles), entao "as 50 mais acessadas" nao
# existe; a data e o unico critetio de relevancia disponivel, e significa
# "rotulagem vigente".
ORDEM_DA_JANELA = "effective_time:desc"


def url(nome: str, campo: str, limite: int = LIMITE_BULAS,
        ordem: str = ORDEM_DA_JANELA) -> str:
    busca = quote(f'openfda.{campo}:"{nome}"', safe="")
    sufixo = f"&sort={quote(ordem, safe=':')}" if ordem else ""
    return f"{BASE}?search={busca}{sufixo}&limit={int(limite)}"


def versao_openfda(*, obter_json: Callable) -> str:
    """`meta.last_updated` do proprio endpoint: o openFDA nao tem numero de
    versao, e a data da ultima carga e o que identifica o conteudo."""
    resposta = obter_json(f"{BASE}?limit=1")
    valor = ((resposta or {}).get("meta") or {}).get("last_updated")
    if not valor:
        raise ValueError("openFDA: resposta sem meta.last_updated")
    return str(valor)


def buscar_bulas(nome: str, *, obter_json: Callable,
                 limite: int = LIMITE_BULAS) -> dict:
    """{'campo_do_nome', 'bulas', 'total'}: a primeira busca com resultado.

    `bulas` vazio significa que nenhum dos tres campos de nome achou algo --
    nao que o farmaco seja seguro. `total` e quantas bulas existem para o
    termo, para que `n_bulas_consultadas` possa ser comparado com ele."""
    for campo in CAMPOS_NOME:
        resposta = obter_json(url(nome, campo, limite)) or {}
        bulas = resposta.get("results") or []
        if bulas:
            total = (((resposta.get("meta") or {}).get("results") or {})
                     .get("total", len(bulas)))
            return {"campo_do_nome": campo, "bulas": bulas, "total": total}
    return {"campo_do_nome": "", "bulas": [], "total": 0}


def url_dirigida(nome: str, campo: str, termos: list[str],
                 secoes: tuple[str, ...]) -> str:
    """Consulta DIRIGIDA: "existe alguma bula deste principio ativo com algum
    destes termos em alguma destas secoes?".

    Por que nao basta procurar na janela
    ------------------------------------
    A janela e um recorte; esta busca roda no acervo INTEIRO (265.009 bulas no
    export de 2026-11-25). Medido no ciprofloxacino: 118 das 166 bulas trazem
    `neuropathy` em secao forte, e a resposta dirigida custa 202 KB contra
    27 MB de baixar todas as 166 (`limit=1000`, o maximo da API). Completa e
    ~130x mais barata.

    `meta.results.total` da resposta e o numero de bulas que casam -- o peso da
    evidencia, que uma contagem feita dentro da janela nao sabe dizer: 118 de
    166 no ciprofloxacino contra 1 de 211 na indicacao da capsaicina."""
    termo_ou = " OR ".join(f'{secao}:"{t.strip()}"'
                           for secao in secoes
                           for t in termos if t.strip())
    busca = quote(f'openfda.{campo}:"{nome}" AND ({termo_ou})', safe="")
    return f"{BASE}?search={busca}&limit=1"


def achados_dirigidos(nome: str, campo: str, termos: list[str], forca: str, *,
                      obter_json: Callable) -> dict:
    """{'total', 'achados'} da consulta dirigida de uma forca de secao.

    `total` zero significa que NENHUMA bula do principio ativo casa -- no
    acervo inteiro, nao so na janela. `achados` traz o que a bula devolvida
    justifica, no formato de `achados()`, para o rastro."""
    secoes = SECOES_POR_FORCA[forca]
    if not campo or not [t for t in termos if t.strip()]:
        return {"total": 0, "achados": []}
    resposta = obter_json(url_dirigida(nome, campo, termos, secoes)) or {}
    bulas = resposta.get("results") or []
    total = (((resposta.get("meta") or {}).get("results") or {})
             .get("total", len(bulas)))
    # A bula devolvida tambem traz as outras secoes; so a forca perguntada
    # interessa aqui, para que cada consulta responda uma pergunta.
    lista = [a for a in achados(bulas, termos) if a["forca"] == forca]
    return {"total": total, "achados": lista}


def url_de_indicacao(nome: str, campo: str, termos: list[str]) -> str:
    """Consulta DIRIGIDA: "existe alguma bula deste principio ativo que o
    indique para a condicao?".

    Por que nao basta procurar na janela de 50 bulas
    ------------------------------------------------
    Medido: a capsaicina tem 211 bulas no openFDA, e a UNICA que a indica para
    dor neuropatica (Qutenza) esta fora das 50 que a busca por nome devolve.
    Sem esta consulta, a capsaicina -- que TRATA a condicao -- e reprovada por
    `paresthesia@warnings_and_cautions`, e o grupo de controle que trata passa
    de 0/6 para 1/6 de falso positivo.

    E a mesma falha que tirou o `limit=1` deste modulo: uma bula arbitraria
    entre centenas. Aumentar a janela nao a resolve, so a adia -- 211 hoje,
    mais amanha. Perguntar e O(1) e independe do tamanho da janela."""
    termo_ou = " OR ".join(f'indications_and_usage:"{t.strip()}"'
                           for t in termos if t.strip())
    busca = quote(f'openfda.{campo}:"{nome}" AND ({termo_ou})', safe="")
    return f"{BASE}?search={busca}&limit=1"


def achados_de_indicacao(nome: str, campo: str, termos: list[str], *,
                         obter_json: Callable) -> list[dict]:
    """Achados de INDICACAO vindos da consulta dirigida, no formato de
    `achados`. Lista vazia significa que nenhuma bula do principio ativo indica
    o farmaco para a condicao -- e nao que nao se olhou."""
    return achados_dirigidos(nome, campo, termos, "indicacao",
                             obter_json=obter_json)["achados"]


def _padrao(termo: str) -> re.Pattern:
    return re.compile(r"\b" + re.escape(termo.strip()) + r"\b", re.IGNORECASE)


def trecho(texto: str, casamento: re.Match, janela: int = JANELA_TRECHO) -> str:
    inicio = max(0, casamento.start() - janela)
    fim = min(len(texto), casamento.end() + janela)
    return (("..." if inicio else "") + texto[inicio:fim].strip()
            + ("..." if fim < len(texto) else ""))


def achados(bulas: list[dict], termos: list[str]) -> list[dict]:
    """Uma linha por (termo, secao), com o primeiro trecho que a justifica.

    Agrega todas as bulas: o que importa e se ALGUMA bula do principio ativo
    traz o termo naquela secao. Deduplicado por (termo, secao) -- 50 bulas do
    mesmo generico diriam a mesma coisa 50 vezes."""
    vistos: dict[tuple[str, str], dict] = {}
    grupos = ((SECOES_FORTES, "forte"), (SECOES_FRACAS, "fraca"),
              (SECOES_INDICACAO, "indicacao"))
    padroes = [(t, _padrao(t)) for t in termos if t.strip()]
    for bula in bulas:
        for secoes, forca in grupos:
            for secao in secoes:
                for texto in bula.get(secao) or []:
                    for termo, padrao in padroes:
                        casamento = padrao.search(texto)
                        if casamento and (termo, secao) not in vistos:
                            vistos[(termo, secao)] = {
                                "termo": termo, "secao": secao, "forca": forca,
                                "trecho": trecho(texto, casamento)}
    return [vistos[c] for c in sorted(vistos)]


def veredito(resultado: dict, lista: list[dict], dirigidos: dict | None = None):
    """(veredito, valor) do criterio `alerta_bula`.

    `dirigidos` sao as contagens das consultas dirigidas, {forca: total}. Quando
    vem, elas DECIDEM, e a janela so fornece o rastro: a janela e um recorte de
    50 bulas e a consulta dirigida roda no acervo inteiro. Sem elas (uso direto
    do modulo, ou execucao sem rede para isso), vale o que a janela viu.

    * sem bula nenhuma -> `indeterminado`: ausencia de documento nao e ausencia
      de risco, e e o caso de quase todo produto natural;
    * termo na INDICACAO -> `passa`, mesmo havendo secao forte: o farmaco trata
      a condicao (pregabalina, duloxetina);
    * termo em secao FORTE -> `reprova`;
    * so em secao fraca -> `passa`, com o achado dito no valor: parestesia em
      `adverse_reactions` e quase universal e foi a origem dos QUATRO falsos
      positivos da regra antiga. O achado nao desaparece -- ele esta em
      `alertas_bula.tsv` com `forca = fraca`, para conferencia; o que ele nao
      faz e virar reprovacao quando o criterio esta promovido;
    * bula sem nenhum termo -> `passa`."""
    if not resultado.get("bulas"):
        return "indeterminado", "sem bula no openFDA"
    dirigidos = dirigidos or {}
    bulas_do_nome = ((resultado.get("total") or len(resultado["bulas"])))

    def _conta(forca: str, achados_da_forca: list[dict]) -> int:
        """Quantas bulas sustentam a forca: o acervo quando ha consulta
        dirigida, a janela quando nao ha."""
        if forca in dirigidos:
            return dirigidos[forca]
        return 1 if achados_da_forca else 0

    indicacao = [a for a in lista if a["forca"] == "indicacao"]
    if _conta("indicacao", indicacao):
        termos = ";".join(sorted({a["termo"] for a in indicacao})) or "a condicao"
        return "passa", (f"indicado para {termos}"
                         + (f" ({dirigidos['indicacao']} de {bulas_do_nome} bulas)"
                            if "indicacao" in dirigidos else ""))
    fortes = [a for a in lista if a["forca"] == "forte"]
    if _conta("forte", fortes):
        detalhe = ";".join(f"{a['termo']}@{a['secao']}" for a in fortes)
        return "reprova", (detalhe or "termo em secao forte") + (
            f" ({dirigidos['forte']} de {bulas_do_nome} bulas)"
            if "forte" in dirigidos else "")
    fracas = [a for a in lista if a["forca"] == "fraca"]
    if fracas:
        return "passa", ("so em secao fraca (nao conta; 4/6 de falso positivo "
                         "medido): "
                         + ";".join(f"{a['termo']}@{a['secao']}" for a in fracas))
    return "passa", f"{len(resultado['bulas'])} bula(s), nenhum termo"
