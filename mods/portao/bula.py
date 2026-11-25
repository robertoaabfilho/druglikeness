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


def url(nome: str, campo: str, limite: int = LIMITE_BULAS) -> str:
    busca = quote(f'openfda.{campo}:"{nome}"', safe="")
    return f"{BASE}?search={busca}&limit={int(limite)}"


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


def veredito(resultado: dict, lista: list[dict]):
    """(veredito, valor) do criterio `alerta_bula`.

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
    indicacao = [a for a in lista if a["forca"] == "indicacao"]
    if indicacao:
        return "passa", ("indicado para " +
                         ";".join(sorted({a["termo"] for a in indicacao})))
    fortes = [a for a in lista if a["forca"] == "forte"]
    if fortes:
        return "reprova", ";".join(f"{a['termo']}@{a['secao']}" for a in fortes)
    fracas = [a for a in lista if a["forca"] == "fraca"]
    if fracas:
        return "passa", ("so em secao fraca (nao conta; 4/6 de falso positivo "
                         "medido): "
                         + ";".join(f"{a['termo']}@{a['secao']}" for a in fracas))
    return "passa", f"{len(resultado['bulas'])} bula(s), nenhum termo"
