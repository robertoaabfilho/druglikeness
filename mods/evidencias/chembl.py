# -*- coding: utf-8 -*-
"""
chembl.py

Ponte entre tres espacos de identificadores: o Open Targets nomeia alvos por
gene Ensembl e UniProt, o PP4CKD os indexa por ChEMBL ID, e o arquivo de
evidencia medida usa nome preferido mais organismo.

Regras
------
* Um UniProt devolve varios alvos do ChEMBL (a proteina, grupos de
  seletividade, complexos). `SINGLE PROTEIN` vem primeiro; devolver um grupo
  onde se queria a proteina mudaria em silencio o que a previsao significa.
* Por nome, nunca se escolhe sozinho: "noradrenaline transporter" existe em
  humano, camundongo e rato. Mais de um candidato levanta `AlvoAmbiguo`, que
  leva os candidatos para o chamador mostrar.
* Ligante: InChIKey exato primeiro; depois o bloco de conectividade (14
  primeiros caracteres), porque o ChEMBL indexa por forma salina e so o bloco
  reconhece cloridrato e base livre como o mesmo composto. No bloco, so se
  aceita candidato com a camada de estereo IGUAL a pedida, ou quando a pedida
  ou a do candidato e plana (UHFFFAOYSA, estereo nao informado) e o bloco tem
  um unico candidato. Enantiomeros compartilham o bloco e costumam ser farmacos
  diferentes: estereo diferente devolve None, nunca outro estereoisomero.
* Lista que a API devolveu truncada (total_count maior que a pagina) nao
  serve de base para decidir: levanta `ListaTruncada`.
* Alvo por UniProt consulta primeiro so `SINGLE PROTEIN` no servidor, para que
  a paginacao nunca corte a proteina; so sem resultado volta sem o filtro.
* Resolucao por item; nada aqui materializa o conjunto de ligantes.
* `obter_json` e injetavel (a tarefa 9 passa uma versao com cache). O padrao
  usa `requests`; os testes passam um falso que le respostas gravadas.
"""

from typing import Callable
from urllib.parse import quote

BASE = "https://www.ebi.ac.uk/chembl/api/data/"
SINGLE = "SINGLE PROTEIN"


class ListaTruncada(Exception):
    """A API devolveu so parte dos resultados; decidir sobre ela seria chutar."""


class RespostaInesperada(ValueError):
    """Corpo sem a chave esperada (ex.: erro da API): nao e "zero resultados"."""


def erro_de_item(e: Exception) -> bool:
    """O erro e de UM item (lista truncada, corpo inesperado, ou resposta HTTP
    de erro como o 500 do ChEMBL), e nao de rede fora do ar: o chamador registra
    e segue. Excecao sem resposta HTTP (timeout, DNS) nao e de item."""
    return (isinstance(e, (ListaTruncada, RespostaInesperada))
            or (isinstance(e, OSError) and getattr(e, "response", None) is not None))


def _sem_truncar(resposta: dict, chave: str, consulta: str) -> list[dict]:
    if chave not in resposta:
        raise RespostaInesperada(f"{consulta}: resposta sem {chave!r}")
    itens = resposta[chave]
    total = (resposta.get("page_meta") or {}).get("total_count")
    if total is not None and total > len(itens):
        raise ListaTruncada(f"{consulta}: {total} resultados, pagina trouxe {len(itens)}")
    return itens


class AlvoAmbiguo(Exception):
    """Mais de um alvo serve para o nome pedido; o chamador decide."""

    def __init__(self, mensagem: str, candidatos: list[dict]):
        super().__init__(mensagem)
        self.candidatos = candidatos


def _obter_json_requests(url: str) -> dict:
    import requests
    resposta = requests.get(url, timeout=60)
    resposta.raise_for_status()
    return resposta.json()


def _url(caminho: str, **filtros: str) -> str:
    consulta = "&".join(f"{k}={quote(str(v), safe='')}" for k, v in filtros.items())
    return BASE + caminho + ("?" + consulta if consulta else "")


def alvos_por_uniprot(acc: str, organismo: str, *,
                      obter_json: Callable[[str], dict] = _obter_json_requests
                      ) -> list[dict]:
    """Alvos do ChEMBL que contem o UniProt `acc`, `SINGLE PROTEIN` primeiro.
    Organismo vazio nao filtra."""
    base = {"target_components__accession": acc}
    if organismo:
        base["organism"] = organismo
    # Filtro no servidor: com varios alvos por UniProt, a pagina poderia
    # cortar a proteina antes de qualquer ordenacao local.
    consulta = f"alvos de {acc}"
    alvos = _sem_truncar(obter_json(_url("target.json", **base, target_type=SINGLE,
                                         limit="5")), "targets", consulta)
    if alvos:
        return alvos
    alvos = _sem_truncar(obter_json(_url("target.json", **base, limit="50")),
                         "targets", consulta)
    return sorted(alvos, key=lambda a: a.get("target_type") != SINGLE)


def resolver_alvo_por_nome(nome: str, organismo: str, *,
                           obter_json: Callable[[str], dict] = _obter_json_requests
                           ) -> dict:
    """O unico alvo que casa com `nome` (e `organismo`, se dado).
    Zero candidatos: LookupError. Mais de um: AlvoAmbiguo."""
    filtros = {"pref_name__icontains": nome}
    if organismo:
        filtros["organism__iexact"] = organismo
    filtros["limit"] = "20"
    achados = _sem_truncar(obter_json(_url("target.json", **filtros)),
                           "targets", f"alvos para {nome!r}")
    if not achados:
        raise LookupError(f"nenhum alvo no ChEMBL para {nome!r} ({organismo or 'qualquer organismo'})")
    if len(achados) > 1:
        # Com organismo dado, uma unica proteina isolada desempata; sem ele,
        # nunca: ortologos sao alvos diferentes.
        proteinas = [a for a in achados if a.get("target_type") == SINGLE]
        if organismo and len(proteinas) == 1:
            return proteinas[0]
        raise AlvoAmbiguo(f"{len(achados)} alvos para {nome!r}: escolha pelo ChEMBL ID",
                          achados)
    return achados[0]


def _moleculas(filtro: str, valor: str, obter_json, limite: int) -> list[dict]:
    url = _url("molecule.json", **{filtro: valor, "limit": str(limite)})
    return _sem_truncar(obter_json(url), "molecules", f"moleculas {valor}")


PLANO = "UHFFFAOYSA"


def _estereo(chave: str) -> str:
    return chave.partition("-")[2].split("-")[0]


def _chave(molecula: dict) -> str:
    return ((molecula.get("molecule_structures") or {})
            .get("standard_inchi_key") or "")


def resolver_ligante(inchikey: str, inchikey_parental: str, *,
                     obter_json: Callable[[str], dict] = _obter_json_requests
                     ) -> str | None:
    """ChEMBL ID do ligante, ou None. Exato (InChIKey, depois parental) e,
    se falhar, bloco de conectividade."""
    chaves = [c.strip().upper() for c in (inchikey, inchikey_parental) if c and c.strip()]
    for chave in chaves:
        achados = _moleculas("molecule_structures__standard_inchi_key", chave,
                             obter_json, 5)
        if achados:
            return achados[0]["molecule_chembl_id"]
    for chave in chaves:
        bloco = chave.partition("-")[0]
        if len(bloco) != 14:
            continue
        achados = _moleculas("molecule_structures__standard_inchi_key__istartswith",
                             bloco, obter_json, 20)
        pedido = _estereo(chave)
        iguais = [m for m in achados if _estereo(_chave(m)) == pedido]
        if iguais:
            return iguais[0]["molecule_chembl_id"]
        # Leniencia so quando um dos lados nao informa estereo, e so se for
        # inequivoco; estereo diferente e outra molecula.
        # Mesmo assim so se o bloco inteiro tem um unico candidato: havendo
        # outros estereoisomeros, o plano e so mais uma molecula.
        if len(achados) == 1 and PLANO in (pedido, _estereo(_chave(achados[0]))):
            return achados[0]["molecule_chembl_id"]
    return None


def uniprot_do_alvo(chembl_id: str, *,
                    obter_json: Callable[[str], dict] = _obter_json_requests
                    ) -> str | None:
    """UniProt do alvo. Complexo ou grupo (mais de um componente) devolve None:
    nao ha um unico acesso que o represente."""
    alvo = obter_json(BASE + f"target/{quote(chembl_id)}.json")
    componentes = alvo.get("target_components") or []
    if len(componentes) != 1:
        return None
    return componentes[0].get("accession") or None


def classe_de_proteina(chembl_id: str, *,
                       obter_json: Callable[[str], dict] = _obter_json_requests
                       ) -> str | None:
    """Nome da classe de proteina (a "familia" do colapso), ou None."""
    acc = uniprot_do_alvo(chembl_id, obter_json=obter_json)
    if not acc:
        return None
    componentes = obter_json(_url("target_component.json",
                                  accession=acc, limit="5")).get("target_components", [])
    for comp in componentes:
        for classe in comp.get("protein_classifications") or []:
            ident = classe.get("protein_classification_id")
            if ident is None:
                continue
            nome = obter_json(BASE + f"protein_classification/{ident}.json").get("pref_name")
            if nome:
                return nome
    return None


def dados_do_alvo(chembl_id: str, *,
                  obter_json: Callable[[str], dict] = _obter_json_requests) -> dict:
    """Nome, organismo e tipo do alvo, para casar com o arquivo de ensaios.
    Corpo sem nome ou organismo (ex.: erro da API): RespostaInesperada."""
    alvo = obter_json(BASE + f"target/{quote(chembl_id)}.json")
    if not alvo.get("pref_name") or not alvo.get("organism"):
        raise RespostaInesperada(f"alvo {chembl_id}: resposta sem pref_name/organism")
    return {"nome": alvo["pref_name"], "organismo": alvo["organism"],
            "tipo": alvo.get("target_type") or ""}
