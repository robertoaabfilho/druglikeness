# -*- coding: utf-8 -*-
"""
reactome.py

Expansao por via: a partir dos alvos semente, traz as proteinas que dividem
uma via do Reactome com eles.

Regras
------
* So dados de ANOTACAO do Reactome (CC0): mapeamento UniProt -> vias e
  participantes de cada via. Os arquivos de Functional Interaction nao podem
  ser redistribuidos sem licenca da Pathway Solutions e dependiam de conteudo
  KEGG; nada aqui os le.
* Teto por via (padrao 150 proteinas): vias maiores sao concentradores
  genericos ("G alpha (i) signalling events", 317 proteinas) e sao descartadas
  inteiras; as especificas ficam.
* Um alvo expandido so entra se existir no espaco do modelo
  (`alvos_expandidos_no_modelo`); sem isso o preditor nao o pontua.
* Toda adicao sai com `origem = reactome`, `via_origem` e `alvo_semente`.
* Falha por item (`ListaTruncada`, `RespostaInesperada`, resposta HTTP de erro) e registrada em
  `relatorio["problemas"]` e a execucao segue.
"""

from typing import Callable

from mods.evidencias.chembl import (
    ListaTruncada, RespostaInesperada, alvos_por_uniprot, erro_de_item)

BASE = "https://reactome.org/ContentService/"
ESPECIE = "9606"
ORGANISMO = "Homo sapiens"


def _obter_json_requests(url: str):
    import requests
    resposta = requests.get(url, timeout=60)
    if resposta.status_code == 404:
        return []  # UniProt sem via no Reactome: resposta legitima, nao erro
    resposta.raise_for_status()
    return resposta.json()


def _lista(resposta, consulta: str) -> list:
    if not isinstance(resposta, list):
        raise RespostaInesperada(f"{consulta}: esperava lista, veio {type(resposta).__name__}")
    return resposta


def vias_do_uniprot(acc: str, *, obter_json: Callable = _obter_json_requests
                    ) -> list[dict]:
    """Vias humanas do Reactome que contem `acc`: `stId` e `nome`."""
    url = f"{BASE}data/mapping/UniProt/{acc}/pathways?species={ESPECIE}"
    vias = _lista(obter_json(url), f"vias de {acc}")
    return [{"stId": v["stId"], "nome": v.get("displayName", "")} for v in vias]


def proteinas_da_via(st_id: str, *, obter_json: Callable = _obter_json_requests
                     ) -> list[str]:
    """UniProts (sem sufixo de isoforma) dos participantes da via, sem repeticao
    e na ordem em que o Reactome os lista. Moleculas pequenas ficam de fora."""
    url = f"{BASE}data/participants/{st_id}/referenceEntities"
    achados: dict[str, None] = {}
    for e in _lista(obter_json(url), f"participantes de {st_id}"):
        if e.get("databaseName") == "UniProt" and e.get("identifier"):
            achados[e["identifier"].split("-")[0]] = None
    return list(achados)


def expandir_vias(sementes: list[dict], *, teto: int = 150,
                  obter_json: Callable = _obter_json_requests,
                  relatorio: dict | None = None) -> list[dict]:
    """Proteinas que dividem via (de ate `teto` proteinas) com as sementes.
    Cada semente traz `uniprot` e `alvo_chembl_id`. Uma linha por UniProt
    (as sementes nao voltam): `via_origem` e a primeira via que a trouxe e
    `alvo_semente` junta, com ';', todas as sementes que a alcancaram.
    `relatorio` recebe `vias_vistas`, `vias_descartadas` ({stId: n}) e
    `problemas`."""
    rel = relatorio if relatorio is not None else {}
    rel.setdefault("vias_vistas", 0)
    rel.setdefault("vias_descartadas", {})
    rel.setdefault("problemas", [])
    proprias = {s["uniprot"] for s in sementes if s.get("uniprot")}
    cache: dict[str, list[str] | None] = {}
    achados: dict[str, dict] = {}
    for semente in sementes:
        acc = semente.get("uniprot")
        if not acc:
            continue
        try:
            vias = vias_do_uniprot(acc, obter_json=obter_json)
        except (ListaTruncada, OSError, ValueError) as e:
            if not erro_de_item(e):
                raise
            rel["problemas"].append(f"{semente.get('alvo_chembl_id')}: {e}")
            continue
        for via in vias:
            st_id = via["stId"]
            if st_id not in cache:
                try:
                    cache[st_id] = proteinas_da_via(st_id, obter_json=obter_json)
                except (ListaTruncada, OSError, ValueError) as e:
                    if not erro_de_item(e):
                        raise
                    rel["problemas"].append(f"{st_id}: {e}")
                    cache[st_id] = None
                    continue
                rel["vias_vistas"] += 1
                if len(cache[st_id]) > teto:
                    rel["vias_descartadas"][st_id] = len(cache[st_id])
            proteinas = cache[st_id]
            if proteinas is None or len(proteinas) > teto:
                continue
            for prot in proteinas:
                if prot in proprias:
                    continue
                if prot not in achados:
                    achados[prot] = {"uniprot": prot, "origem": "reactome",
                                     "via_origem": st_id, "via_nome": via["nome"],
                                     "alvo_semente": semente["alvo_chembl_id"]}
                else:
                    sementes_ja = achados[prot]["alvo_semente"].split(";")
                    if semente["alvo_chembl_id"] not in sementes_ja:
                        achados[prot]["alvo_semente"] += ";" + semente["alvo_chembl_id"]
    return list(achados.values())


def alvos_expandidos_no_modelo(expandidos: list[dict], rotulos: list[dict], *,
                               obter_json_chembl: Callable[[str], dict],
                               ja_presentes: set[str] = frozenset(),
                               relatorio: dict | None = None,
                               organismo: str = ORGANISMO) -> list[dict]:
    """Segundo filtro da expansao: so entra o que existe no espaco do modelo.
    `ja_presentes` sao ChEMBL IDs que a lista ja tem (usuario, Open Targets):
    nao se repetem. `relatorio` recebe a contagem de cada corte."""
    rel = relatorio if relatorio is not None else {}
    for k in ("proteinas", "sem_alvo_chembl", "fora_do_modelo", "ja_na_lista",
              "no_modelo"):
        rel.setdefault(k, 0)
    rel.setdefault("problemas", [])
    modelo = {r["chembl_id"]: r for r in rotulos}
    saida = []
    for e in expandidos:
        rel["proteinas"] += 1
        try:
            alvos = alvos_por_uniprot(e["uniprot"], organismo,
                                      obter_json=obter_json_chembl)
        except (ListaTruncada, OSError, ValueError) as erro:
            if not erro_de_item(erro):
                raise
            rel["problemas"].append(f"{e['uniprot']}: {erro}")
            continue
        if not alvos:
            rel["sem_alvo_chembl"] += 1
            continue
        alvo = alvos[0]
        chembl_id = alvo["target_chembl_id"]
        if chembl_id not in modelo:
            rel["fora_do_modelo"] += 1
            continue
        if chembl_id in ja_presentes:
            rel["ja_na_lista"] += 1
            continue
        rel["no_modelo"] += 1
        saida.append({**e, "alvo_chembl_id": chembl_id,
                      "alvo": alvo.get("pref_name") or "",
                      "organismo": alvo.get("organism") or organismo,
                      "tipo_alvo": alvo.get("target_type") or "",
                      "no_modelo": True,
                      "n_moleculas_treino": modelo[chembl_id]["n_ativos"]})
    return saida
