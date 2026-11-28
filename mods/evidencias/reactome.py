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
* Um alvo expandido entra mesmo FORA do espaco do modelo, marcado
  `no_modelo = False` (`alvos_expandidos_no_modelo`). O preditor nao o pontua
  -- e e isso que a marca diz --, mas ele recebe evidencia MEDIDA do ChEMBL, e
  medida ranqueia acima de qualquer previsao. Antes ele era descartado, e era
  a maior perda do enriquecimento.
* Toda adicao sai com `origem = reactome`, `via_origem` (TODAS as vias que a
  trouxeram, separadas por ';', como `familia_colapsados`) e `alvo_semente`
  (idem). Uma proteina trazida por tres vias nao se confunde com a marginal.
* Falha por item (`ListaTruncada`, `RespostaInesperada`, resposta HTTP de erro) e registrada em
  `relatorio["problemas"]` e a execucao segue.
"""

from typing import Callable

from mods.evidencias.esquema import acrescentar_observacao
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


def _juntar(linha: dict, campo: str, valor: str) -> None:
    """Acrescenta `valor` a lista ';' do campo, sem repetir."""
    if valor not in linha[campo].split(";"):
        linha[campo] += ";" + valor


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
    (as sementes nao voltam): `via_origem` junta, com ';', todas as vias (dentro
    do teto) que a trouxeram e `alvo_semente` todas as sementes que a alcancaram.
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
                    _juntar(achados[prot], "via_origem", st_id)
                    _juntar(achados[prot], "via_nome", via["nome"])
                    _juntar(achados[prot], "alvo_semente", semente["alvo_chembl_id"])
    return list(achados.values())


def alvos_expandidos_no_modelo(expandidos: list[dict], rotulos: list[dict], *,
                               obter_json_chembl: Callable[[str], dict],
                               ja_presentes: set[str] = frozenset(),
                               relatorio: dict | None = None,
                               organismo: str = ORGANISMO) -> list[dict]:
    """Segundo filtro da expansao: resolve cada proteina num alvo do ChEMBL.

    **Alvo fora do espaco do modelo ENTRA, marcado `no_modelo = False`**, em
    vez de ser descartado -- essa era a maior perda do enriquecimento, e o
    motivo do descarte nao se sustenta. Alvo fora do modelo nao recebe
    evidencia PREVISTA (`casar_alvos` o omite, corretamente), mas recebe
    evidencia MEDIDA, porque `medido.py` usa os rotulos so para dar NOME ao
    alvo. E medida ranqueia em `1a`/`1b`, ACIMA de qualquer previsao:
    descartar esses alvos jogava fora a classe de evidencia mais forte do
    pipeline. O custo e mais consulta ao ChEMBL, e tempo nao e restricao.

    Corrige tambem uma assimetria: declarar protegia e enriquecer nao. Alvo do
    usuario fora do modelo ficava, alvo do Open Targets fora do modelo ficava,
    e so o do Reactome saia.

    O que CONTINUA filtrando: proteina sem alvo no ChEMBL (sem ID nao se liga
    a evidencia medida nem ao rotulo do modelo), alvo que a lista ja tem
    (`ja_presentes`: usuario, Open Targets) e segunda proteina que cai no
    mesmo alvo (grupo, complexo), que vira rastro na linha que ja existe.

    `relatorio` recebe a contagem de cada caso. `fora_do_modelo` passa a
    significar **entraram marcados**, nao "foram descartados"; `entraram` e o
    total que saiu desta funcao."""
    rel = relatorio if relatorio is not None else {}
    for k in ("proteinas", "sem_alvo_chembl", "fora_do_modelo", "ja_na_lista",
              "no_modelo", "entraram"):
        rel.setdefault(k, 0)
    rel.setdefault("problemas", [])
    modelo = {r["chembl_id"]: r for r in rotulos}
    saida = []
    por_id: dict[str, dict] = {}
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
        if chembl_id in ja_presentes:
            rel["ja_na_lista"] += 1
            continue
        if chembl_id in por_id:
            # Outra proteina caiu no mesmo alvo do ChEMBL (grupo, complexo):
            # uma linha so, com o rastro das duas.
            rel["mesmo_alvo_de_outra_proteina"] = rel.get("mesmo_alvo_de_outra_proteina", 0) + 1
            for campo in ("via_origem", "via_nome", "alvo_semente"):
                for v in e[campo].split(";"):
                    _juntar(por_id[chembl_id], campo, v)
            continue
        rotulo = modelo.get(chembl_id)
        rel["no_modelo" if rotulo else "fora_do_modelo"] += 1
        rel["entraram"] += 1
        linha = {**e, "alvo_chembl_id": chembl_id,
                 "alvo": alvo.get("pref_name") or "",
                 "organismo": alvo.get("organism") or organismo,
                 "tipo_alvo": alvo.get("target_type") or "",
                 "no_modelo": bool(rotulo),
                 # Vazio, nao zero: zero diria "treinado com zero moleculas".
                 # A distincao importa para o canal negativo, que exige >= 50.
                 "n_moleculas_treino": rotulo["n_ativos"] if rotulo else ""}
        if not rotulo:
            # Sem isto escrito na linha, zero previsoes para este alvo
            # pareceria "sem atividade" em vez de "nao foi pontuado".
            linha["observacao"] = acrescentar_observacao(
                e.get("observacao"),
                "fora do espaco do modelo: NAO recebe evidencia prevista; "
                "recebe evidencia medida (ChEMBL), que ranqueia acima de "
                "qualquer previsao")
        saida.append(linha)
        por_id[chembl_id] = saida[-1]
    return saida

