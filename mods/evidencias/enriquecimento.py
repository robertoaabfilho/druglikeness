# -*- coding: utf-8 -*-
"""
enriquecimento.py

Orquestra a segunda saida da etapa 1: lista do usuario + alvos do Open Targets
para os sintomas + expansao por via do Reactome.

Regras
------
* So acrescenta. As linhas do usuario entram primeiro, na ordem dada, e nada
  depois mexe nelas: um alvo do Open Targets ou do Reactome que ja esta na
  lista nao e repetido.
* O EFO de cada sintoma vem do chamador (`efos`); sintoma sem EFO fica de fora
  e e relatado, nunca adivinhado.
* Sementes do Reactome: os alvos do usuario e os do Open Targets que existem
  no modelo e tem um unico UniProt.
* O corte do Open Targets (padrao 0.4 na CLI) escolhe pelo escore do proprio
  alvo. Alvo que o Reactome alcanca e o Open Targets associa ABAIXO do corte
  continua `reactome` (promove-lo faria `origem` depender de a via te-lo
  alcancado) e ganha escore e `observacao`.
* Falha por item vira `observacao`/`problemas`; so erro de rede aborta.
"""

from typing import Callable

from mods.evidencias import chembl
from mods.evidencias.chembl import (
    AlvoAmbiguo, ListaTruncada, RespostaInesperada, erro_de_item)
from mods.evidencias.open_targets import (
    alvos_associados, alvos_do_chembl, colapsar_familia)
from mods.evidencias.esquema import acrescentar_observacao
from mods.evidencias.reactome import alvos_expandidos_no_modelo, expandir_vias

_ERROS_POR_ITEM = (ListaTruncada, RespostaInesperada, AlvoAmbiguo, LookupError,
                   OSError)


def _de_item(e: Exception) -> bool:
    return isinstance(e, (AlvoAmbiguo, LookupError)) or erro_de_item(e)


def _juntar_associacao(linha: dict, outra: dict) -> None:
    """Acumula sintomas e EFOs (';') e guarda o maior escore de associacao."""
    for campo in ("sintoma", "efo_id"):
        atuais = [v for v in str(linha.get(campo) or "").split(";") if v]
        if outra.get(campo) and outra[campo] not in atuais:
            linha[campo] = ";".join(atuais + [outra[campo]])
    if outra.get("score") is not None and outra["score"] > (linha.get("score") or 0):
        linha["score"] = outra["score"]


def _linha_do_usuario(a: dict, no_modelo: dict, obter_json) -> dict:
    linha = {"alvo_chembl_id": a.get("alvo_chembl_id") or "",
             "alvo": a.get("alvo") or "", "organismo": a.get("organismo") or "",
             "tipo_alvo": "", "origem": "usuario", "no_modelo": False,
             "n_moleculas_treino": "", "observacao": ""}
    try:
        if not linha["alvo_chembl_id"]:
            achado = chembl.resolver_alvo_por_nome(
                linha["alvo"], linha["organismo"], obter_json=obter_json)
            linha["alvo_chembl_id"] = achado["target_chembl_id"]
        dados = chembl.dados_do_alvo(linha["alvo_chembl_id"], obter_json=obter_json)
        linha.update(alvo=dados["nome"], organismo=dados["organismo"],
                     tipo_alvo=dados["tipo"])
    except _ERROS_POR_ITEM as e:
        if not _de_item(e):
            raise
        linha["observacao"] = f"nao resolvido no ChEMBL: {e}"
    modelo = no_modelo.get(linha["alvo_chembl_id"])
    if modelo:
        linha.update(no_modelo=True, n_moleculas_treino=modelo["n_ativos"])
    return linha


def reactome_fora_do_modelo(relatorio: dict) -> int:
    """Quantos alvos o Reactome trouxe de FORA do espaco do modelo.

    E a chave que mede o efeito da decisao 2, e por isso tem funcao propria: o
    contador mora aninhado em `relatorio["reactome_modelo"]`, e quem o lia por
    `relatorio["fora_do_modelo"]` recebia 0 em toda execucao real -- um zero
    indistinguivel de "a mudanca nao fez nada"."""
    return int((relatorio.get("reactome_modelo") or {}).get("fora_do_modelo", 0))


def enriquecer(alvos_usuario: list[dict], efos: dict[str, str], rotulos: list[dict],
               *, corte_associacao: float, teto_via: int,
               obter_json_chembl: Callable, obter_json_open_targets: Callable,
               obter_json_reactome: Callable, relatorio: dict | None = None
               ) -> list[dict]:
    """Linhas de `alvos_enriquecidos.tsv`. `efos` mapeia sintoma -> EFO ID.
    `relatorio` recebe contagens (`open_targets`, `reactome`) e `problemas`."""
    rel = relatorio if relatorio is not None else {}
    problemas = rel.setdefault("problemas", [])
    no_modelo = {r["chembl_id"]: r for r in rotulos}

    linhas = [_linha_do_usuario(a, no_modelo, obter_json_chembl) for a in alvos_usuario]
    por_id = {l["alvo_chembl_id"]: l for l in linhas if l["alvo_chembl_id"]}
    presentes = set(por_id)

    candidatos: list[dict] = []
    # UniProt -> melhor associacao ABAIXO do corte: alvo que o Reactome alcanca
    # e que o Open Targets tambem associa ao sintoma nao pode parecer so "via".
    abaixo: dict[str, dict] = {}
    for sintoma, efo in efos.items():
        try:
            todos = alvos_associados(efo, corte=0.0,
                                     obter_json=obter_json_open_targets)
        except (ListaTruncada, OSError, ValueError) as e:
            if not erro_de_item(e):
                raise
            problemas.append(f"Open Targets, {sintoma!r} ({efo}): {e}")
            continue
        for c in todos:
            if c["score"] >= corte_associacao:
                candidatos.append({**c, "sintoma": sintoma})
            elif c["uniprot"] and c["score"] > abaixo.get(c["uniprot"], {}).get("score", -1):
                abaixo[c["uniprot"]] = {**c, "sintoma": sintoma}
    rel["open_targets_candidatos"] = len(candidatos)

    resolvidos = alvos_do_chembl(candidatos, rotulos, obter_json_chembl=obter_json_chembl)
    ot: list[dict] = []
    for r in resolvidos:
        if r["alvo_chembl_id"] and r["alvo_chembl_id"] in presentes:
            # Ja na lista (do usuario ou de outro sintoma): nao repete, mas o
            # rastro de cada sintoma que o implica fica (so campos de
            # associacao; `origem` e o resto da linha nao mudam).
            _juntar_associacao(por_id[r["alvo_chembl_id"]], r)
            continue
        if r["alvo_chembl_id"]:
            presentes.add(r["alvo_chembl_id"])
            por_id[r["alvo_chembl_id"]] = r
            try:
                r["classe"] = chembl.classe_de_proteina(
                    r["alvo_chembl_id"], obter_json=obter_json_chembl)
            except _ERROS_POR_ITEM as e:
                if not _de_item(e):
                    raise
                r["classe"] = None
                r["observacao"] = (r["observacao"] + "; " if r["observacao"] else "") + f"classe: {e}"
        ot.append(r)
    rel_familia: dict = {}
    ot = colapsar_familia(ot, relatorio=rel_familia)
    rel["dobrados_por_familia"] = rel_familia["dobrados_por_familia"]
    rel["dobrados_fora_do_modelo"] = list(rel_familia["dobrados_fora_do_modelo"])
    # O alvo DOBRADO saiu da lista, entao ele nao esta "ja presente": deixa-lo
    # em `presentes` fazia o Reactome recusa-lo como "ja na lista" e ele
    # desaparecia do pipeline inteiro -- justamente a classe que a decisao 2
    # protege, alvo fora do modelo que ainda recebe evidencia medida.
    ficaram = {l["alvo_chembl_id"] for l in ot if l.get("alvo_chembl_id")}
    for ident in rel_familia["dobrados_fora_do_modelo"]:
        if ident not in ficaram:
            presentes.discard(ident)
            por_id.pop(ident, None)
    rel["open_targets"] = len(ot)
    linhas.extend(ot)

    sementes = []
    for l in linhas:
        if l["alvo_chembl_id"] and l["no_modelo"]:
            try:
                acc = chembl.uniprot_do_alvo(l["alvo_chembl_id"], obter_json=obter_json_chembl)
            except (ListaTruncada, OSError, ValueError) as e:
                if not erro_de_item(e):
                    raise
                problemas.append(f"UniProt de {l['alvo_chembl_id']}: {e}")
                continue
            if acc:
                sementes.append({"uniprot": acc, "alvo_chembl_id": l["alvo_chembl_id"]})
    rel["sementes_reactome"] = len(sementes)
    rel_vias: dict = {}
    expandidos = expandir_vias(sementes, teto=teto_via, obter_json=obter_json_reactome,
                               relatorio=rel_vias)
    rel_modelo: dict = {}
    novos = alvos_expandidos_no_modelo(
        expandidos, rotulos, obter_json_chembl=obter_json_chembl,
        ja_presentes=presentes, relatorio=rel_modelo)
    problemas.extend(rel_vias["problemas"] + rel_modelo["problemas"])
    rel["reactome_vias"] = {k: v for k, v in rel_vias.items() if k != "problemas"}
    rel["reactome_modelo"] = {k: v for k, v in rel_modelo.items() if k != "problemas"}
    for n in novos:
        a = abaixo.get(n["uniprot"])
        if a:
            n.update(efo_id=a["efo_id"], sintoma=a["sintoma"], score=a["score"])
            # ACRESCENTA, nao substitui: com `=` esta linha apagava a
            # explicacao de que o alvo esta fora do modelo, e zero previsoes
            # para ele voltaria a parecer "sem atividade".
            n["observacao"] = acrescentar_observacao(
                n.get("observacao"),
                f"tambem associado no Open Targets (escore {a['score']:.3f}, "
                f"abaixo do corte {corte_associacao})")
            rel["reactome_tambem_open_targets"] = rel.get("reactome_tambem_open_targets", 0) + 1
    rel["reactome"] = len(novos)
    linhas.extend(novos)
    return linhas
