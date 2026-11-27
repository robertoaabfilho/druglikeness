# -*- coding: utf-8 -*-
"""
manifesto.py

Manifesto de uma execucao: o que foi rodado, com quais versoes e quanto saiu.
JSON, legivel por humano; escrito junto da tabela de evidencias.
"""

import json
from pathlib import Path

from mods.evidencias.atomico import escrever_texto


def texto_do_manifesto(*, parametros: dict, versoes: dict,
                       contagens: dict) -> str:
    """O JSON do manifesto, para quem precisa prepara-lo antes de publicar as
    outras saidas (tabela sem manifesto nao e citavel)."""
    manifesto = {"parametros": parametros, "versoes": versoes,
                 "contagens": contagens}
    return json.dumps(manifesto, indent=2, ensure_ascii=False) + "\n"


def escrever_manifesto(caminho, *, parametros: dict, versoes: dict,
                       contagens: dict) -> None:
    """Escrita atomica: o manifesto anterior (chave da repeticao) sobrevive a
    uma falha no meio da escrita."""
    escrever_texto(caminho, texto_do_manifesto(
        parametros=parametros, versoes=versoes, contagens=contagens))


def ler_manifesto(caminho) -> dict:
    return json.loads(Path(caminho).read_text(encoding="utf-8"))


# Versoes que, se mudarem, fazem a repeticao dar outro arquivo.
#
# `pipeline_commit` entrou porque uma mudanca de REGRA basta: a decisao 2 do
# espectro de alvos passou a deixar entrar na lista o alvo que o Reactome
# alcanca fora do espaco do modelo, e com isso os mesmos parametros, os mesmos
# MD5 e as mesmas versoes de servico dao outra lista. Sem esta chave, repetir
# um manifesto de antes da mudanca era aceito em silencio -- e a repeticao
# existe para RECUSAR execucao que divergiria.
_VERSOES_DA_REPETICAO = ("open_targets", "chembl", "reactome",
                         "pipeline_commit",
                         "checkpoint_md5", "rotulos_md5", "evidencias_chembl_md5",
                         "ligantes_md5", "alvos_md5", "sintomas_md5")
# Parametros que mudam o conteudo da saida (caminhos de saida nao entram).
_PARAMETROS_DA_REPETICAO = ("corte", "limite", "sem_rede", "corte_associacao",
                            "teto_via", "efo")


def verificar_manifesto(anterior: dict, atual: dict, *,
                        versoes=_VERSOES_DA_REPETICAO,
                        parametros=_PARAMETROS_DA_REPETICAO) -> list[str]:
    """Divergencias entre o manifesto gravado e a execucao de agora; lista vazia
    significa que repetir e seguro. So o que o anterior registrou e comparado.

    `versoes` e `parametros` sao as chaves que importam para a repeticao. Os
    padroes sao os da etapa 1; a etapa 2 passa os seus (politica, percentil,
    teto de bula), porque repetir com outra politica da outro arquivo."""
    problemas = []
    for grupo, chaves in (("versoes", versoes),
                          ("parametros", parametros)):
        antes, agora = anterior.get(grupo, {}), atual.get(grupo, {})
        for chave in chaves:
            if chave == "pipeline_commit" and not (
                    antes.get(chave) and agora.get(chave)):
                # Commit desconhecido (copia sem git, git ausente) e
                # INDETERMINADO, nao divergencia: recusar por falta de dado
                # trancaria a repeticao em toda arvore sem historico. A mesma
                # regra do resto do pipeline -- ausencia nao vira veredito.
                continue
            if chave in antes and antes[chave] != agora.get(chave):
                if chave.endswith("_md5"):
                    problemas.append(f"{chave}: o arquivo mudou desde o manifesto "
                                     f"(MD5 {antes[chave]!r} -> {agora.get(chave)!r})")
                    continue
                if chave == "pipeline_commit":
                    problemas.append(
                        f"pipeline_commit: o manifesto foi gravado por outro "
                        f"commit do pipeline ({str(antes[chave])[:7]} -> "
                        f"{str(agora.get(chave))[:7]}); as REGRAS podem ter "
                        "mudado, e a saida pode diferir com os mesmos "
                        "parametros")
                    continue
                problemas.append(f"{chave}: manifesto registra {antes[chave]!r}, "
                                 f"agora {agora.get(chave)!r}")
    return problemas


_ENTRADAS_HASHEADAS = ("ligantes_md5", "alvos_md5", "sintomas_md5")


_CONDICIONAIS_ETAPA1 = (("checkpoint", "checkpoint_md5"),
                        ("rotulos", "rotulos_md5"),
                        ("evidencias_chembl", "evidencias_chembl_md5"))


def hashes_de_entrada_ausentes(anterior: dict, *,
                               entradas=_ENTRADAS_HASHEADAS,
                               condicionais=_CONDICIONAIS_ETAPA1) -> list[str]:
    """Chaves de MD5 que o manifesto deveria ter e nao tem: as entradas sempre,
    e as condicionais quando o manifesto diz que foram usadas. Sem elas a
    repeticao nao tem como conferir e deve recusar. Os padroes sao os da etapa
    1; a etapa 2 passa os seus."""
    versoes = anterior.get("versoes", {})
    parametros = anterior.get("parametros", {})
    exigidas = list(entradas)
    for parametro, chave in condicionais:
        if parametros.get(parametro):
            exigidas.append(chave)
    return [c for c in exigidas if not versoes.get(c)]
