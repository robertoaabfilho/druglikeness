# -*- coding: utf-8 -*-
"""
manifesto.py

Manifesto de uma execucao: o que foi rodado, com quais versoes e quanto saiu.
JSON, legivel por humano; escrito junto da tabela de evidencias.
"""

import json
from pathlib import Path


def escrever_manifesto(caminho, *, parametros: dict, versoes: dict,
                       contagens: dict) -> None:
    manifesto = {"parametros": parametros, "versoes": versoes,
                 "contagens": contagens}
    Path(caminho).write_text(
        json.dumps(manifesto, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8")


def ler_manifesto(caminho) -> dict:
    return json.loads(Path(caminho).read_text(encoding="utf-8"))


# Versoes que, se mudarem, fazem a repeticao dar outro arquivo.
_VERSOES_DA_REPETICAO = ("open_targets", "chembl", "reactome",
                         "checkpoint_md5", "rotulos_md5", "evidencias_chembl_md5",
                         "ligantes_md5", "alvos_md5", "sintomas_md5")
# Parametros que mudam o conteudo da saida (caminhos de saida nao entram).
_PARAMETROS_DA_REPETICAO = ("corte", "limite", "sem_rede", "corte_associacao",
                            "teto_via", "efo")


def verificar_manifesto(anterior: dict, atual: dict) -> list[str]:
    """Divergencias entre o manifesto gravado e a execucao de agora; lista vazia
    significa que repetir e seguro. So o que o anterior registrou e comparado."""
    problemas = []
    for grupo, chaves in (("versoes", _VERSOES_DA_REPETICAO),
                          ("parametros", _PARAMETROS_DA_REPETICAO)):
        antes, agora = anterior.get(grupo, {}), atual.get(grupo, {})
        for chave in chaves:
            if chave in antes and antes[chave] != agora.get(chave):
                if chave.endswith("_md5"):
                    problemas.append(f"{chave}: o arquivo mudou desde o manifesto "
                                     f"(MD5 {antes[chave]!r} -> {agora.get(chave)!r})")
                    continue
                problemas.append(f"{chave}: manifesto registra {antes[chave]!r}, "
                                 f"agora {agora.get(chave)!r}")
    return problemas


_ENTRADAS_HASHEADAS = ("ligantes_md5", "alvos_md5", "sintomas_md5")


def hashes_de_entrada_ausentes(anterior: dict) -> list[str]:
    """Chaves de MD5 que o manifesto deveria ter e nao tem: as tres entradas
    sempre, e checkpoint/rotulos/evidencias medidas quando o manifesto diz que
    foram usados. Sem elas a repeticao nao tem como conferir e deve recusar."""
    versoes = anterior.get("versoes", {})
    parametros = anterior.get("parametros", {})
    exigidas = list(_ENTRADAS_HASHEADAS)
    for parametro, chave in (("checkpoint", "checkpoint_md5"), ("rotulos", "rotulos_md5"),
                             ("evidencias_chembl", "evidencias_chembl_md5")):
        if parametros.get(parametro):
            exigidas.append(chave)
    return [c for c in exigidas if not versoes.get(c)]
