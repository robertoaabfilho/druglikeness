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
                         "checkpoint_md5", "rotulos_md5", "evidencias_chembl_md5")
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
                problemas.append(f"{chave}: manifesto registra {antes[chave]!r}, "
                                 f"agora {agora.get(chave)!r}")
    return problemas
