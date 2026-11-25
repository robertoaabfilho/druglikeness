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
