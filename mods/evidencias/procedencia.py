# -*- coding: utf-8 -*-
"""
procedencia.py

De onde veio esta execucao: commit do pipeline, repositorio, versoes do
ambiente e MD5 das entradas. Sao os campos que fazem o manifesto servir de
citacao, e valem igual para toda etapa -- a etapa 2 nao pode registrar a
procedencia de um jeito diferente da etapa 1, ou dois manifestos da mesma
analise discordariam sobre qual codigo rodou.
"""

import hashlib
import os
import subprocess
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent  # mods/


def _git(args: list[str], cwd: Path) -> str | None:
    try:
        r = subprocess.run(["git", *args], capture_output=True, text=True,
                           timeout=10, cwd=cwd)
    except (OSError, subprocess.SubprocessError):
        return None
    return r.stdout.strip() if r.returncode == 0 else None


def commit_do_pipeline() -> str:
    """Hash do commit atual; 'desconhecido' sem git (ex.: pacote sem .git).

    Duas coisas que um hash sozinho esconderia:
    * arvore suja -- o sufixo '+sujo' diz que o codigo que rodou NAO e o do
      commit, e um hash limpo de uma arvore modificada e uma afirmacao falsa
      de reprodutibilidade;
    * pacote embutido em outro repositorio -- o hash seria o do repositorio de
      fora. Vem acompanhado do caminho da raiz (`pipeline_repositorio`) e do
      aviso quando essa raiz nao e a do proprio pacote.
    """
    commit = _git(["rev-parse", "HEAD"], RAIZ)
    if not commit:
        return "desconhecido"
    sujo = _git(["status", "--porcelain"], RAIZ)
    if sujo is None:
        return commit + "+estado_desconhecido"
    return commit + ("+sujo" if sujo else "")


def repositorio_do_pipeline() -> dict:
    """Onde o commit foi lido, para que um pacote embutido fique visivel."""
    raiz = _git(["rev-parse", "--show-toplevel"], RAIZ)
    esperada = RAIZ.parent  # a raiz do proprio druglikeness
    dados = {"pipeline_repositorio": raiz or "nenhum"}
    if raiz and Path(raiz).resolve() != esperada.resolve():
        dados["pipeline_repositorio_aviso"] = (
            f"o commit e do repositorio {raiz}, que nao e a raiz do pacote "
            f"({esperada}): o pacote esta embutido em outro repositorio e o "
            "hash nao identifica o codigo desta etapa")
    return dados


def versoes_do_ambiente() -> dict:
    """Python e RDKit: o manifesto e a citacao, e o RDKit decide a quimica
    (InChIKey parental, fingerprint). Caminho nao resolvido nao e problema --
    os MD5 prendem o conteudo, que e mais forte."""
    versoes = {"python": sys.version.split()[0]}
    try:
        import rdkit
        versoes["rdkit"] = getattr(rdkit, "__version__", "desconhecida")
    except ImportError:
        versoes["rdkit"] = "ausente"
    return versoes


def md5(caminho: str) -> str:
    h = hashlib.md5(usedforsecurity=False)
    with open(caminho, "rb") as f:
        for bloco in iter(lambda: f.read(1 << 20), b""):
            h.update(bloco)
    return h.hexdigest()


def md5_ou_none(caminho):
    """MD5 de um caminho opcional. Arquivo ilegivel aqui nao e o erro a
    relatar: a etapa que o usa o faz, com a mensagem dela."""
    try:
        return md5(os.path.expanduser(caminho)) if caminho else None
    except OSError:
        return None
