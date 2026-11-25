#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
coletar_evidencias.py

Etapa 1: le ligantes, alvos e sintomas, valida tudo e escreve a tabela de
evidencias (`evidencias.tsv`) e o manifesto da execucao (`manifesto.json`).

Nesta versao nenhuma fonte de evidencia esta ligada ainda: a tabela sai so com
o cabecalho. As fontes entram nas proximas tarefas.

Validacao
---------
Toda a validacao roda ANTES de abrir o arquivo de saida; entrada ruim nao
deixa arquivo nenhum. Problemas `erro` de linha descartam a linha e sao
relatados; a execucao so para quando uma entrada fica sem nenhum registro
utilizavel (coluna obrigatoria ausente, arquivo vazio, todas as linhas ruins).
Problemas `aviso` sao relatados e a execucao continua.
"""

import argparse
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

# Permite `python mods/coletar_evidencias.py`: a raiz do repositorio entra no
# caminho para que `mods.evidencias.*` resolva (com `python -m` ja resolve).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from mods.evidencias.entrada import ler_alvos, ler_ligantes, ler_sintomas
from mods.evidencias.esquema import EscritorEvidencias
from mods.evidencias.manifesto import escrever_manifesto


def _commit_do_pipeline() -> str:
    """Hash do commit atual; 'desconhecido' sem git (ex.: pacote sem .git)."""
    try:
        r = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                           text=True, timeout=10,
                           cwd=Path(__file__).resolve().parent)
        if r.returncode == 0 and r.stdout.strip():
            return r.stdout.strip()
    except (OSError, subprocess.SubprocessError):
        pass
    return "desconhecido"


def _relatar(nome: str, problemas) -> None:
    for p in problemas:
        print(f"[{p.gravidade}] {nome}, linha {p.linha}, {p.campo}: {p.mensagem}",
              file=sys.stderr)


def _inteiro_positivo(texto: str) -> int:
    try:
        valor = int(texto)
    except ValueError:
        valor = 0
    if valor < 1:
        raise argparse.ArgumentTypeError(
            f"--limite deve ser um inteiro >= 1, recebido '{texto}'")
    return valor


def _argumentos(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--ligantes", required=True, help="TSV com nome e smiles")
    ap.add_argument("--alvos", required=True, help="TSV com alvo_chembl_id, ou alvo + organismo")
    ap.add_argument("--sintomas", required=True, help="um sintoma por linha")
    ap.add_argument("--saida", required=True, help="evidencias.tsv a escrever")
    ap.add_argument("--manifesto", help="padrao: manifesto.json ao lado da saida")
    ap.add_argument("--corte", type=float, default=0.1,
                    help="corte de probabilidade do modelo (padrao 0.1)")
    ap.add_argument("--limite", type=_inteiro_positivo, help="usar so os N primeiros ligantes")
    ap.add_argument("--sem-rede", action="store_true",
                    help="pula enriquecimento e evidencia obtida por API")
    return ap.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _argumentos(argv)

    try:
        ligantes, p_lig = ler_ligantes(args.ligantes, args.limite)
        alvos, p_alvo = ler_alvos(args.alvos)
        sintomas = ler_sintomas(args.sintomas)
    except (OSError, UnicodeDecodeError) as e:
        print(f"[erro] nao foi possivel ler a entrada: {e}", file=sys.stderr)
        return 2
    _relatar("ligantes", p_lig)
    _relatar("alvos", p_alvo)

    vazias = [n for n, itens in (("ligantes", ligantes), ("alvos", alvos)) if not itens]
    for n in vazias:
        print(f"[erro] {n}: nenhum registro utilizavel; execucao interrompida",
              file=sys.stderr)
    if vazias:
        return 1
    if not sintomas:
        print("[aviso] sintomas: nenhum sintoma na entrada; NENHUM enriquecimento "
              "sera feito (so os alvos informados serao usados)", file=sys.stderr)

    saida = Path(args.saida)
    manifesto = Path(args.manifesto) if args.manifesto else saida.parent / "manifesto.json"

    saida.parent.mkdir(parents=True, exist_ok=True)
    manifesto.parent.mkdir(parents=True, exist_ok=True)
    n_linhas = 0
    with EscritorEvidencias(saida):
        pass  # as fontes de evidencia entram nas proximas tarefas

    descartadas = sum(p.gravidade == "erro" for p in p_lig + p_alvo)
    escrever_manifesto(
        manifesto,
        parametros={"corte": args.corte, "limite": args.limite,
                    "sem_rede": args.sem_rede, "ligantes": args.ligantes,
                    "alvos": args.alvos, "sintomas": args.sintomas},
        versoes={"pipeline_commit": _commit_do_pipeline(),
                 "data": datetime.now(timezone.utc).isoformat(timespec="seconds")},
        contagens={"ligantes": len(ligantes), "alvos": len(alvos),
                   "sintomas": len(sintomas),
                   "problemas_ligantes": len(p_lig), "problemas_alvos": len(p_alvo),
                   "linhas_descartadas": descartadas,
                   "linhas_evidencia": n_linhas})
    if descartadas:
        print(f"[ATENCAO] execucao PARCIAL: {descartadas} linha(s) de entrada "
              "descartada(s) por erro (detalhes acima)", file=sys.stderr)
    print(f"{n_linhas} linhas de evidencia em {saida}; manifesto em {manifesto}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
