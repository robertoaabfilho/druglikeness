#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
coletar_evidencias.py

Etapa 1: le ligantes, alvos e sintomas, valida tudo e escreve a tabela de
evidencias (`evidencias.tsv`) e o manifesto da execucao (`manifesto.json`).

Fontes ligadas: evidencia prevista (PP4CKD), com --checkpoint e --rotulos;
evidencia medida (ChEMBL), com --evidencias-chembl.
Sem elas a tabela sai so com o cabecalho. As demais fontes entram nas proximas
tarefas.

Validacao
---------
Toda a validacao roda ANTES de abrir o arquivo de saida; entrada ruim nao
deixa arquivo nenhum. Problemas `erro` de linha descartam a linha e sao
relatados; a execucao so para quando uma entrada fica sem nenhum registro
utilizavel (coluna obrigatoria ausente, arquivo vazio, todas as linhas ruins).
Problemas `aviso` sao relatados e a execucao continua.
"""

import argparse
import hashlib
import importlib.util
import os
import pickle
import subprocess
import sys
import tempfile
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


def _md5(caminho: str) -> str:
    h = hashlib.md5(usedforsecurity=False)
    with open(caminho, "rb") as f:
        for bloco in iter(lambda: f.read(1 << 20), b""):
            h.update(bloco)
    return h.hexdigest()


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
    ap.add_argument("--checkpoint", help="pesos .pt do PP4CKD; liga a evidencia prevista (use junto com --rotulos)")
    ap.add_argument("--rotulos", help="target_labels.tsv do PP4CKD (use junto com --checkpoint; um nao funciona sem o outro)")
    ap.add_argument("--evidencias-chembl", help="evidencias_chembl.csv (ensaios medidos); liga a evidencia medida")
    ap.add_argument("--sem-rede", action="store_true",
                    help="pula enriquecimento e evidencia obtida por API")
    return ap.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _argumentos(argv)
    if args.checkpoint and not args.rotulos:
        print("[erro] --checkpoint e --rotulos devem ser dados juntos "
              "(--rotulos sozinho serve so para nomear alvos dados por ID)",
              file=sys.stderr)
        return 2

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
    contagens_modelo: dict = {}
    rotulos = None
    if args.rotulos:
        # Tudo que pode falhar sem tocar na saida falha ANTES de abri-la.
        faltam = [pacote for modulo, pacote in (("torch", "torch"), ("rdkit", "rdkit"))
                  if importlib.util.find_spec(modulo) is None]
        if faltam:
            print(f"[erro] --checkpoint/--rotulos precisam de {', '.join(faltam)}: "
                  f"pip install {' '.join(faltam)}", file=sys.stderr)
            return 2
        try:
            from mods.prever_alvos import carregar_rotulos
            rotulos = carregar_rotulos(args.rotulos)
        except (OSError, RuntimeError, ValueError) as e:
            print(f"[erro] rotulos: {e}", file=sys.stderr)
            return 2
        if args.checkpoint:
            from mods.evidencias.predito import linhas_preditas

    if args.evidencias_chembl and not Path(args.evidencias_chembl).expanduser().is_file():
        print(f"[erro] --evidencias-chembl: arquivo nao encontrado: "
              f"{args.evidencias_chembl}", file=sys.stderr)
        return 2
    contagens_medido: dict = {}
    if args.evidencias_chembl:
        from mods.evidencias.medido import validar_cabecalho
        try:
            validar_cabecalho(args.evidencias_chembl)
        except (OSError, ValueError, UnicodeDecodeError) as e:
            print(f"[erro] evidencia medida: {e}", file=sys.stderr)
            return 2

    # Escrita atomica: a tabela nasce num temporario NO MESMO diretorio (mesmo
    # sistema de arquivos) e so substitui o destino ao fim de uma execucao
    # bem-sucedida. Qualquer falha deixa o evidencias.tsv anterior intacto.
    tmp = None
    etapa = "evidencia prevista"
    try:
        fd, nome_tmp = tempfile.mkstemp(dir=saida.parent, prefix=saida.name + ".",
                                        suffix=".tmp")
        os.close(fd)
        tmp = Path(nome_tmp)
        with EscritorEvidencias(tmp) as escritor:
            if args.checkpoint:
                for linha in linhas_preditas(
                        ligantes, alvos, checkpoint=args.checkpoint,
                        rotulos=rotulos, corte=args.corte,
                        contagens=contagens_modelo):
                    escritor.escrever(linha)
                    n_linhas += 1
            if args.evidencias_chembl:
                etapa = "evidencia medida"
                from mods.evidencias import chembl
                from mods.evidencias.medido import linhas_medidas
                for linha in linhas_medidas(
                        args.evidencias_chembl, ligantes, alvos, rotulos=rotulos,
                        resolver_alvo=None if args.sem_rede else chembl.dados_do_alvo,
                        resolver_ligante=None if args.sem_rede else chembl.resolver_ligante,
                        contagens=contagens_medido):
                    escritor.escrever(linha)
                    n_linhas += 1
            # as demais fontes de evidencia entram nas proximas tarefas
        # mkstemp cria 0600; o arquivo publicado deve ter o modo de uma escrita
        # normal (0666 menos a umask), igual ao do manifesto.
        mascara = os.umask(0)
        os.umask(mascara)
        try:
            os.chmod(tmp, 0o666 & ~mascara)
        except OSError:
            pass  # ajuste de cortesia: nunca custa o resultado (vfat, FUSE, SMB)
        os.replace(tmp, saida)
    except (OSError, RuntimeError, ValueError, ImportError, EOFError,
            pickle.UnpicklingError) as e:
        # Execucao falha nao deixa tabela nova nem manifesto.
        print(f"[erro] {etapa}: {e}", file=sys.stderr)
        return 2
    finally:
        if tmp is not None:
            tmp.unlink(missing_ok=True)  # so o temporario, nunca o arquivo do usuario

    for problema in contagens_medido.get("problemas", []):
        print(f"[aviso] evidencia medida, {problema}", file=sys.stderr)
    if (contagens_medido.get("linhas_lidas") and
            not contagens_medido.get("linhas_emitidas")):
        print(f"[aviso] evidencia medida: {contagens_medido['linhas_lidas']} linhas "
              "lidas de --evidencias-chembl e NENHUMA casou com os ligantes e alvos "
              "pedidos; confira o arquivo, os InChIKeys e os nomes/organismos dos alvos "
              "(zero linhas nao significa 'sem evidencia')", file=sys.stderr)
    ilegiveis = contagens_modelo.get("smiles_ilegiveis", 0)
    if ilegiveis:
        print(f"[aviso] {ilegiveis} ligante(s) com SMILES ilegivel para o modelo: "
              "NAO foram pontuados (zero linhas deles nao significa 'sem atividade')",
              file=sys.stderr)

    descartadas = sum(p.gravidade == "erro" for p in p_lig + p_alvo)
    escrever_manifesto(
        manifesto,
        parametros={"corte": args.corte, "limite": args.limite,
                    "sem_rede": args.sem_rede,
                    "checkpoint": args.checkpoint,
                    "rotulos": args.rotulos,
                    "evidencias_chembl": args.evidencias_chembl,
                    "evidencias_chembl_md5": (_md5(os.path.expanduser(args.evidencias_chembl))
                                              if args.evidencias_chembl else None),
                    "rotulos_md5": _md5(args.rotulos) if args.rotulos else None,
                    "ligantes": args.ligantes,
                    "alvos": args.alvos, "sintomas": args.sintomas},
        versoes={"pipeline_commit": _commit_do_pipeline(),
                 "data": datetime.now(timezone.utc).isoformat(timespec="seconds")},
        contagens={"ligantes": len(ligantes), "alvos": len(alvos),
                   "sintomas": len(sintomas),
                   "problemas_ligantes": len(p_lig), "problemas_alvos": len(p_alvo),
                   "linhas_descartadas": descartadas,
                   "smiles_ilegiveis_modelo": ilegiveis,
                   "linhas_evidencia": n_linhas,
                   "medido_problemas": len(contagens_medido.get("problemas", [])),
                   "medido_linhas_lidas": contagens_medido.get("linhas_lidas", 0),
                   "medido_linhas_emitidas": contagens_medido.get("linhas_emitidas", 0)})
    if descartadas:
        print(f"[ATENCAO] execucao PARCIAL: {descartadas} linha(s) de entrada "
              "descartada(s) por erro (detalhes acima)", file=sys.stderr)
    print(f"{n_linhas} linhas de evidencia em {saida}; manifesto em {manifesto}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
