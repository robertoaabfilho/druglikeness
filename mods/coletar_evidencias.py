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
from functools import partial
from pathlib import Path

# Permite `python mods/coletar_evidencias.py`: a raiz do repositorio entra no
# caminho para que `mods.evidencias.*` resolva (com `python -m` ja resolve).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from mods.evidencias.entrada import ler_alvos, ler_ligantes, ler_sintomas
from mods.evidencias.esquema import EscritorEvidencias
from mods.evidencias.manifesto import (
    escrever_manifesto, hashes_de_entrada_ausentes, ler_manifesto,
    verificar_manifesto)


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


def _abrir_cache(args, saida: Path, anterior: dict | None):
    """Cache de rede da execucao: (Cache, versoes dos servicos).
    Execucao normal: consulta as versoes ao vivo e usa a pasta delas.
    Repeticao (`anterior` = manifesto lido): so leitura, na pasta que o
    manifesto aponta; nunca toca a rede."""
    from mods.evidencias import cache as c
    if anterior is None:
        raiz = Path(args.cache) if args.cache else saida.parent / "cache_rede"
        cache, versoes = c.abrir_ao_vivo(raiz)
    else:
        par = anterior.get("parametros", {})
        pasta = par.get("cache_pasta")
        raiz = args.cache or par.get("cache_raiz")
        if not pasta or not raiz:
            raise c.VersaoIndeterminada(
                "o manifesto nao aponta para um cache de rede; so manifestos de "
                "execucoes que usaram a rede (com cache) podem ser repetidos")
        cache, versoes = c.abrir_para_repeticao(raiz, pasta)
    return cache, versoes, str(Path(raiz).resolve()), cache.diretorio.name


def _enriquecer(args, alvos, sintomas, efos, rotulos, contagens, cache) -> list[dict]:
    """Sintoma -> EFO (so os fixados com --efo; os demais listam candidatos e
    ficam de fora) e depois Open Targets + Reactome."""
    from mods.evidencias import open_targets
    from mods.evidencias.enriquecimento import enriquecer
    validos: dict[str, str] = {}
    for sintoma in sintomas:
        if sintoma in efos:
            validos[sintoma] = efos[sintoma]
            continue
        cands = open_targets.buscar_efo(sintoma, obter_json=cache.obter_json)
        lista = "; ".join(f"{c['id']} ({c['name']})" for c in cands[:5]) or "nenhum"
        print(f"[aviso] sintoma {sintoma!r} SEM EFO fixado: nao enriquecido. "
              f"Candidatos: {lista}. Escolha com --efo '{sintoma}=<EFO_ID>'",
              file=sys.stderr)
    for sintoma in efos.keys() - set(sintomas):
        print(f"[aviso] --efo para {sintoma!r}, que nao esta em --sintomas: ignorado",
              file=sys.stderr)
    return enriquecer(
        alvos, validos, rotulos, corte_associacao=args.corte_associacao,
        teto_via=args.teto_via, obter_json_chembl=cache.obter_json,
        obter_json_open_targets=cache.obter_json,
        obter_json_reactome=cache.obter_json, relatorio=contagens)


def _argumentos(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--ligantes", required=True, help="TSV com nome e smiles")
    ap.add_argument("--alvos", required=True, help="TSV com alvo_chembl_id, ou alvo + organismo")
    ap.add_argument("--sintomas", required=True, help="um sintoma por linha")
    ap.add_argument("--saida", required=True, help="evidencias.tsv a escrever")
    ap.add_argument("--manifesto", metavar="ARQUIVO",
                    help="MODO DE REPETICAO: manifesto.json de uma execucao anterior. "
                         "Recusa rodar se a versao do Open Targets/ChEMBL/Reactome ou o "
                         "MD5 do checkpoint/rotulos/evidencias diferirem do registrado, ou "
                         "se um parametro mudou; serve toda resposta de rede do cache "
                         "que ele aponta (sem rede; falta no cache e erro). O novo "
                         "manifesto sai ao lado da saida e nunca sobrescreve este")
    ap.add_argument("--cache", metavar="DIR",
                    help="raiz do cache de rede (padrao: cache_rede/ ao lado da saida; em "
                         "repeticao, a do manifesto). Subpasta por versao dos servicos")
    ap.add_argument("--corte", type=float, default=0.1,
                    help="corte de probabilidade do modelo (padrao 0.1)")
    ap.add_argument("--limite", type=_inteiro_positivo, help="usar so os N primeiros ligantes")
    ap.add_argument("--checkpoint", help="pesos .pt do PP4CKD; liga a evidencia prevista (use junto com --rotulos)")
    ap.add_argument("--rotulos", help="target_labels.tsv do PP4CKD (use junto com --checkpoint; um nao funciona sem o outro)")
    ap.add_argument("--evidencias-chembl", help="evidencias_chembl.csv (ensaios medidos); liga a evidencia medida")
    ap.add_argument("--alvos-enriquecidos",
                    help="alvos_enriquecidos.tsv a escrever: alvos do usuario mais os do "
                         "Open Targets (por sintoma) e da expansao por via do Reactome; "
                         "exige --rotulos")
    ap.add_argument("--efo", action="append", default=[], metavar="SINTOMA=EFO_ID",
                    help="fixa o EFO de um sintoma (repetivel); sintoma sem --efo NAO e "
                         "enriquecido, so tem seus candidatos listados")
    ap.add_argument("--corte-associacao", type=float, default=0.4,
                    help="escore minimo de associacao do Open Targets (padrao 0.4, "
                         "ajustavel; vai no manifesto). Escolhe pelo escore do proprio "
                         "alvo, nao por acaso de via: o 'origem' nao deve depender de o "
                         "Reactome ter alcancado o alvo. Na dor neuropatica: 30 alvos em "
                         ">=0.5, 44 em >=0.4 (canais de sodio e adrenoceptores alfa-2)")
    ap.add_argument("--teto-via", type=int, default=150,
                    help="vias do Reactome com mais proteinas que isto sao ignoradas (padrao 150)")
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

    efos: dict[str, str] = {}
    for par in args.efo:
        sintoma, _, efo = par.partition("=")
        if not sintoma.strip() or not efo.strip():
            print(f"[erro] --efo espera 'SINTOMA=EFO_ID', recebido {par!r}", file=sys.stderr)
            return 2
        efos[sintoma.strip()] = efo.strip()
    if args.alvos_enriquecidos and not args.rotulos:
        print("[erro] --alvos-enriquecidos precisa de --rotulos (o espaco de alvos "
              "do modelo decide o que a expansao pode acrescentar)", file=sys.stderr)
        return 2
    if args.alvos_enriquecidos and args.sem_rede:
        print("[erro] --alvos-enriquecidos precisa de rede; incompativel com --sem-rede",
              file=sys.stderr)
        return 2
    if efos and not args.alvos_enriquecidos:
        print("[aviso] --efo ignorado: sem --alvos-enriquecidos nao ha enriquecimento",
              file=sys.stderr)

    saida = Path(args.saida)
    anterior = None
    manifesto = saida.parent / "manifesto.json"
    if args.manifesto:
        try:
            anterior = ler_manifesto(args.manifesto)
        except (OSError, ValueError) as e:
            print(f"[erro] --manifesto: nao foi possivel ler {args.manifesto}: {e}",
                  file=sys.stderr)
            return 2
        ausentes = hashes_de_entrada_ausentes(anterior)
        if ausentes:
            print(f"[erro] repeticao recusada: o manifesto {args.manifesto} e anterior "
                  f"ao registro de MD5 das entradas e nao permite conferir os arquivos "
                  f"de entrada (ausentes: {', '.join(ausentes)}). Gere um manifesto novo "
                  "com uma execucao normal", file=sys.stderr)
            return 2
        if manifesto.resolve() == Path(args.manifesto).resolve():
            manifesto = saida.parent / "manifesto_repeticao.json"

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

    contagens_enriquecimento: dict = {}
    linhas_enriquecidas: list[dict] = []
    # Cache e versoes ANTES de qualquer outra chamada de rede: sem versao
    # conhecida (ou, em repeticao, com versao divergente) nada roda.
    cache = None
    versoes_rede: dict = {}
    pasta_cache = raiz_cache = None
    usa_rede = bool(args.alvos_enriquecidos) or (
        bool(args.evidencias_chembl) and not args.sem_rede)
    try:
        if usa_rede:
            cache, versoes_rede, raiz_cache, pasta_cache = _abrir_cache(
                args, saida, anterior)
    except (OSError, RuntimeError, ValueError) as e:
        print(f"[erro] cache de rede: {e}", file=sys.stderr)
        return 2
    def md5_ou_none(caminho):
        # Arquivo ilegivel aqui nao e o erro a relatar: a etapa que o usa o faz.
        try:
            return _md5(os.path.expanduser(caminho)) if caminho else None
        except OSError:
            return None
    md5s = {"checkpoint_md5": md5_ou_none(args.checkpoint),
            "rotulos_md5": md5_ou_none(args.rotulos),
            "evidencias_chembl_md5": md5_ou_none(args.evidencias_chembl),
            "ligantes_md5": md5_ou_none(args.ligantes),
            "alvos_md5": md5_ou_none(args.alvos),
            "sintomas_md5": md5_ou_none(args.sintomas)}
    parametros_de_saida = {"corte": args.corte, "limite": args.limite,
                           "sem_rede": args.sem_rede,
                           "corte_associacao": args.corte_associacao,
                           "teto_via": args.teto_via, "efo": dict(sorted(efos.items()))}
    if anterior is not None:
        divergencias = verificar_manifesto(
            anterior, {"versoes": {**versoes_rede, **md5s},
                       "parametros": parametros_de_saida})
        if divergencias:
            print("[erro] repeticao recusada, a execucao difere do manifesto:\n  "
                  + "\n  ".join(divergencias), file=sys.stderr)
            return 2
    if args.alvos_enriquecidos:
        # Antes de tocar nas saidas: erro de rede aqui nao deixa tabela nenhuma.
        try:
            linhas_enriquecidas = _enriquecer(args, alvos, sintomas, efos, rotulos,
                                              contagens_enriquecimento, cache)
        except (OSError, RuntimeError, ValueError) as e:
            print(f"[erro] enriquecimento: {e}", file=sys.stderr)
            return 2

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
                        resolver_alvo=None if args.sem_rede else partial(
                            chembl.dados_do_alvo, obter_json=cache.obter_json),
                        resolver_ligante=None if args.sem_rede else partial(
                            chembl.resolver_ligante, obter_json=cache.obter_json),
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

    if args.alvos_enriquecidos:
        from mods.evidencias.open_targets import escrever_alvos_enriquecidos
        try:
            Path(args.alvos_enriquecidos).parent.mkdir(parents=True, exist_ok=True)
            escrever_alvos_enriquecidos(args.alvos_enriquecidos, linhas_enriquecidas)
        except OSError as e:
            print(f"[erro] alvos enriquecidos: {e}", file=sys.stderr)
            return 2
        for problema in contagens_enriquecimento.get("problemas", []):
            print(f"[aviso] enriquecimento, {problema}", file=sys.stderr)

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
    # Derivado do que a execucao fez: a evidencia so usou o enriquecido se
    # --alvos aponta para o proprio arquivo enriquecido.
    enriquecido_usado = bool(args.alvos_enriquecidos) and (
        Path(args.alvos).expanduser().resolve()
        == Path(args.alvos_enriquecidos).expanduser().resolve())
    escrever_manifesto(
        manifesto,
        parametros={**parametros_de_saida,
                    "checkpoint": args.checkpoint,
                    "rotulos": args.rotulos,
                    "evidencias_chembl": args.evidencias_chembl,
                    "evidencias_chembl_md5": md5s["evidencias_chembl_md5"],
                    "rotulos_md5": md5s["rotulos_md5"],
                    "cache_raiz": raiz_cache, "cache_pasta": pasta_cache,
                    "repeticao_de": args.manifesto,
                    "ligantes": args.ligantes,
                    "alvos": args.alvos, "sintomas": args.sintomas},
        versoes={**versoes_rede, **md5s,
                 "pipeline_commit": _commit_do_pipeline(),
                 "data": datetime.now(timezone.utc).isoformat(timespec="seconds")},
        contagens={"ligantes": len(ligantes), "alvos": len(alvos),
                   "sintomas": len(sintomas),
                   "problemas_ligantes": len(p_lig), "problemas_alvos": len(p_alvo),
                   "linhas_descartadas": descartadas,
                   "smiles_ilegiveis_modelo": ilegiveis,
                   "linhas_evidencia": n_linhas,
                   "medido_problemas": len(contagens_medido.get("problemas", [])),
                   "medido_linhas_lidas": contagens_medido.get("linhas_lidas", 0),
                   "alvos_enriquecidos_open_targets": contagens_enriquecimento.get("open_targets", 0),
                   "alvos_enriquecidos_escritos": len(linhas_enriquecidas) if args.alvos_enriquecidos else 0,
                   "alvos_enriquecidos_usados_na_evidencia": enriquecido_usado,
                   "alvos_enriquecidos_reactome": contagens_enriquecimento.get("reactome", 0),
                   "medido_linhas_emitidas": contagens_medido.get("linhas_emitidas", 0)})
    if descartadas:
        print(f"[ATENCAO] execucao PARCIAL: {descartadas} linha(s) de entrada "
              "descartada(s) por erro (detalhes acima)", file=sys.stderr)
    if args.alvos_enriquecidos and not enriquecido_usado:
        proximo = (f"python mods/coletar_evidencias.py --ligantes {args.ligantes} "
                   f"--alvos {args.alvos_enriquecidos} --sintomas {args.sintomas} "
                   f"--saida <nova_evidencias.tsv>"
                   + (f" --checkpoint {args.checkpoint}" if args.checkpoint else "")
                   + (f" --rotulos {args.rotulos}" if args.rotulos else "")
                   + (f" --evidencias-chembl {args.evidencias_chembl}"
                      if args.evidencias_chembl else ""))
        print(f"[ATENCAO] {len(linhas_enriquecidas)} alvos gravados em "
              f"{args.alvos_enriquecidos}, mas o evidencias.tsv desta execucao NAO os "
              "usou (so os alvos de --alvos). Revise o arquivo e rode a coleta com ele:\n  "
              + proximo, file=sys.stderr)
    print(f"{n_linhas} linhas de evidencia em {saida}; manifesto em {manifesto}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
