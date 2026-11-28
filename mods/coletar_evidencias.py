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

Codigos de saida
----------------
0 execucao completa; 3 execucao PARCIAL (alguma linha de entrada descartada ou
alguma falha de item: a saida existe, mas esta incompleta e o manifesto diz
isso); 1 nenhuma entrada utilizavel; 2 erro (nada e publicado).

Publicacao das saidas
---------------------
As tres saidas (tabela, lista enriquecida, manifesto) sao preparadas em
temporarios e so entao publicadas, juntas. Tabela sem manifesto nao e citavel,
e era o que sobrava quando a escrita do manifesto falhava depois da troca.
"""

import argparse
import importlib.util
import pickle
import sys
from datetime import datetime, timezone
from functools import partial
from pathlib import Path

# Permite `python mods/coletar_evidencias.py`: a raiz do repositorio entra no
# caminho para que `mods.evidencias.*` resolva (com `python -m` ja resolve).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from mods.evidencias import atomico, procedencia
from mods.evidencias.cache import FaltaNoCache, VersaoIndeterminada
from mods.evidencias.entrada import (
    ORIGENS_DE_ENRIQUECIMENTO, ler_alvos, ler_ligantes, ler_sintomas)
from mods.evidencias.esquema import EscritorEvidencias
from mods.evidencias.enriquecimento import reactome_fora_do_modelo
from mods.evidencias.manifesto import (
    hashes_de_entrada_ausentes, ler_manifesto, texto_do_manifesto,
    verificar_manifesto)
from mods.evidencias.resumo import teto_do_modelo

# Codigos de saida (ver docstring do modulo).
SAIDA_OK = 0
SAIDA_SEM_ENTRADA = 1
SAIDA_ERRO = 2
SAIDA_PARCIAL = 3


# Procedencia (commit, repositorio, versoes, MD5) e compartilhada com a etapa
# 2: duas etapas da mesma analise nao podem registrar quem rodou de formas
# diferentes. Os aliases abaixo existem porque este modulo ja expunha esses
# nomes.
_git = procedencia._git
_commit_do_pipeline = procedencia.commit_do_pipeline
_repositorio_do_pipeline = procedencia.repositorio_do_pipeline
_versoes_do_ambiente = procedencia.versoes_do_ambiente
_md5 = procedencia.md5


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
        return SAIDA_ERRO

    try:
        ligantes, p_lig = ler_ligantes(args.ligantes, args.limite)
        alvos, p_alvo = ler_alvos(args.alvos)
        sintomas = ler_sintomas(args.sintomas)
    except (OSError, UnicodeDecodeError) as e:
        print(f"[erro] nao foi possivel ler a entrada: {e}", file=sys.stderr)
        return SAIDA_ERRO
    _relatar("ligantes", p_lig)
    _relatar("alvos", p_alvo)

    vazias = [n for n, itens in (("ligantes", ligantes), ("alvos", alvos)) if not itens]
    for n in vazias:
        print(f"[erro] {n}: nenhum registro utilizavel; execucao interrompida",
              file=sys.stderr)
    if vazias:
        return SAIDA_SEM_ENTRADA
    if not sintomas:
        print("[aviso] sintomas: nenhum sintoma na entrada; NENHUM enriquecimento "
              "sera feito (so os alvos informados serao usados)", file=sys.stderr)

    efos: dict[str, str] = {}
    for par in args.efo:
        sintoma, _, efo = par.partition("=")
        if not sintoma.strip() or not efo.strip():
            print(f"[erro] --efo espera 'SINTOMA=EFO_ID', recebido {par!r}", file=sys.stderr)
            return SAIDA_ERRO
        efos[sintoma.strip()] = efo.strip()
    if args.alvos_enriquecidos and not args.rotulos:
        print("[erro] --alvos-enriquecidos precisa de --rotulos (o espaco de alvos "
              "do modelo decide o que a expansao pode acrescentar)", file=sys.stderr)
        return SAIDA_ERRO
    if args.alvos_enriquecidos and args.sem_rede:
        print("[erro] --alvos-enriquecidos precisa de rede; incompativel com --sem-rede",
              file=sys.stderr)
        return SAIDA_ERRO
    if args.alvos_enriquecidos:
        destino_enr = Path(args.alvos_enriquecidos).expanduser().resolve()
        if destino_enr == Path(args.alvos).expanduser().resolve():
            print("[erro] --alvos e --alvos-enriquecidos apontam para o mesmo "
                  f"arquivo ({args.alvos_enriquecidos}): a execucao leria a lista "
                  "revisada e a reescreveria por cima da sua revisao. Grave a lista "
                  "nova em outro caminho", file=sys.stderr)
            return SAIDA_ERRO
        if destino_enr == Path(args.saida).expanduser().resolve():
            print("[erro] --saida e --alvos-enriquecidos apontam para o mesmo "
                  f"arquivo ({args.saida}): sao duas saidas diferentes e uma "
                  "apagaria a outra", file=sys.stderr)
            return SAIDA_ERRO
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
            return SAIDA_ERRO
        ausentes = hashes_de_entrada_ausentes(anterior)
        if ausentes:
            print(f"[erro] repeticao recusada: o manifesto {args.manifesto} e anterior "
                  f"ao registro de MD5 das entradas e nao permite conferir os arquivos "
                  f"de entrada (ausentes: {', '.join(ausentes)}). Gere um manifesto novo "
                  "com uma execucao normal", file=sys.stderr)
            return SAIDA_ERRO
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
            return SAIDA_ERRO
        try:
            from mods.modelo import carregar_rotulos
            rotulos = carregar_rotulos(args.rotulos)
        except (OSError, RuntimeError, ValueError) as e:
            print(f"[erro] rotulos: {e}", file=sys.stderr)
            return SAIDA_ERRO
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
    except (OSError, RuntimeError, ValueError, VersaoIndeterminada) as e:
        print(f"[erro] cache de rede: {e}", file=sys.stderr)
        return SAIDA_ERRO
    md5_ou_none = procedencia.md5_ou_none
    # O dispositivo entra no manifesto porque a GPU nao reproduz (36 de 312
    # colunas diferem entre duas execucoes de GPU, 0 de 312 na CPU). Esta
    # etapa roda a previsao sempre na CPU, e o manifesto PROVA qual foi -- sem
    # isso, duas tabelas diferentes seriam indistinguiveis, e o modo de
    # repeticao nao detectaria a diferenca.
    md5s = {"dispositivo": "cpu",
            "checkpoint_md5": md5_ou_none(args.checkpoint),
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
            return SAIDA_ERRO
    if args.alvos_enriquecidos:
        # Antes de tocar nas saidas: erro de rede aqui nao deixa tabela nenhuma.
        try:
            linhas_enriquecidas = _enriquecer(args, alvos, sintomas, efos, rotulos,
                                              contagens_enriquecimento, cache)
        except (OSError, RuntimeError, ValueError, FaltaNoCache) as e:
            print(f"[erro] enriquecimento: {e}", file=sys.stderr)
            return SAIDA_ERRO

    if args.evidencias_chembl and not Path(args.evidencias_chembl).expanduser().is_file():
        print(f"[erro] --evidencias-chembl: arquivo nao encontrado: "
              f"{args.evidencias_chembl}", file=sys.stderr)
        return SAIDA_ERRO
    contagens_medido: dict = {}
    if args.evidencias_chembl:
        from mods.evidencias.medido import validar_cabecalho
        try:
            validar_cabecalho(args.evidencias_chembl)
        except (OSError, ValueError, UnicodeDecodeError) as e:
            print(f"[erro] evidencia medida: {e}", file=sys.stderr)
            return SAIDA_ERRO

    # Escrita atomica: cada saida nasce num temporario NO MESMO diretorio
    # (mesmo sistema de arquivos) e so toma o lugar do destino quando TODAS
    # estao prontas. Qualquer falha deixa os arquivos anteriores intactos, e
    # nunca se publica uma tabela sem o manifesto que a torna citavel.
    pendentes: list[tuple[Path, Path]] = []
    publicados: list[Path] = []
    # Quais alvos receberam evidencia MEDIDA. So os IDs, nunca as linhas: numa
    # execucao real sao ~23 milhoes de linhas, e o teto do modelo precisa
    # apenas de "este alvo recebeu medida, sim ou nao".
    alvos_com_medida: set[str] = set()
    etapa = "evidencia prevista"
    try:
        tabela = atomico.temporario(saida)
        pendentes.append((tabela, saida))
        with EscritorEvidencias(tabela) as escritor:
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
                    ident = str(linha.get("alvo_chembl_id") or "").strip()
                    if ident:
                        alvos_com_medida.add(ident)
            # as demais fontes de evidencia entram nas proximas tarefas

        for problema in contagens_modelo.get("problemas", []):
            print(f"[aviso] evidencia prevista, {problema}", file=sys.stderr)
        if args.alvos_enriquecidos:
            for problema in contagens_enriquecimento.get("problemas", []):
                print(f"[aviso] enriquecimento, {problema}", file=sys.stderr)
        for problema in contagens_medido.get("problemas", []):
            print(f"[aviso] evidencia medida, {problema}", file=sys.stderr)
        # deduplicacoes deliberadas: informadas, mas nao tornam a saida parcial
        for nome_fonte, c in (("evidencia prevista", contagens_modelo),
                              ("evidencia medida", contagens_medido)):
            for aviso in c.get("avisos", []):
                print(f"[info] {nome_fonte}, {aviso}", file=sys.stderr)
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

        if args.alvos_enriquecidos:
            etapa = "alvos enriquecidos"
            from mods.evidencias.open_targets import texto_dos_alvos_enriquecidos
            destino_enriquecidos = Path(args.alvos_enriquecidos)
            destino_enriquecidos.parent.mkdir(parents=True, exist_ok=True)
            arquivo = atomico.temporario(destino_enriquecidos)
            pendentes.append((arquivo, destino_enriquecidos))
            arquivo.write_text(texto_dos_alvos_enriquecidos(linhas_enriquecidas),
                               encoding="utf-8")

        etapa = "manifesto"
        descartadas = sum(p.gravidade == "erro" for p in p_lig + p_alvo)
        falhas_de_item = {
            "previsto": list(contagens_modelo.get("problemas", [])),
            "enriquecimento": list(contagens_enriquecimento.get("problemas", [])),
            "medido": list(contagens_medido.get("problemas", []))}
        parcial = bool(descartadas) or any(falhas_de_item.values())
        # Derivado do que a execucao usou, nao do que a linha de comando pediu:
        # a evidencia usou enriquecimento se a LISTA DE ALVOS lida carrega
        # procedencia de enriquecimento. Passar --alvos-enriquecidos nao faz a
        # evidencia usa-la, e usa-la sem passar a flag era registrado como
        # 'false' -- o manifesto afirmava o contrario da verdade.
        origens_dos_alvos: dict[str, int] = {}
        for a in alvos:
            origem = a.get("origem_alvo") or "(sem origem)"
            origens_dos_alvos[origem] = origens_dos_alvos.get(origem, 0) + 1
        enriquecido_usado = any(o in ORIGENS_DE_ENRIQUECIMENTO for o in origens_dos_alvos)
        arquivo = atomico.temporario(manifesto)
        pendentes.append((arquivo, manifesto))
        arquivo.write_text(texto_do_manifesto(
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
            versoes={**versoes_rede, **md5s, **_versoes_do_ambiente(),
                     **_repositorio_do_pipeline(),
                     "pipeline_commit": _commit_do_pipeline(),
                     "data": datetime.now(timezone.utc).isoformat(timespec="seconds")},
            contagens={"ligantes": len(ligantes), "alvos": len(alvos),
                       "sintomas": len(sintomas),
                       "problemas_ligantes": len(p_lig), "problemas_alvos": len(p_alvo),
                       "linhas_descartadas": descartadas,
                       "execucao_parcial": parcial,
                       "smiles_ilegiveis_modelo": ilegiveis,
                       "linhas_evidencia": n_linhas,
                       "medido_problemas": len(contagens_medido.get("problemas", [])),
                       "falhas_de_item": falhas_de_item,
                       "medido_linhas_lidas": contagens_medido.get("linhas_lidas", 0),
                       "alvos_enriquecidos_open_targets": contagens_enriquecimento.get("open_targets", 0),
                       "alvos_enriquecidos_escritos": len(linhas_enriquecidas) if args.alvos_enriquecidos else 0,
                       "origem_dos_alvos_da_evidencia": origens_dos_alvos,
                       "alvos_enriquecidos_usados_na_evidencia": enriquecido_usado,
                       "alvos_enriquecidos_reactome": contagens_enriquecimento.get("reactome", 0),
                       "medido_linhas_emitidas": contagens_medido.get("linhas_emitidas", 0),
                       # Quanto do espectro de alvos o modelo atual NAO cobre:
                       # e o numero que justifica -- ou dispensa -- retreinar o
                       # preditor, medido a cada execucao em vez de citado de
                       # uma execucao antiga. A lista e a que a EVIDENCIA usou.
                       "teto_do_modelo": teto_do_modelo(
                           alvos, alvos_com_medida, n_ligantes=len(ligantes),
                           rotulos=rotulos,
                           dobrados=contagens_enriquecimento.get(
                               "dobrados_por_familia", 0)),
                       # E o teto da lista que ESTA execucao escreveu, quando
                       # ela enriqueceu: na execucao que faz o enriquecimento,
                       # o teto acima mede a lista PRE-enriquecimento, e quem
                       # citasse so aquele numero concluiria que o teto e zero.
                       "teto_da_lista_enriquecida": teto_do_modelo(
                           linhas_enriquecidas, alvos_com_medida,
                           n_ligantes=len(ligantes), rotulos=rotulos,
                           dobrados=contagens_enriquecimento.get(
                               "dobrados_por_familia", 0))
                       if linhas_enriquecidas else None,
                       "reactome_fora_do_modelo_entraram":
                           reactome_fora_do_modelo(contagens_enriquecimento),
                       "dobrados_por_familia_fora_do_modelo":
                           contagens_enriquecimento.get(
                               "dobrados_fora_do_modelo", [])}),
            encoding="utf-8")

        # So aqui algum arquivo do usuario e substituido, e todos de uma vez.
        etapa = "publicacao"
        # Confere TODOS os destinos antes de trocar QUALQUER um: um
        # travamento no meio do laco deixaria tabela sem manifesto.
        atomico.conferir_publicavel([d for _, d in pendentes])
        for arquivo, destino in pendentes:
            atomico.publicar(arquivo, destino)
            publicados.append(destino)
    except (OSError, RuntimeError, ValueError, ImportError, EOFError,
            pickle.UnpicklingError, FaltaNoCache) as e:
        if etapa == "publicacao":
            # Troca de nome que falha no meio: dizer exatamente o que ficou de
            # fora, em vez de deixar o usuario supor que a saida esta completa.
            faltaram = [str(d) for _, d in pendentes if d not in publicados]
            print(f"[erro] publicacao incompleta: {e}. Publicado(s): "
                  f"{', '.join(str(d) for d in publicados) or 'nenhum'}; NAO "
                  f"publicado(s): {', '.join(faltaram)}", file=sys.stderr)
        else:
            # Execucao falha nao deixa saida nova nenhuma.
            print(f"[erro] {etapa}: {e}", file=sys.stderr)
        return SAIDA_ERRO
    finally:
        # So os temporarios; o que foi publicado ja nao esta aqui.
        for arquivo, _ in pendentes:
            arquivo.unlink(missing_ok=True)

    if anterior is not None:
        antes = anterior.get("contagens", {}).get("falhas_de_item")
        agora = falhas_de_item
        if antes is None:
            print("[aviso] repeticao: o manifesto nao registra as falhas de item da "
                  "execucao original; nao ha como saber se esta as reproduziu",
                  file=sys.stderr)
        elif antes != agora:
            print("[ATENCAO] repeticao: as falhas de item diferem das da execucao "
                  f"original ({sum(map(len, antes.values()))} registradas, "
                  f"{sum(map(len, agora.values()))} agora); a saida pode diferir",
                  file=sys.stderr)
        elif any(antes.values()):
            print(f"[aviso] repeticao reproduziu uma execucao que ja tinha "
                  f"{sum(map(len, antes.values()))} falha(s) de item (registradas no "
                  "manifesto): a saida repete essas ausencias, nao as resolve",
                  file=sys.stderr)
    if parcial:
        n_falhas = sum(map(len, falhas_de_item.values()))
        print(f"[ATENCAO] execucao PARCIAL (codigo {SAIDA_PARCIAL}): "
              f"{descartadas} linha(s) de entrada descartada(s) por erro e "
              f"{n_falhas} falha(s) de item (detalhes acima). A saida existe, mas "
              "esta incompleta; o manifesto registra isso em execucao_parcial",
              file=sys.stderr)
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
    return SAIDA_PARCIAL if parcial else SAIDA_OK


if __name__ == "__main__":
    sys.exit(main())
