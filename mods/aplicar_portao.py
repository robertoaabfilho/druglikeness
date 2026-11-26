#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
aplicar_portao.py

Etapa 2: le os ligantes da etapa 1, calcula propriedade e ADMET, julga cada
criterio e escreve o veredito -- `propriedades.tsv`, `portao.tsv` e o
`manifesto.json` da execucao, todos dentro de `--saida`.

O nome "portao" descreve a POSICAO, nao a severidade. Medido (ver
`docs/superpowers/specs/2026-11-25-etapa2-portao-design.md`): reprovar por
PAINS apagaria a amitriptilina, que trata dor neuropatica e e o positivo de
referencia do pipeline; por BRENK, 65% da base. Entao o portao marca, anota e
explica; por padrao so SMILES ilegivel e alerta forte de bula reprovam. Quem
quiser estreita-lo promove um criterio com `--reprovar`, e o resumo diz antes
quantos ligantes isso custaria.

Criterios e politica
--------------------
Cada criterio devolve um veredito INTRINSECO (`passa`, `marca`, `reprova`,
`indeterminado`) e a politica decide o que fazer com os adversos.
`indeterminado` NUNCA tira um ligante, com politica nenhuma: ausencia de dado
nao e evidencia de problema. Criterio que a execucao nao pediu (bula sem
`--termos-alerta`) sai `indeterminado`, nunca `passa` -- "nao perguntei" e
"sem alerta" nao sao a mesma afirmacao.

Codigos de saida
----------------
0 execucao completa; 3 PARCIAL (alguma linha de entrada descartada: a saida
existe, mas esta incompleta e o manifesto diz isso); 1 nenhum ligante
utilizavel; 2 erro (nada e publicado).

Publicacao das saidas
---------------------
Cada saida nasce num temporario no mesmo diretorio e so toma o lugar do
destino quando TODAS estao prontas. Tabela sem manifesto nao e citavel.
"""

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

# Permite `python mods/aplicar_portao.py`: a raiz do repositorio entra no
# caminho para que `mods.portao.*` resolva (com `python -m` ja resolve).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from mods.evidencias import atomico, procedencia
from mods.evidencias.entrada import ler_ligantes
from mods.evidencias.manifesto import texto_do_manifesto
from mods.portao import admet, bula, decisao, esquema, propriedades

SAIDA_OK = 0
SAIDA_SEM_ENTRADA = 1
SAIDA_ERRO = 2
SAIDA_PARCIAL = 3

NOME_PROPRIEDADES = "propriedades.tsv"
NOME_PORTAO = "portao.tsv"
NOME_MANIFESTO = "manifesto.json"
NOME_ALERTAS_BULA = "alertas_bula.tsv"

# Criterios julgados por esta tarefa. Bula e expressao entram nas proximas; ate
# la saem `indeterminado` com o motivo escrito, e nao `passa`.
MOTIVO_SEM_BULA = "bula nao consultada (sem --termos-alerta)"


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


def _lista(valor: str | None) -> list[str]:
    """`--reprovar a,b` e `--reprovar a --reprovar b` dao o mesmo resultado."""
    return [p.strip() for p in (valor or "").split(",") if p.strip()]


def _veredito_estrutural(linha: dict):
    """(veredito, valor) do criterio `alerta_estrutural`.

    MARCA, nunca reprova por padrao: medido, PAINS pega a amitriptilina e
    BRENK pega 65% da base. O valor traz o catalogo e os filtros NOMEADOS --
    "BRENK = 1" nao se confere, "BRENK: nitro_group" se confere."""
    achados = []
    for coluna in propriedades.CATALOGOS:
        nomes = linha.get(coluna) or ""
        if nomes:
            achados.append(f"{coluna.removeprefix('alerta_')}: {nomes}")
    if not achados:
        return "passa", "nenhum catalogo casou"
    return "marca", "; ".join(achados)


PREFIXO_ADMET = "admet_"


def _mapa_de_colunas_admet(colunas_admet: list[str]) -> dict[str, str]:
    """{nome do ADMET-AI: nome publicado na tabela}.

    O ADMET-AI 2.0.1 publica uma coluna `tpsa`, que e o nome que este
    repositorio ja usa para a TPSA do RDKit. Sem tratamento, o cabecalho sai
    com `tpsa` DUAS vezes e o valor do RDKit e silenciosamente substituido pelo
    da outra ferramenta -- duas colunas respondendo a mesma pergunta sem
    procedencia, que e justamente o que `propriedades.py` proibe.

    A colidente e publicada com o prefixo `admet_`; as outras 103 mantem o nome
    ORIGINAL, que e o que as torna citaveis no artigo que cita a ferramenta. O
    manifesto registra cada renomeacao, para que ninguem precise adivinhar de
    onde veio `admet_tpsa`."""
    base = set(esquema.COLUNAS_PROPRIEDADES)
    originais = set(colunas_admet)
    mapa: dict[str, str] = {}
    for coluna in colunas_admet:
        if coluna not in base:
            mapa[coluna] = coluna
            continue
        publicado = PREFIXO_ADMET + coluna
        if publicado in base or publicado in originais:
            raise ValueError(
                f"coluna {coluna!r} do ADMET-AI colide com o esquema e o nome "
                f"de desvio {publicado!r} tambem ja existe; renomeie a coluna "
                "do esquema antes de seguir")
        mapa[coluna] = publicado
    return mapa


def _colunas_de_propriedades(mapa_admet: dict[str, str]) -> list[str]:
    """Identidade, descritores e alertas; depois as colunas do ADMET-AI;
    `procedencia` por ultimo.

    As 104 colunas nao estao fixadas no codigo: vem da instalacao. Sem o
    pacote nao ha como nomea-las, e inventa-las faria o repositorio mentir
    sobre a versao instalada -- entao o cabecalho simplesmente sai sem elas e o
    manifesto registra `admet_ai: ausente`."""
    base = [c for c in esquema.COLUNAS_PROPRIEDADES if c != "procedencia"]
    return base + list(mapa_admet.values()) + ["procedencia"]


def _argumentos(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--ligantes", required=True,
                    help="TSV da etapa 1 (nome, smiles, ...)")
    ap.add_argument("--saida", required=True, metavar="DIR",
                    help=f"diretorio que recebe {NOME_PROPRIEDADES}, "
                         f"{NOME_PORTAO} e {NOME_MANIFESTO}")
    ap.add_argument("--limite", type=_inteiro_positivo,
                    help="usar so os N primeiros ligantes")
    ap.add_argument("--percentil-toxicidade", type=float, default=90.0,
                    metavar="P",
                    help="percentil contra os aprovados do DrugBank acima do qual "
                         "o criterio `toxicidade` marca (padrao 90). Percentil, e "
                         "nao probabilidade crua: 0,67 de DILI nao quer dizer nada "
                         "sozinho -- a aspirina tem 0,67 e esta aprovada")
    ap.add_argument("--workers", type=int, default=0, metavar="N",
                    help="num_workers do ADMET-AI (padrao 0). Medido: 1.600 "
                         "moleculas saem em 8,6 s com 0 e em 384 s com 8, porque "
                         "sao dez modelos pagando o nascimento de N processos")
    ap.add_argument("--termos-alerta", action="append", default=[], metavar="T[,T]",
                    help="termos da condicao a buscar na bula (repetivel, aceita "
                         "virgula). SEM eles a checagem de bula NAO roda e o "
                         "criterio `alerta_bula` sai `indeterminado`")
    ap.add_argument("--limite-bula", type=_inteiro_positivo, default=500,
                    metavar="N",
                    help="teto de ligantes consultados no openFDA (padrao 500). "
                         "A base tem 108.935 ligantes e o openFDA e um servico "
                         "publico com limite de taxa; o que passar do teto sai "
                         "`indeterminado` com esse motivo escrito, nunca "
                         "'sem alerta'")
    ap.add_argument("--sem-rede", action="store_true",
                    help="nao consulta a bula; o criterio sai `indeterminado`")
    ap.add_argument("--cache", metavar="DIR",
                    help="raiz do cache de rede (padrao: cache_rede/ dentro de "
                         "--saida). Subpasta por versao dos servicos")
    ap.add_argument("--reprovar", action="append", default=[], metavar="C[,C]",
                    help="promove criterio(s) a reprovacao (repetivel). "
                         "Os que existem: " + ", ".join(decisao.CRITERIOS))
    ap.add_argument("--marcar", action="append", default=[], metavar="C[,C]",
                    help="rebaixa criterio(s) a somente marcar (repetivel)")
    return ap.parse_args(argv)


def _politica(args):
    promovidos = [c for v in args.reprovar for c in _lista(v)]
    rebaixados = [c for v in args.marcar for c in _lista(v)]
    repetidos = sorted(set(promovidos) & set(rebaixados))
    if repetidos:
        raise decisao.CriterioDesconhecido(
            f"criterio em --reprovar e --marcar ao mesmo tempo: "
            f"{', '.join(repetidos)}; escolha um")
    return decisao.montar_politica(promovidos, rebaixados)


def _abrir_cache(args, destino: Path):
    """Cache de rede da execucao: (Cache, versoes dos servicos, raiz, pasta).

    `com_openfda=True` acrescenta a data de carga do openFDA a chave da pasta:
    um cache gravado sem ela nao pode servir respostas de bula."""
    from mods.evidencias import cache as c
    raiz = Path(args.cache) if args.cache else destino / "cache_rede"
    cache, versoes = c.abrir_ao_vivo(raiz, com_openfda=True)
    return cache, versoes, str(Path(raiz).resolve()), cache.diretorio.name


def _consultar_bula(nome: str, termos: list[str], *, obter_json):
    """(veredito, valor, achados) do criterio `alerta_bula` para um ligante.

    A indicacao vem de consulta DIRIGIDA, nao da janela de 50 bulas: medido, a
    capsaicina tem 211 bulas e a unica que a indica para dor neuropatica esta
    fora das 50 devolvidas pela busca por nome. Sem a consulta dirigida ela e
    reprovada -- e ela TRATA a condicao, o que levaria o falso positivo do
    grupo de controle de 0/6 a 1/6."""
    resultado = bula.buscar_bulas(nome, obter_json=obter_json)
    achados = bula.achados(resultado["bulas"], termos)
    indicacao = bula.achados_de_indicacao(
        nome, resultado["campo_do_nome"], termos, obter_json=obter_json)
    # Deduplicado por (termo, secao): a consulta dirigida pode repetir um
    # achado que a janela ja trouxe.
    por_chave = {(a["termo"], a["secao"]): a for a in achados + indicacao}
    lista = [por_chave[c] for c in sorted(por_chave)]
    veredito, valor = bula.veredito(resultado, lista)
    return veredito, valor, lista, resultado


def _preparar_saida(args):
    """Diretorio de saida pronto, ou uma mensagem de erro.

    Falha ANTES de qualquer calculo: apontar `--saida` para um arquivo e erro
    de digitacao frequente, e descobri-lo depois de rodar o ADMET-AI na base
    inteira custaria a execucao."""
    destino = Path(args.saida).expanduser()
    if destino.exists() and not destino.is_dir():
        return None, (f"--saida deve ser um diretorio; {destino} e um arquivo "
                      "que existe e nao sera sobrescrito")
    if Path(args.ligantes).expanduser().resolve() in (
            destino / NOME_PROPRIEDADES, destino / NOME_PORTAO):
        return None, "--saida sobrescreveria o proprio arquivo de --ligantes"
    try:
        destino.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        return None, f"nao foi possivel criar {destino}: {e}"
    return destino, None


def _vereditos_do_ligante(linha_propriedades: dict, problema: str | None,
                          linha_admet: dict | None, percentil: float,
                          bula_veredito=None) -> dict:
    """{criterio: (veredito, valor)} para os criterios de LIGANTE.

    `expressao` nao entra: e criterio de ALVO, e `portao.tsv` tem uma linha por
    ligante. Ele vive em `alvos_anotados.tsv` (proxima tarefa)."""
    if problema:
        # Rede de seguranca: a leitura da etapa 1 ja descarta SMILES ilegivel,
        # com erro de linha e execucao parcial. Se um SMILES atravessar a
        # leitura e quebrar aqui, o ligante fica na saida reprovado, com o
        # motivo -- nunca desaparece em silencio.
        return {"smiles": ("reprova", problema),
                "alerta_bula": bula_veredito or ("indeterminado", MOTIVO_SEM_BULA),
                "alerta_estrutural": ("indeterminado", "sem molecula"),
                "toxicidade": ("indeterminado", "sem molecula")}
    return {
        "smiles": ("passa", ""),
        "alerta_bula": bula_veredito or ("indeterminado", MOTIVO_SEM_BULA),
        "alerta_estrutural": _veredito_estrutural(linha_propriedades),
        "toxicidade": admet.veredito_de_toxicidade(linha_admet, percentil),
    }


def _resumo(contador, politica, destino: Path) -> None:
    """O numero que o desenho exige: quantos ligantes CADA criterio reprovaria
    se fosse promovido. O botao e medido antes de ser apertado, em vez de o
    usuario promover um criterio e descobrir depois que perdeu a
    amitriptilina."""
    resumo = contador.resumo(politica)
    print(f"\n{resumo['itens']} ligante(s): {resumo['aprovados']} aprovado(s), "
          f"{resumo['reprovados']} reprovado(s). Por criterio:", file=sys.stderr)
    for criterio, dados in resumo["por_criterio"].items():
        adversos, regra = dados["adversos"], dados["politica"]
        efeito = ("reprovou" if regra == decisao.REPROVAR
                  else "reprovaria se promovido (--reprovar " + criterio + ")")
        print(f"  {criterio} [{regra}]: {adversos} {efeito}; "
              f"marca {dados['marca']}, reprova {dados['reprova']}, "
              f"passa {dados['passa']}, indeterminado {dados['indeterminado']}",
              file=sys.stderr)
    print(f"Saidas em {destino}", file=sys.stderr)


def main(argv: list[str] | None = None) -> int:
    args = _argumentos(argv)

    try:
        politica = _politica(args)
    except decisao.CriterioDesconhecido as e:
        print(f"[erro] {e}", file=sys.stderr)
        return SAIDA_ERRO

    destino, erro = _preparar_saida(args)
    if erro:
        print(f"[erro] {erro}", file=sys.stderr)
        return SAIDA_ERRO

    try:
        # `manter_ilegiveis`: um ligante que desaparece da tabela de vereditos
        # nao foi reprovado, foi apagado -- anti-objetivo declarado do desenho.
        # Ele fica, com as propriedades vazias e o criterio `smiles` reprovando.
        ligantes, problemas = ler_ligantes(args.ligantes, args.limite,
                                           manter_ilegiveis=True)
    except (OSError, UnicodeDecodeError) as e:
        print(f"[erro] nao foi possivel ler --ligantes: {e}", file=sys.stderr)
        return SAIDA_ERRO
    _relatar("ligantes", problemas)
    if not ligantes:
        print("[erro] ligantes: nenhum registro utilizavel; execucao "
              "interrompida", file=sys.stderr)
        return SAIDA_SEM_ENTRADA

    # ADMET-AI: tudo que pode falhar sem tocar nas saidas falha antes de abri-las.
    tem_admet = admet.disponivel()
    colunas_admet: list[str] = []
    versao_admet = "ausente"
    if tem_admet:
        try:
            colunas_admet = admet.colunas(args.workers)
            versao_admet = admet.versao()
        except (ImportError, OSError, RuntimeError, ValueError) as e:
            print(f"[erro] ADMET-AI instalado mas nao utilizavel: {e}",
                  file=sys.stderr)
            return SAIDA_ERRO
    else:
        print("[aviso] pacote `admet_ai` ausente: a etapa roda sem ADMET. As "
              "colunas de ADMET nao saem e o criterio `toxicidade` fica "
              "`indeterminado` -- que NUNCA reprova. Instale com: "
              "pip install admet-ai", file=sys.stderr)

    termos = [t for v in args.termos_alerta for t in _lista(v)]
    if termos and args.sem_rede:
        print("[erro] --termos-alerta precisa de rede (a bula vem do openFDA); "
              "incompativel com --sem-rede", file=sys.stderr)
        return SAIDA_ERRO
    # Cache e versoes ANTES de abrir qualquer saida: sem versao conhecida nao
    # ha como garantir que o cache vale, e falhar aqui nao deixa arquivo algum.
    cache = None
    versoes_rede: dict = {}
    raiz_cache = pasta_cache = None
    if termos:
        try:
            cache, versoes_rede, raiz_cache, pasta_cache = _abrir_cache(args, destino)
        except Exception as e:   # VersaoIndeterminada, OSError, ValueError
            print(f"[erro] cache de rede: {e}", file=sys.stderr)
            return SAIDA_ERRO

    contador = decisao.Contador(decisao.CRITERIOS_LIGANTE)
    consultados = 0
    falhas_de_bula: list[str] = []
    try:
        mapa_admet = _mapa_de_colunas_admet(colunas_admet)
    except ValueError as e:
        print(f"[erro] colunas do ADMET-AI: {e}", file=sys.stderr)
        return SAIDA_ERRO
    renomeadas = {o: p for o, p in mapa_admet.items() if o != p}
    if renomeadas:
        print("[info] coluna(s) do ADMET-AI publicada(s) com prefixo para nao "
              "colidir com o esquema (o manifesto registra): "
              + ", ".join(f"{o} -> {p}" for o, p in sorted(renomeadas.items())),
              file=sys.stderr)
    colunas_propriedades = _colunas_de_propriedades(mapa_admet)
    colunas_portao = esquema.colunas_portao(decisao.CRITERIOS_LIGANTE)
    procedencia_base = f"RDKit {propriedades.versao_rdkit()}"
    if tem_admet:
        procedencia_base += f"; ADMET-AI {versao_admet}"

    pendentes: list[tuple[Path, Path]] = []
    publicados: list[Path] = []
    etapa = "propriedades"
    try:
        tmp_propriedades = atomico.temporario(destino / NOME_PROPRIEDADES)
        pendentes.append((tmp_propriedades, destino / NOME_PROPRIEDADES))
        tmp_portao = atomico.temporario(destino / NOME_PORTAO)
        pendentes.append((tmp_portao, destino / NOME_PORTAO))
        if termos:
            tmp_bula = atomico.temporario(destino / NOME_ALERTAS_BULA)
            pendentes.append((tmp_bula, destino / NOME_ALERTAS_BULA))
        else:
            # Tabela vazia seria pior que ausente: "sem alerta" e "nao
            # perguntei" nao sao a mesma afirmacao.
            tmp_bula = None

        from contextlib import ExitStack
        with ExitStack() as pilha:
            esc_prop = pilha.enter_context(
                esquema.EscritorTSV(tmp_propriedades, colunas_propriedades))
            esc_portao = pilha.enter_context(
                esquema.EscritorTSV(tmp_portao, colunas_portao))
            esc_bula = pilha.enter_context(esquema.EscritorTSV(
                tmp_bula, esquema.COLUNAS_ALERTAS_BULA)) if tmp_bula else None
            procedencia_bula = f"openFDA {versoes_rede.get('openfda', '')}".strip()
            # Em lotes porque a base real tem 108.935 ligantes e
            # `propriedades.tsv` tem 120 colunas: so o lote corrente fica em
            # memoria. A ordem da entrada e preservada -- a repeticao exige
            # saida deterministica.
            for inicio in range(0, len(ligantes), admet.LOTE):
                lote = ligantes[inicio:inicio + admet.LOTE]
                previsoes: dict[str, dict] = {}
                if tem_admet:
                    previsoes = admet.prever(
                        [l.get("smiles", "") for l in lote],
                        num_workers=args.workers)
                for ligante in lote:
                    linha, problema = propriedades.linha_de_propriedades(ligante)
                    # SMILES ausente do resultado e SMILES que o ADMET-AI
                    # descartou; nao se inventa valor para ele.
                    linha_admet = previsoes.get(ligante.get("smiles", ""))
                    if linha_admet:
                        linha.update({publicado: linha_admet.get(original)
                                      for original, publicado
                                      in mapa_admet.items()})
                    linha["procedencia"] = procedencia_base
                    esc_prop.escrever(linha)

                    bula_veredito = None
                    if termos:
                        if consultados >= args.limite_bula:
                            bula_veredito = ("indeterminado",
                                             f"acima do teto de bula "
                                             f"(--limite-bula {args.limite_bula})")
                        else:
                            consultados += 1
                            try:
                                veredito_b, valor_b, achados_b, resultado_b = (
                                    _consultar_bula(linha["ligante"], termos,
                                                    obter_json=cache.obter_json))
                            except Exception as e:
                                # Falha de UM item nao derruba a execucao; ela
                                # fica PARCIAL e o manifesto diz qual item
                                # faltou. Rede caida nunca vira "sem alerta".
                                falhas_de_bula.append(
                                    f"{linha['ligante']}: {type(e).__name__}: {e}")
                                bula_veredito = ("indeterminado",
                                                 f"consulta a bula falhou: {e}")
                            else:
                                bula_veredito = (veredito_b, valor_b)
                                for achado in achados_b:
                                    esc_bula.escrever({
                                        "ligante": linha["ligante"],
                                        "ligante_inchikey_parental": linha[
                                            "ligante_inchikey_parental"],
                                        "n_bulas_consultadas": len(
                                            resultado_b["bulas"]),
                                        "campo_do_nome": resultado_b[
                                            "campo_do_nome"],
                                        "procedencia": procedencia_bula,
                                        **achado})
                    vereditos = _vereditos_do_ligante(
                        linha, problema, linha_admet, args.percentil_toxicidade,
                        bula_veredito)
                    identidade = {
                        "ligante": linha["ligante"],
                        "ligante_inchikey_parental": linha[
                            "ligante_inchikey_parental"]}
                    linha_portao = decisao.linha_de_portao(
                        identidade, vereditos, politica)
                    esc_portao.escrever(linha_portao)
                    contador.contar(vereditos, linha_portao["aprovado"])

        etapa = "manifesto"
        descartadas = sum(p.gravidade == "erro" for p in problemas)
        ilegiveis = sum(1 for l in ligantes if not l.get("inchikey_parental"))
        parcial = bool(descartadas) or bool(falhas_de_bula)
        resumo = contador.resumo(politica)
        tmp_manifesto = atomico.temporario(destino / NOME_MANIFESTO)
        pendentes.append((tmp_manifesto, destino / NOME_MANIFESTO))
        tmp_manifesto.write_text(texto_do_manifesto(
            parametros={"ligantes": args.ligantes,
                        "politica": politica,
                        "percentil_toxicidade": args.percentil_toxicidade,
                        "workers": args.workers,
                        "limite": args.limite,
                        "termos_alerta": termos,
                        "limite_bula": args.limite_bula,
                        "sem_rede": args.sem_rede,
                        "cache_raiz": raiz_cache,
                        "cache_pasta": pasta_cache},
            versoes={**versoes_rede, "admet_ai": versao_admet,
                     "ligantes_md5": procedencia.md5_ou_none(args.ligantes),
                     **procedencia.versoes_do_ambiente(),
                     **procedencia.repositorio_do_pipeline(),
                     "pipeline_commit": procedencia.commit_do_pipeline(),
                     "data": datetime.now(timezone.utc).isoformat(
                         timespec="seconds")},
            contagens={"ligantes": len(ligantes),
                       "problemas_ligantes": len(problemas),
                       "linhas_descartadas": descartadas,
                       "smiles_ilegiveis_mantidos": ilegiveis,
                       "bulas_consultadas": consultados,
                       "falhas_de_bula": falhas_de_bula,
                       "execucao_parcial": parcial,
                       "colunas_de_admet": len(mapa_admet),
                       "colunas_de_admet_renomeadas": renomeadas,
                       **resumo}), encoding="utf-8")

        # So aqui algum arquivo do usuario e substituido, e todos de uma vez.
        etapa = "publicacao"
        for arquivo, caminho in pendentes:
            atomico.publicar(arquivo, caminho)
            publicados.append(caminho)
    except (OSError, RuntimeError, ValueError, ImportError) as e:
        if etapa == "publicacao":
            faltaram = [str(d) for _, d in pendentes if d not in publicados]
            print(f"[erro] publicacao incompleta: {e}. Publicado(s): "
                  f"{', '.join(str(d) for d in publicados) or 'nenhum'}; NAO "
                  f"publicado(s): {', '.join(faltaram)}", file=sys.stderr)
        else:
            print(f"[erro] {etapa}: {e}", file=sys.stderr)
        return SAIDA_ERRO
    finally:
        for arquivo, _ in pendentes:
            arquivo.unlink(missing_ok=True)

    for falha in falhas_de_bula:
        print(f"[aviso] bula, {falha}", file=sys.stderr)
    if ilegiveis:
        # Codigo 0, nao 3: o portao DECIDIU sobre esses ligantes (reprovados
        # pelo criterio `smiles`, com o motivo na linha). Execucao parcial e
        # quando algo se perde, e aqui nada se perdeu.
        print(f"[aviso] {ilegiveis} ligante(s) com SMILES que o RDKit nao le: "
              "estao na saida com as propriedades vazias e reprovados pelo "
              "criterio `smiles` (veja a coluna `motivos`), nao apagados",
              file=sys.stderr)
    _resumo(contador, politica, destino)
    if parcial:
        print(f"[ATENCAO] execucao PARCIAL (codigo {SAIDA_PARCIAL}): "
              f"{descartadas} linha(s) de entrada descartada(s) por erro e "
              f"{len(falhas_de_bula)} falha(s) de consulta a bula (detalhes "
              "acima). As saidas existem, mas estao incompletas; o manifesto "
              "registra isso em execucao_parcial", file=sys.stderr)
    return SAIDA_PARCIAL if parcial else SAIDA_OK


if __name__ == "__main__":
    sys.exit(main())
