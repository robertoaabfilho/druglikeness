#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ranquear.py

Fase 3: ordena compostos pela evidencia de alvo da etapa 1, em faixas
lexicograficas, e publica `ranking.tsv`, `nao_ranqueados.tsv` e o
`manifesto.json` dentro de `--saida`.

Desenho: `docs/superpowers/specs/2026-11-26-fase3-ranqueamento-design.md`.

Por que a ordem nao sai de semelhanca quimica
---------------------------------------------
Ordenar por semelhanca fisico-quimica com as drogas de referencia do sintoma e
a abordagem obvia -- e e cega a alvo. Medido no controle de 12 farmacos: os dois QED mais altos sao ciprofloxacino (0,893) e
linezolida (0,892), AMBOS do grupo que CAUSA dor neuropatica, enquanto os que
tratam ficam entre 0,539 e 0,814. Nenhum ajuste de peso corrige isso, porque a
informacao que separa os dois grupos -- em que alvo o composto age -- nao
entrava na conta.

As faixas, e por que nao ha constante nenhuma aqui
--------------------------------------------------
    1a  MEDIDA em alvo do usuario          -> max pchembl
    1b  MEDIDA em alvo do enriquecimento   -> max pchembl
    2a  so PREVISTA, alvo do usuario       -> max p x (1 - p_negativo)
    2b  so PREVISTA, alvo do enriquecimento
    3   nenhuma evidencia afirmativa       -> NAO ranqueia

A procedencia e sub-faixa em vez de peso, e a classe de evidencia e a chave
externa: evidencia MEDIDA num alvo que o pesquisador nao listou e uma
descoberta, e rebaixa-la sob a previsao do alvo ja conhecido esconderia o unico
resultado que o enriquecimento tem para dar.

Codigos de saida
----------------
0 completa; 3 PARCIAL; 1 nenhum composto utilizavel; 2 erro (nada publicado).
"""

import argparse
import sys
from pathlib import Path
from typing import Callable, NamedTuple

# Permite `python mods/ranquear.py`: a raiz do repositorio entra no caminho.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from mods.evidencias import atomico, procedencia
from mods.evidencias.manifesto import (hashes_de_entrada_ausentes,
                                       ler_manifesto, texto_do_manifesto,
                                       verificar_manifesto)
from mods.evidencias.resumo import MINIMO_TREINO
from mods.portao.esquema import EscritorTSV, ler_tsv
from mods.ranque import esquema, negativo
from mods.ranque.faixa import (ORIGENS_DO_USUARIO, chave, classificar,
                               ORIGENS_DE_ENRIQUECIMENTO)

SAIDA_OK = 0
SAIDA_SEM_ENTRADA = 1
SAIDA_ERRO = 2
SAIDA_PARCIAL = 3

NOME_RANKING = "ranking.tsv"
NOME_NAO_RANQUEADOS = "nao_ranqueados.tsv"
NOME_MANIFESTO = "manifesto.json"


INDETERMINADO = "indeterminado"
LOTE = 512
NOME_MANIFESTO_REPETICAO = "manifesto_repeticao.json"

# Chaves que, se mudarem, fazem a repeticao dar outro arquivo. As da etapa 1
# nao servem: esta fase nao consulta servico nenhum, e depende do PP4CKD, que
# a etapa 2 nao conhece.
VERSOES_DA_REPETICAO = ("rdkit", "torch", "pp4ckd_md5", "rotulos_md5",
                        "evidencias_md5", "portao_md5", "alvos_md5",
                        "ligantes_md5")
# `rotulos_md5` entra porque o `target_labels.tsv` decide organismo, `n_ativos`
# e `pref_name` -- a mascara do canal negativo inteira. Sem ele, repetir com
# outro espaco de alvos era ACEITO e publicava outro escore como se fosse o
# mesmo. A etapa 1 ja listava `rotulos_md5` entre as chaves dela.
# `dispositivo` NAO entra: ele so e conhecido depois de abrir o modelo, e a
# verificacao tem de recusar antes de escrever qualquer coisa. Divergencia de
# dispositivo AVISA, como na etapa 2, onde repetir uma execucao de GPU avisa em
# vez de prometer igualdade bit a bit.
PARAMETROS_DA_REPETICAO = ("canal_negativo", "minimo_treino", "lote",
                           "limite", "portao")
# MD5 que o manifesto TEM de trazer para a repeticao poder conferir a entrada.
ENTRADAS_HASHEADAS = ("evidencias_md5",)
# As condicionais sao lidas dos PARAMETROS do manifesto anterior (e o que
# `hashes_de_entrada_ausentes` compara), por isso a chave e o nome do
# parametro registrado, nao o da flag.
CONDICIONAIS = (("portao", "portao_md5"), ("alvos", "alvos_md5"),
                ("canal_negativo", "pp4ckd_md5"),
                ("canal_negativo", "rotulos_md5"),
                ("canal_negativo", "ligantes_md5"))


class Modelo(NamedTuple):
    """O PP4CKD como esta fase o usa, com a costura que os testes substituem.

    `prever` recebe uma lista de fingerprints e devolve a matriz de
    probabilidades (um composto por linha, um alvo por coluna). A rede fica
    atras dela para que o canal possa ser testado sem o checkpoint real, que
    nao esta no repositorio."""
    rotulos: list[dict]
    cfg: dict
    dispositivo: str
    prever: Callable
    checkpoint: str


def _inteiro_positivo(texto: str) -> int:
    valor = int(texto)
    if valor <= 0:
        raise argparse.ArgumentTypeError("tem de ser maior que zero")
    return valor


def analisar(argv=None):
    ap = argparse.ArgumentParser(
        prog="ranquear.py",
        description="Fase 3: ordena compostos pela evidencia de alvo.")
    ap.add_argument("--evidencias", required=True, metavar="TSV",
                    help="evidencias.tsv da etapa 1")
    ap.add_argument("--saida", required=True, metavar="DIR",
                    help="DIRETORIO que recebe as tabelas e o manifesto")
    ap.add_argument("--portao", metavar="TSV",
                    help="portao.tsv da etapa 2. Reprovado nao ranqueia, mas "
                         "sai em " + NOME_NAO_RANQUEADOS + " com o criterio. "
                         "Sem isto, a coluna diz 'nao consultado' -- nunca "
                         "'aprovado' por omissao")
    ap.add_argument("--alvos", metavar="TSV",
                    help="alvos_enriquecidos.tsv da etapa 1")
    ap.add_argument("--limite", type=_inteiro_positivo, metavar="N",
                    help="usar so os N primeiros compostos")
    ap.add_argument("--ligantes", metavar="TSV",
                    help="ligantes.tsv da etapa 1. O canal negativo precisa do "
                         "SMILES, que `evidencias.tsv` nao traz -- sem isto o "
                         "canal nao roda")
    ap.add_argument("--checkpoint", metavar="PT",
                    help="pesos do PP4CKD. Sem eles nao ha canal negativo, e o "
                         "manifesto registra `canal_negativo: false`")
    ap.add_argument("--rotulos", metavar="TSV",
                    help="target_labels.tsv do PP4CKD: diz quais alvos existem "
                         "e quantas moleculas treinaram cada um")
    ap.add_argument("--manifesto-etapa1", metavar="JSON",
                    help="manifesto.json da etapa 1. RECUSA rodar se o MD5 do "
                         "--checkpoint diferir do que ele registra: os dois "
                         "lados de `p x (1 - p_negativo)` tem de vir do mesmo "
                         "modelo, senao a ordem depende de quais arquivos "
                         "estavam em disco")
    ap.add_argument("--minimo-treino", type=_inteiro_positivo,
                    default=MINIMO_TREINO, metavar="N",
                    help=f"alvo com menos de N moleculas ativas nao entra no "
                         f"canal negativo (padrao {MINIMO_TREINO}). Requisito "
                         "1: esparso e ausencia de evidencia, nao evidencia de "
                         "atividade")
    ap.add_argument("--lote", type=_inteiro_positivo, default=LOTE,
                    metavar="N",
                    help=f"compostos por lote no forward (padrao {LOTE}). A "
                         "matriz completa seriam 3,3 GB; so o agregado por "
                         "composto e guardado")
    ap.add_argument("--gpu", action="store_true",
                    help="roda o PP4CKD na GPU. O padrao e CPU porque a GPU "
                         "nao reproduz -- medido na etapa 2: 36 de 312 colunas "
                         "diferem entre duas execucoes de GPU, e 0 de 312 na CPU")
    ap.add_argument("--manifesto", metavar="JSON",
                    help="MODO DE REPETICAO: manifesto.json de uma execucao "
                         "anterior. RECUSA rodar se um MD5 de entrada, a "
                         "versao do RDKit ou do PP4CKD, ou um parametro que "
                         "muda a saida diferirem do registrado. O manifesto "
                         f"novo sai como {NOME_MANIFESTO_REPETICAO} e nunca "
                         "sobrescreve o original")
    ap.add_argument("--sem-canal-negativo", action="store_true",
                    help="desliga o canal negativo; a faixa prevista ordena "
                         "pelo p sozinho. Medido, o canal leva a AUC de 0,750 "
                         "para 0,806 no controle de 12 -- mas n = 12 e os seis "
                         "causadores sao todos antibioticos, entao a regra e "
                         "desligavel de proposito")
    return ap.parse_args(argv)


def _erro_de_checkpoint(checkpoint, manifesto_etapa1) -> str | None:
    """Mensagem quando o checkpoint nao e o que a etapa 1 usou, ou None.

    `p_afirmativo` vem de `evidencias.tsv`, publicado pela etapa 1; o
    `p_negativo` vem do forward desta fase. Se sairem de checkpoints
    diferentes, a multiplicacao combina dois modelos num numero so, e a ordem
    passa a depender de qual par de arquivos estava em disco."""
    if not (checkpoint and manifesto_etapa1):
        return None
    # A existencia do checkpoint e conferida AQUI, e nao na varredura geral de
    # arquivos, porque so este guarda precisa ler o arquivo. Sem o guarda, quem
    # reporta um checkpoint ausente e o carregador, com a mensagem dele.
    if not Path(checkpoint).is_file():
        return f"--checkpoint: arquivo nao encontrado: {checkpoint}"
    try:
        registrado = (ler_manifesto(manifesto_etapa1)
                      .get("versoes", {}).get("checkpoint_md5"))
    except (OSError, ValueError) as e:
        return f"--manifesto-etapa1 nao pode ser lido: {e}"
    if not registrado:
        return ("--manifesto-etapa1 nao registra `checkpoint_md5`; so "
                "manifestos de execucoes que usaram os pesos podem servir de "
                "referencia para o canal negativo")
    atual = procedencia.md5(str(checkpoint))
    if atual != registrado:
        return ("o --checkpoint nao e o que a etapa 1 usou: md5 atual "
                f"{atual}, md5 registrado {registrado}. Os dois lados de "
                "`p x (1 - p_negativo)` tem de vir do mesmo modelo")
    return None


def _validar(args) -> str | None:
    """Mensagem de erro, ou None. Roda ANTES de abrir qualquer saida."""
    quer_canal = bool(args.checkpoint) and not args.sem_canal_negativo
    if quer_canal and not args.ligantes:
        return ("--checkpoint pede o canal negativo, que precisa do SMILES: "
                "informe --ligantes (o evidencias.tsv nao traz essa coluna)")
    if quer_canal and not args.rotulos:
        return "--checkpoint sem --rotulos: nao se sabe quais alvos existem"
    if quer_canal and not args.manifesto_etapa1:
        return ("--checkpoint sem --manifesto-etapa1: sem ele nao ha como "
                "conferir que `p_afirmativo` (da etapa 1) e `p_negativo` "
                "(desta fase) vem do MESMO modelo, e a ordem passaria a "
                "depender de quais arquivos estavam em disco")
    erro = _erro_de_checkpoint(args.checkpoint, args.manifesto_etapa1)
    if erro:
        return erro
    for rotulo, caminho in (("--evidencias", args.evidencias),
                            ("--portao", args.portao),
                            ("--ligantes", args.ligantes),
                            ("--rotulos", args.rotulos),
                            ("--manifesto-etapa1", args.manifesto_etapa1),
                            ("--manifesto", args.manifesto),
                            ("--alvos", args.alvos)):
        if caminho and not Path(caminho).is_file():
            return f"{rotulo}: arquivo nao encontrado: {caminho}"
    destino = Path(args.saida)
    if destino.exists() and not destino.is_dir():
        return (f"--saida tem de ser um diretorio; {args.saida} e um arquivo "
                "que existe")
    # Por NOME de saida, e nao por diretorio: o manifesto da etapa 1 chama-se
    # `manifesto.json`, que e exatamente um dos arquivos que esta fase publica,
    # e ele e a chave de repeticao DELA.
    publicados = {(destino / n).resolve()
                  for n in (NOME_RANKING, NOME_NAO_RANQUEADOS,
                            NOME_MANIFESTO_REPETICAO if args.manifesto
                            else NOME_MANIFESTO)}
    for rotulo, caminho in (("--evidencias", args.evidencias),
                            ("--portao", args.portao),
                            ("--alvos", args.alvos),
                            ("--ligantes", args.ligantes),
                            ("--manifesto-etapa1", args.manifesto_etapa1),
                            ("--manifesto", args.manifesto)):
        if caminho and Path(caminho).resolve() in publicados:
            return (f"{rotulo} aponta para um arquivo que esta execucao "
                    f"publicaria ({Path(caminho).name}); a entrada seria "
                    "sobrescrita")
    return None


def _por_composto(caminho, limite=None):
    """({inchikey: (nome, [linhas])}, linhas descartadas).

    A contagem de descartadas e devolvida, nao impressa: descartar linha em
    silencio com codigo 0 e o defeito que o codigo 3 existe para evitar --
    automacao nao le `stderr`."""
    compostos: dict[str, tuple[str, list[dict]]] = {}
    descartadas = 0
    for linha in ler_tsv(caminho):
        chave_composto = (linha.get("ligante_inchikey_parental") or "").strip()
        nome = (linha.get("ligante") or "").strip()
        if not chave_composto:
            # Sem a forma parental nao ha chave de juncao: a linha nao tem a
            # quem pertencer.
            descartadas += 1
            continue
        if chave_composto not in compostos:
            if limite is not None and len(compostos) >= limite:
                continue
            compostos[chave_composto] = (nome, [])
        compostos[chave_composto][1].append(linha)
    return compostos, descartadas


def _portao(caminho):
    """{inchikey: (aprovado, motivos)} de `portao.tsv`, ou None sem o arquivo."""
    if not caminho:
        return None
    por_composto = {}
    for linha in ler_tsv(caminho):
        chave_composto = (linha.get("ligante_inchikey_parental") or "").strip()
        if chave_composto:
            por_composto[chave_composto] = (
                (linha.get("aprovado") or "").strip(),
                (linha.get("motivos") or "").strip())
    return por_composto


def _versao_do_torch(args) -> str | None:
    """A versao do torch, ou None quando esta execucao nao usa o modelo.

    Ela muda o forward, logo muda a saida, e `versoes_do_ambiente()` registra
    so python e rdkit. Importada aqui dentro para que quem roda com
    `--sem-canal-negativo` nao pague o import."""
    if args.sem_canal_negativo or not args.checkpoint:
        return None
    try:
        import torch
        return torch.__version__
    except ImportError:
        return None


def _erro_de_largura(saidas_do_checkpoint: int, n_rotulos: int) -> str | None:
    """Mensagem quando a largura do checkpoint nao e o numero de rotulos.

    A coluna *i* de um nao e o alvo *i* do outro: pontuar assim atribuiria a
    probabilidade de um alvo ao nome de outro, e o canal negativo mascararia
    as colunas erradas. A etapa 1 ja aplica este guarda
    (`predito.py`); o tarball do Zenodo esta truncado na origem, entao o par
    divergente e o caso ESPERADO, nao a excecao."""
    if saidas_do_checkpoint == n_rotulos:
        return None
    return (f"o --checkpoint tem {saidas_do_checkpoint} saidas e os --rotulos "
            f"tem {n_rotulos} alvos: a coluna i de um nao e o alvo i do outro. "
            "Troque os dois juntos")


def _parametros(args, canal_negativo, dispositivo) -> dict:
    """Os parametros que mudam o CONTEUDO da saida. Caminhos de saida nao
    entram: repetir para outra pasta e repetir."""
    return {"canal_negativo": canal_negativo,
            "dispositivo": dispositivo,
            "alvos": bool(args.alvos),
            "minimo_treino": args.minimo_treino,
            "lote": args.lote,
            "limite": args.limite,
            "portao": bool(args.portao)}


def _motivo_do_portao(aprovado: str, motivos: str) -> list[str]:
    """O motivo da reprovacao do portao, ou lista vazia.

    Uma funcao, e nao duas escritas no lugar, porque o motivo nao pode depender
    do CAMINHO pelo qual o composto saiu: quem nao tem evidencia e foi
    reprovado tem de dizer as duas coisas, igual a quem tem evidencia em alvo
    nao afirmativo."""
    if aprovado != "nao":
        return []
    return [f"reprovado no portao ({motivos})" if motivos
            else "reprovado no portao"]


def _recusar(problemas: list[str]) -> None:
    print("[ERRO] a repeticao foi recusada; nada foi escrito:",
          file=sys.stderr)
    for problema in problemas:
        print(f"  - {problema}", file=sys.stderr)


def _abrir_modelo(args) -> Modelo | None:
    """O PP4CKD, ou None quando esta execucao nao tem canal negativo.

    O import de `torch` fica aqui dentro: uma execucao com
    `--sem-canal-negativo` nao paga por ele, e os testes substituem esta
    funcao em vez de montar um checkpoint de 27 MB."""
    if args.sem_canal_negativo or not args.checkpoint:
        return None
    import torch

    from mods.modelo import carregar_modelo, carregar_rotulos

    rotulos = carregar_rotulos(args.rotulos)
    dispositivo = "cuda" if args.gpu else "cpu"
    rede, cfg, dispositivo, nome = carregar_modelo(
        None, "ECFP4", dispositivo, checkpoint=args.checkpoint)
    saidas = next((m.out_features for m in reversed(list(rede.modules()))
                   if hasattr(m, "out_features")), len(rotulos))
    erro = _erro_de_largura(saidas, len(rotulos))
    if erro:
        raise ValueError(erro)

    def prever(fingerprints):
        import numpy as np
        lote = torch.tensor(np.array(fingerprints), dtype=torch.float32)
        with torch.no_grad():
            return torch.sigmoid(rede(lote.to(dispositivo))).cpu().numpy()

    return Modelo(rotulos=rotulos, cfg=cfg, dispositivo=dispositivo,
                  prever=prever, checkpoint=nome)


def _ligantes_de_entrada(caminho) -> dict[str, tuple[str, str]]:
    """{inchikey parental: (nome, smiles)} de `ligantes.tsv`.

    Serve a dois propositos: o SMILES alimenta o canal negativo, e a lista de
    chaves e a unica forma de saber quem a etapa 1 conhecia e `evidencias.tsv`
    nao menciona."""
    if not caminho:
        return {}
    por_composto: dict[str, tuple[str, str]] = {}
    for linha in ler_tsv(caminho):
        chave_composto = (linha.get("inchikey_parental")
                          or linha.get("ligante_inchikey_parental") or "").strip()
        if chave_composto and chave_composto not in por_composto:
            por_composto[chave_composto] = (
                (linha.get("nome") or linha.get("ligante") or "").strip(),
                (linha.get("smiles") or "").strip())
    return por_composto


def _alvos_afirmativos(args) -> tuple[set[str], set[str]]:
    """(IDs, nomes) dos alvos AFIRMATIVOS desta execucao.

    Lidos da LISTA DECLARADA -- `--alvos` quando ha, e o `evidencias.tsv`
    INTEIRO, sem o corte de `--limite`.

    Antes isto saia do conjunto de compostos ja truncado, e o efeito era a
    inversao que o requisito 2 proibe: o alvo que o usuario PEDIU virava
    antialvo assim que nenhum composto incluido tivesse linha nele, e a
    atividade no alvo pretendido passava a penalizar o composto. O
    `evidencias.tsv` tambem ja vem filtrado pelo `--corte` da etapa 1, entao
    deduzir a declaracao do subproduto nunca foi seguro.
    """
    afirmativas = set(ORIGENS_DO_USUARIO) | set(ORIGENS_DE_ENRIQUECIMENTO)
    ids: set[str] = set()
    nomes: set[str] = set()

    def juntar(linha, coluna_origem):
        if (linha.get(coluna_origem) or "").strip() not in afirmativas:
            return
        id_alvo = (linha.get("alvo_chembl_id") or "").strip().upper()
        if id_alvo:
            ids.add(id_alvo)
        nome = (linha.get("alvo") or "").strip().casefold()
        if nome:
            nomes.add(nome)

    if args.alvos:
        for linha in ler_tsv(args.alvos):
            # A lista enriquecida chama a coluna `origem`; a de evidencias,
            # `origem_alvo`. Aceitar as duas e o que `entrada.py` ja faz.
            juntar(linha, "origem" if "origem" in linha else "origem_alvo")
    for linha in ler_tsv(args.evidencias):
        juntar(linha, "origem_alvo")
    return ids, nomes


def _indices_afirmativos(ids: set[str], nomes: set[str], rotulos) -> list[int]:
    """Indices, no modelo, dos alvos afirmativos."""
    indices = []
    for i, r in enumerate(rotulos):
        if (r.get("chembl_id") or "").strip().upper() in ids:
            indices.append(i)
        elif (r.get("nome") or "").strip().casefold() in nomes:
            indices.append(i)
    return indices


def _canal_por_composto(modelo, compostos, smiles_por_composto,
                        afirmativos: list[int], minimo_treino: int, lote: int):
    """{inchikey: (p_negativo ou None, nome do alvo)}.

    So o agregado por composto e guardado: a matriz inteira seriam
    108.935 x 7.676 (~3,3 GB em float32), e ela nao e necessaria."""
    from mods.modelo import fingerprint_de

    mascara = negativo.mascara(modelo.rotulos, afirmativos,
                               minimo_treino=minimo_treino)
    resultado: dict[str, tuple[float | None, str]] = {}
    sem_fingerprint = 0
    pendentes: list[tuple[str, object]] = []

    def despachar():
        if not pendentes:
            return
        matriz = modelo.prever([fp for _, fp in pendentes])
        for (inchikey, _), linha in zip(pendentes, matriz):
            resultado[inchikey] = (
                negativo.agregar(linha, mascara),
                negativo.alvo_do_maximo(modelo.rotulos, linha, mascara))
        pendentes.clear()

    for inchikey in compostos:
        fp = fingerprint_de(smiles_por_composto.get(inchikey, ""), modelo.cfg)
        if fp is None:
            # Sem fingerprint nao ha forward. `indeterminado`, nunca zero:
            # zero afirmaria "nenhuma atividade indesejada". E conta como falha
            # de item: o canal foi PEDIDO e nao aconteceu para este composto.
            resultado[inchikey] = (None, "")
            sem_fingerprint += 1
            continue
        pendentes.append((inchikey, fp))
        if len(pendentes) >= lote:
            despachar()
    despachar()
    return resultado, sem_fingerprint


def _veredito(portao, inchikey) -> tuple[str, str]:
    if portao is None:
        return esquema.PORTAO_NAO_CONSULTADO, ""
    if inchikey not in portao:
        return esquema.AUSENTE_DO_PORTAO, ""
    return portao[inchikey]


def _texto(valor) -> str:
    """Numero numa celula, sem zeros a direita inuteis.

    Seis casas FIXAS achatavam para `0` tudo abaixo de 5e-7, e varios
    compostos ordenados entre si mostravam o mesmo `0` -- a coluna que explica
    a ordem deixava de explica-la. Seis digitos SIGNIFICATIVOS preservam a
    magnitude; a notacao exponencial aparece so quando ela e necessaria para
    nao mentir sobre o valor."""
    if valor is None or valor == "":
        return ""
    numero = float(valor)
    if numero == 0.0:
        return "0"
    texto = f"{numero:.6g}"
    return texto.rstrip("0").rstrip(".") if "." in texto and "e" not in texto \
        else texto


def main(argv=None) -> int:
    args = analisar(argv)
    erro = _validar(args)
    if erro:
        print(f"[ERRO] {erro}", file=sys.stderr)
        return SAIDA_ERRO

    try:
        compostos, linhas_descartadas = _por_composto(
            args.evidencias, args.limite)
        portao = _portao(args.portao)
    except (OSError, ValueError) as e:
        print(f"[ERRO] {e}", file=sys.stderr)
        return SAIDA_ERRO

    if not compostos:
        print("[ERRO] nenhum composto utilizavel em " + args.evidencias,
              file=sys.stderr)
        return SAIDA_SEM_ENTRADA

    versoes = procedencia.versoes_do_ambiente() | {
        "pipeline_commit": procedencia.commit_do_pipeline(),
        "evidencias_md5": procedencia.md5(args.evidencias),
        "portao_md5": procedencia.md5_ou_none(args.portao),
        "alvos_md5": procedencia.md5_ou_none(args.alvos),
        "ligantes_md5": procedencia.md5_ou_none(args.ligantes),
        "rotulos_md5": procedencia.md5_ou_none(args.rotulos),
        "pp4ckd_md5": procedencia.md5_ou_none(args.checkpoint),
        "torch": _versao_do_torch(args)}

    anterior = None
    if args.manifesto:
        try:
            anterior = ler_manifesto(args.manifesto)
        except (OSError, ValueError) as e:
            print(f"[ERRO] --manifesto nao pode ser lido: {e}",
                  file=sys.stderr)
            return SAIDA_ERRO
        # `verificar_manifesto` so compara chaves PRESENTES no anterior, entao
        # um manifesto truncado ou editado a mao "repetiria" sem conferir nada.
        ausentes = hashes_de_entrada_ausentes(
            anterior, entradas=ENTRADAS_HASHEADAS, condicionais=CONDICIONAIS)
        if ausentes:
            _recusar([f"o manifesto nao registra {c}; ele nao pode servir de "
                      "referencia para a repeticao" for c in ausentes])
            return SAIDA_ERRO
        # Primeira passada: tudo o que NAO depende do modelo. Assim um MD5 que
        # mudou recusa antes de pagar o carregamento dos pesos -- a licao da
        # etapa 2, onde a validacao rodava depois de carregar dez modelos.
        problemas = verificar_manifesto(
            anterior, {"versoes": versoes,
                       "parametros": _parametros(args, None, None)},
            versoes=VERSOES_DA_REPETICAO,
            parametros=tuple(c for c in PARAMETROS_DA_REPETICAO
                             if c != "canal_negativo"))
        if problemas:
            _recusar(problemas)
            return SAIDA_ERRO

    try:
        modelo = (None if args.sem_canal_negativo or not args.checkpoint
                  else _abrir_modelo(args))
        if anterior is not None:
            problemas = verificar_manifesto(
                anterior,
                {"parametros": _parametros(args, modelo is not None, None)},
                versoes=(), parametros=("canal_negativo",))
            if problemas:
                _recusar(problemas)
                return SAIDA_ERRO
            antes = (anterior.get("parametros") or {}).get("dispositivo")
            agora = modelo.dispositivo if modelo else None
            if antes != agora:
                print(f"[AVISO] a execucao original rodou em {antes!r} e esta "
                      f"roda em {agora!r}: a igualdade bit a bit nao e "
                      "prometida entre dispositivos", file=sys.stderr)
        ligantes_de_entrada = _ligantes_de_entrada(args.ligantes)
        canal = {}
        sem_fingerprint = 0
        if modelo is not None:
            ids, nomes = _alvos_afirmativos(args)
            canal, sem_fingerprint = _canal_por_composto(
                modelo, compostos,
                {k: smiles for k, (_, smiles) in ligantes_de_entrada.items()},
                _indices_afirmativos(ids, nomes, modelo.rotulos),
                args.minimo_treino, args.lote)
    except (OSError, ValueError, RuntimeError) as e:
        print(f"[ERRO] {e}", file=sys.stderr)
        return SAIDA_ERRO

    ranqueados = []
    fora = []
    # Ligante que a etapa 1 conhece e `evidencias.tsv` nao menciona: a etapa 1
    # so escreve LINHAS DE EVIDENCIA, entao quem nao tem evidencia alguma nao
    # aparece la. Sem isto ele sumia -- nem ranqueado nem listado --, que e o
    # caso mais extremo de faixa 3. Sem `--ligantes` nao ha como sabe-lo, e
    # inventar a contagem seria pior que nao te-la.
    for inchikey, (nome_de_entrada, _) in ligantes_de_entrada.items():
        if inchikey in compostos:
            continue
        aprovado, motivos = _veredito(portao, inchikey)
        fora.append({"ligante": nome_de_entrada,
                     "ligante_inchikey_parental": inchikey,
                     "motivo": "; ".join(
                         _motivo_do_portao(aprovado, motivos)
                         + ["sem evidencia na etapa 1"]),
                     "faixa": "3", "portao_aprovado": aprovado,
                     "portao_motivos": motivos,
                     "procedencia": procedencia_do_ranking()})

    for inchikey, (nome, linhas) in compostos.items():
        p_negativo, alvo_negativo = canal.get(inchikey, (None, ""))
        resultado = classificar(linhas, p_negativo=p_negativo)
        aprovado, motivos = _veredito(portao, inchikey)
        motivos_de_fora = _motivo_do_portao(aprovado, motivos)
        if resultado["faixa"] == "3":
            motivos_de_fora.append(esquema.MOTIVO_FAIXA_3)
        if motivos_de_fora:
            fora.append({"ligante": nome,
                         "ligante_inchikey_parental": inchikey,
                         "motivo": "; ".join(motivos_de_fora),
                         "faixa": resultado["faixa"],
                         "portao_aprovado": aprovado,
                         "portao_motivos": motivos,
                         "procedencia": procedencia_do_ranking()})
            continue
        ranqueados.append((chave(resultado, inchikey), nome, inchikey,
                           resultado, aprovado, motivos,
                           p_negativo, alvo_negativo))

    ranqueados.sort(key=lambda item: item[0])

    destino = Path(args.saida)
    try:
        destino.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        # Codigo 2, nao 1: 1 significa "nenhuma entrada utilizavel", e uma
        # automacao leria a causa errada.
        print(f"[ERRO] --saida nao pode ser criada: {e}", file=sys.stderr)
        return SAIDA_ERRO
    parcial = bool(linhas_descartadas or sem_fingerprint)
    nome_do_manifesto = (NOME_MANIFESTO_REPETICAO if args.manifesto
                         else NOME_MANIFESTO)
    caminhos = {nome: destino / nome
                for nome in (NOME_RANKING, NOME_NAO_RANQUEADOS,
                             nome_do_manifesto)}
    temporarios = {nome: atomico.temporario(c) for nome, c in caminhos.items()}
    try:
        with EscritorTSV(temporarios[NOME_RANKING],
                         esquema.COLUNAS_RANKING) as esc:
            for posicao, (_, nome, inchikey, resultado, aprovado, motivos,
                          p_negativo, alvo_negativo) in \
                    enumerate(ranqueados, start=1):
                esc.escrever({
                    "posicao": posicao, "ligante": nome,
                    "ligante_inchikey_parental": inchikey,
                    "faixa": resultado["faixa"],
                    "valor": _texto(resultado["valor"]),
                    "alvo": resultado["alvo"],
                    "p_afirmativo": _texto(resultado["p_afirmativo"]),
                    "pchembl": _texto(resultado["pchembl"]),
                    "p_negativo": (INDETERMINADO if p_negativo is None
                                   else _texto(p_negativo)),
                    "alvo_negativo": alvo_negativo,
                    "medida_sem_pchembl": resultado["medida_sem_pchembl"],
                    "portao_aprovado": aprovado, "portao_motivos": motivos,
                    "observacao": resultado["observacao"],
                    "procedencia": procedencia_do_ranking()})
        with EscritorTSV(temporarios[NOME_NAO_RANQUEADOS],
                         esquema.COLUNAS_NAO_RANQUEADOS) as esc:
            for linha in fora:
                esc.escrever(linha)
        # Escrita direta no temporario: a publicacao atomica das tres saidas
        # juntas e feita abaixo, e `atomico.escrever_texto` faria um segundo
        # temporario para um arquivo que ainda nem e o destino.
        Path(temporarios[nome_do_manifesto]).write_text(texto_do_manifesto(
            parametros=_parametros(
                args, modelo is not None,
                modelo.dispositivo if modelo else None) | {
                    "execucao_parcial": parcial},
            versoes=versoes,
            contagens={"compostos": len(compostos),
                       "ranqueados": len(ranqueados),
                       "nao_ranqueados": len(fora),
                       "linhas_sem_forma_parental": linhas_descartadas,
                       "sem_fingerprint": sem_fingerprint}), encoding="utf-8")
        for nome, tmp in temporarios.items():
            atomico.publicar(tmp, caminhos[nome])
    except BaseException:
        for tmp in temporarios.values():
            Path(tmp).unlink(missing_ok=True)
        raise

    print(f"\n{len(compostos)} composto(s): {len(ranqueados)} ranqueado(s), "
          f"{len(fora)} fora.")
    if parcial:
        # Codigo 3 e nao 0: automacao nao le `stderr`, e uma execucao que
        # perdeu linhas era indistinguivel de uma completa.
        print(f"[AVISO] execucao PARCIAL: {linhas_descartadas} linha(s) sem "
              f"forma parental, {sem_fingerprint} composto(s) sem fingerprint "
              "(canal negativo indeterminado neles)", file=sys.stderr)
    print(f"Saidas em {destino}")
    return SAIDA_PARCIAL if parcial else SAIDA_OK


def procedencia_do_ranking() -> str:
    """Quem produziu a linha. Fixo: esta fase nao consulta servico nenhum --
    ela ordena o que a etapa 1 mediu e previu."""
    return "fase3/ranquear.py"


if __name__ == "__main__":
    sys.exit(main())
