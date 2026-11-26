#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ranquear.py

Fase 3: ordena compostos pela evidencia de alvo da etapa 1, em faixas
lexicograficas, e publica `ranking.tsv`, `nao_ranqueados.tsv` e o
`manifesto.json` dentro de `--saida`.

Desenho: `docs/superpowers/specs/2026-11-26-fase3-ranqueamento-design.md`.

O que esta fase existe para corrigir
------------------------------------
O ranqueamento antigo (`ranking_multibase.py`) ordena por semelhanca
fisico-quimica com as drogas de referencia do sintoma, e e cego a alvo. Medido
no controle de 12 farmacos: os dois QED mais altos sao ciprofloxacino (0,893) e
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
from mods.evidencias.manifesto import ler_manifesto, texto_do_manifesto
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
    erro = _erro_de_checkpoint(args.checkpoint, args.manifesto_etapa1)
    if erro:
        return erro
    for rotulo, caminho in (("--evidencias", args.evidencias),
                            ("--portao", args.portao),
                            ("--ligantes", args.ligantes),
                            ("--rotulos", args.rotulos),
                            ("--manifesto-etapa1", args.manifesto_etapa1),
                            ("--alvos", args.alvos)):
        if caminho and not Path(caminho).is_file():
            return f"{rotulo}: arquivo nao encontrado: {caminho}"
    destino = Path(args.saida)
    if destino.exists() and not destino.is_dir():
        return (f"--saida tem de ser um diretorio; {args.saida} e um arquivo "
                "que existe")
    for rotulo, caminho in (("--evidencias", args.evidencias),
                            ("--portao", args.portao)):
        if caminho and Path(caminho).resolve().parent == destino.resolve():
            return (f"--saida e o diretorio de {rotulo}; a publicacao "
                    "sobrescreveria a entrada")
    return None


def _por_composto(caminho, limite=None):
    """{inchikey: (nome, [linhas])}, na ordem em que os compostos aparecem."""
    compostos: dict[str, tuple[str, list[dict]]] = {}
    for linha in ler_tsv(caminho):
        chave_composto = (linha.get("ligante_inchikey_parental") or "").strip()
        nome = (linha.get("ligante") or "").strip()
        if not chave_composto:
            # Sem a forma parental nao ha chave de juncao: a etapa 1 ja
            # relatou o caso, e aqui a linha nao tem a quem pertencer.
            continue
        if chave_composto not in compostos:
            if limite is not None and len(compostos) >= limite:
                continue
            compostos[chave_composto] = (nome, [])
        compostos[chave_composto][1].append(linha)
    return compostos


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


def _abrir_modelo(args) -> Modelo | None:
    """O PP4CKD, ou None quando esta execucao nao tem canal negativo.

    O import de `torch` fica aqui dentro: uma execucao com
    `--sem-canal-negativo` nao paga por ele, e os testes substituem esta
    funcao em vez de montar um checkpoint de 27 MB."""
    if args.sem_canal_negativo or not args.checkpoint:
        return None
    import torch

    from mods.prever_alvos import carregar_modelo, carregar_rotulos

    rotulos = carregar_rotulos(args.rotulos)
    dispositivo = "cuda" if args.gpu else "cpu"
    rede, cfg, dispositivo, nome = carregar_modelo(
        None, "ECFP4", dispositivo, checkpoint=args.checkpoint)

    def prever(fingerprints):
        import numpy as np
        lote = torch.tensor(np.array(fingerprints), dtype=torch.float32)
        with torch.no_grad():
            return torch.sigmoid(rede(lote.to(dispositivo))).cpu().numpy()

    return Modelo(rotulos=rotulos, cfg=cfg, dispositivo=dispositivo,
                  prever=prever, checkpoint=nome)


def _smiles_por_composto(caminho) -> dict[str, str]:
    """{inchikey parental: smiles} de `ligantes.tsv`."""
    if not caminho:
        return {}
    por_composto = {}
    for linha in ler_tsv(caminho):
        chave_composto = (linha.get("inchikey_parental")
                          or linha.get("ligante_inchikey_parental") or "").strip()
        smiles = (linha.get("smiles") or "").strip()
        if chave_composto and chave_composto not in por_composto:
            por_composto[chave_composto] = smiles
    return por_composto


def _indices_afirmativos(compostos, rotulos) -> list[int]:
    """Indices, no modelo, dos alvos AFIRMATIVOS desta execucao.

    Afirmativo e o alvo que o usuario pediu ou que o enriquecimento trouxe --
    a procedencia decide, e e por isso que esta lista sai das evidencias e nao
    de um filtro de organismo."""
    afirmativas = set(ORIGENS_DO_USUARIO) | set(ORIGENS_DE_ENRIQUECIMENTO)
    ids = set()
    nomes = set()
    for _, linhas in compostos.values():
        for linha in linhas:
            if (linha.get("origem_alvo") or "").strip() not in afirmativas:
                continue
            id_alvo = (linha.get("alvo_chembl_id") or "").strip().upper()
            if id_alvo:
                ids.add(id_alvo)
            nome = (linha.get("alvo") or "").strip().casefold()
            if nome:
                nomes.add(nome)
    indices = []
    for i, r in enumerate(rotulos):
        if (r.get("chembl_id") or "").strip().upper() in ids:
            indices.append(i)
        elif (r.get("nome") or "").strip().casefold() in nomes:
            indices.append(i)
    return indices


def _canal_por_composto(modelo, compostos, smiles_por_composto,
                        minimo_treino: int, lote: int):
    """{inchikey: (p_negativo ou None, nome do alvo)}.

    So o agregado por composto e guardado: a matriz inteira seriam
    108.935 x 7.676 (~3,3 GB em float32), e ela nao e necessaria."""
    from mods.prever_alvos import fingerprint_de

    mascara = negativo.mascara(
        modelo.rotulos, _indices_afirmativos(compostos, modelo.rotulos),
        minimo_treino=minimo_treino)
    resultado: dict[str, tuple[float | None, str]] = {}
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
            # zero afirmaria "nenhuma atividade indesejada".
            resultado[inchikey] = (None, "")
            continue
        pendentes.append((inchikey, fp))
        if len(pendentes) >= lote:
            despachar()
    despachar()
    return resultado


def _veredito(portao, inchikey) -> tuple[str, str]:
    if portao is None:
        return esquema.PORTAO_NAO_CONSULTADO, ""
    if inchikey not in portao:
        return esquema.AUSENTE_DO_PORTAO, ""
    return portao[inchikey]


def _texto(valor) -> str:
    """Numero numa celula, sem notacao cientifica e sem zeros a esquerda."""
    if valor is None or valor == "":
        return ""
    return f"{float(valor):.6f}".rstrip("0").rstrip(".")


def main(argv=None) -> int:
    args = analisar(argv)
    erro = _validar(args)
    if erro:
        print(f"[ERRO] {erro}", file=sys.stderr)
        return SAIDA_ERRO

    try:
        compostos = _por_composto(args.evidencias, args.limite)
        portao = _portao(args.portao)
    except (OSError, ValueError) as e:
        print(f"[ERRO] {e}", file=sys.stderr)
        return SAIDA_ERRO

    if not compostos:
        print("[ERRO] nenhum composto utilizavel em " + args.evidencias,
              file=sys.stderr)
        return SAIDA_SEM_ENTRADA

    try:
        modelo = _abrir_modelo(args)
        canal = ({} if modelo is None else _canal_por_composto(
            modelo, compostos, _smiles_por_composto(args.ligantes),
            args.minimo_treino, args.lote))
    except (OSError, ValueError, RuntimeError) as e:
        print(f"[ERRO] {e}", file=sys.stderr)
        return SAIDA_ERRO

    ranqueados = []
    fora = []
    for inchikey, (nome, linhas) in compostos.items():
        p_negativo, alvo_negativo = canal.get(inchikey, (None, ""))
        resultado = classificar(linhas, p_negativo=p_negativo)
        aprovado, motivos = _veredito(portao, inchikey)
        motivos_de_fora = []
        if aprovado == "nao":
            motivos_de_fora.append(f"reprovado no portao ({motivos})"
                                   if motivos else "reprovado no portao")
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
    destino.mkdir(parents=True, exist_ok=True)
    caminhos = {nome: destino / nome
                for nome in (NOME_RANKING, NOME_NAO_RANQUEADOS, NOME_MANIFESTO)}
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
        Path(temporarios[NOME_MANIFESTO]).write_text(texto_do_manifesto(
            parametros={"canal_negativo": modelo is not None,
                        "dispositivo": modelo.dispositivo if modelo else None,
                        "minimo_treino": args.minimo_treino,
                        "lote": args.lote,
                        "limite": args.limite,
                        "portao": bool(args.portao)},
            versoes=procedencia.versoes_do_ambiente() | {
                "pipeline_commit": procedencia.commit_do_pipeline(),
                "evidencias_md5": procedencia.md5(args.evidencias),
                "portao_md5": procedencia.md5_ou_none(args.portao),
                "alvos_md5": procedencia.md5_ou_none(args.alvos)},
            contagens={"compostos": len(compostos),
                       "ranqueados": len(ranqueados),
                       "nao_ranqueados": len(fora)}), encoding="utf-8")
        for nome, tmp in temporarios.items():
            atomico.publicar(tmp, caminhos[nome])
    except BaseException:
        for tmp in temporarios.values():
            Path(tmp).unlink(missing_ok=True)
        raise

    print(f"\n{len(compostos)} composto(s): {len(ranqueados)} ranqueado(s), "
          f"{len(fora)} fora.")
    print(f"Saidas em {destino}")
    return SAIDA_OK


def procedencia_do_ranking() -> str:
    """Quem produziu a linha. Fixo: esta fase nao consulta servico nenhum --
    ela ordena o que a etapa 1 mediu e previu."""
    return "fase3/ranquear.py"


if __name__ == "__main__":
    sys.exit(main())
