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

# Permite `python mods/ranquear.py`: a raiz do repositorio entra no caminho.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from mods.evidencias import atomico, procedencia
from mods.evidencias.manifesto import texto_do_manifesto
from mods.portao.esquema import EscritorTSV, ler_tsv
from mods.ranque import esquema
from mods.ranque.faixa import chave, classificar

SAIDA_OK = 0
SAIDA_SEM_ENTRADA = 1
SAIDA_ERRO = 2
SAIDA_PARCIAL = 3

NOME_RANKING = "ranking.tsv"
NOME_NAO_RANQUEADOS = "nao_ranqueados.tsv"
NOME_MANIFESTO = "manifesto.json"


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
    ap.add_argument("--sem-canal-negativo", action="store_true",
                    help="desliga o canal negativo; a faixa prevista ordena "
                         "pelo p sozinho. Medido, o canal leva a AUC de 0,750 "
                         "para 0,806 no controle de 12 -- mas n = 12 e os seis "
                         "causadores sao todos antibioticos, entao a regra e "
                         "desligavel de proposito")
    return ap.parse_args(argv)


def _validar(args) -> str | None:
    """Mensagem de erro, ou None. Roda ANTES de abrir qualquer saida."""
    for rotulo, caminho in (("--evidencias", args.evidencias),
                            ("--portao", args.portao),
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

    ranqueados = []
    fora = []
    for inchikey, (nome, linhas) in compostos.items():
        resultado = classificar(linhas, p_negativo=None)
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
                           resultado, aprovado, motivos))

    ranqueados.sort(key=lambda item: item[0])

    destino = Path(args.saida)
    destino.mkdir(parents=True, exist_ok=True)
    caminhos = {nome: destino / nome
                for nome in (NOME_RANKING, NOME_NAO_RANQUEADOS, NOME_MANIFESTO)}
    temporarios = {nome: atomico.temporario(c) for nome, c in caminhos.items()}
    try:
        with EscritorTSV(temporarios[NOME_RANKING],
                         esquema.COLUNAS_RANKING) as esc:
            for posicao, (_, nome, inchikey, resultado, aprovado, motivos) in \
                    enumerate(ranqueados, start=1):
                esc.escrever({
                    "posicao": posicao, "ligante": nome,
                    "ligante_inchikey_parental": inchikey,
                    "faixa": resultado["faixa"],
                    "valor": _texto(resultado["valor"]),
                    "alvo": resultado["alvo"],
                    "p_afirmativo": _texto(resultado["p_afirmativo"]),
                    "pchembl": _texto(resultado["pchembl"]),
                    "p_negativo": "indeterminado",
                    "alvo_negativo": "",
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
            parametros={"canal_negativo": not args.sem_canal_negativo,
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
