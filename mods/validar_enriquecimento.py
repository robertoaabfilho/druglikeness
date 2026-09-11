#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
validar_enriquecimento.py

Validação retrospectiva por triagem virtual (retrospective virtual
screening), com leave-one-out e decoys -- a metodologia que responde à
pergunta "essa fórmula de ranking é boa, ou uma nota de 80 é só um número
qualquer?".

IDEIA
-----
Uma pontuação isolada não é interpretável sem uma distribuição de
referência para comparar. Em vez disso:

  Para cada droga de referência conhecida (uma indicação real para a
  doença):
    1. Remove essa droga do conjunto de referência -- os PESOS são
       recalculados sem ela (leave-one-out: a droga não pode "vazar"
       informação para a própria validação).
    2. Mistura essa droga (agora tratada como "candidata desconhecida")
       com N moléculas DECOY -- drogas aprovadas para OUTRAS indicações,
       não relacionadas à doença em questão (controle negativo).
    3. Roda o ranking (mesma fórmula do ranking_multibase.py) nesse
       conjunto misto.
    4. Registra em que PERCENTIL do topo a droga conhecida ficou.

  Se a fórmula captura sinal real, droga conhecida deveria ficar
  consistentemente perto do topo (percentil alto) através de várias
  rodadas de leave-one-out -- não apenas às vezes, por acaso.

MÉTRICAS REPORTADAS
--------------------
  - Percentil médio/mediano do topo das drogas conhecidas.
  - Enrichment Factor (EF) nos top 5%/10%/20% -- quantas vezes mais
    frequente que o acaso a droga conhecida aparece nesse corte.
  - ROC-AUC -- droga conhecida vs. decoys, agregado por todas as rodadas
    (via estatística de Mann-Whitney sobre os percentis).

DECOYS
------
Por padrão, busca no ChEMBL uma amostra de drogas aprovadas (fase clínica
máxima 4) para OUTRAS indicações, excluindo qualquer nome que já esteja no
conjunto de referência. Essa é uma escolha pragmática -- outras
metodologias (ex: DUD-E, Mysinger et al. 2012) usam decoys pareados por
propriedade físico-química mas topologicamente diferentes, o que é mais
rigoroso mas exige mais infraestrutura. Isso é uma limitação conhecida
deste script, não uma omissão silenciosa.

USO
---
    python validar_enriquecimento.py --referencia ranking.referencia.csv \
        --n-decoys 50 --metodo perfil --output validacao_enriquecimento.csv

REQUISITOS
----------
    pip install rdkit requests numpy scipy
    farmaco_completo.py e ranking_multibase.py precisam estar na mesma pasta.
"""

import argparse
import csv
import sys
import time
from pathlib import Path
from urllib.parse import quote

import numpy as np
import requests
from scipy.stats import rankdata

from configuracao import carregar_config, get_str, get_int

try:
    from farmaco_completo import processar_lista_druglikeness
except ImportError:
    print("[ERRO] Não foi possível importar 'farmaco_completo.py'. Coloque-o na mesma pasta.", file=sys.stderr)
    sys.exit(1)

try:
    from ranking_multibase import calcular_pesos, gerar_ranking
except ImportError:
    print("[ERRO] Não foi possível importar 'ranking_multibase.py'. Coloque-o na mesma pasta.", file=sys.stderr)
    sys.exit(1)

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
HEADERS_HTTP = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) validar_enriquecimento/1.0",
    "Accept": "application/json",
}


# =========================================================================
# Decoys: drogas aprovadas para OUTRAS indicações (controle negativo)
# =========================================================================

def buscar_decoys_chembl(nomes_excluir: set[str], n_decoys: int, offset_inicial: int = 0) -> list[tuple[str, str]]:
    """
    Busca no ChEMBL uma amostra de moléculas aprovadas (max_phase=4), com
    estrutura conhecida, excluindo qualquer nome já presente no conjunto
    de referência. Retorna [(nome, smiles), ...].
    """
    decoys: list[tuple[str, str]] = []
    offset = offset_inicial
    limite_pagina = 100
    tentativas_sem_novidade = 0

    while len(decoys) < n_decoys and tentativas_sem_novidade < 10:
        params = {
            "max_phase": 4,
            "molecule_structures__isnull": "false",
            "limit": limite_pagina,
            "offset": offset,
            "format": "json",
        }
        try:
            resposta = requests.get(f"{CHEMBL_BASE}/molecule", params=params, headers=HEADERS_HTTP, timeout=60)
        except requests.exceptions.RequestException as erro:
            print(f"[AVISO][ChEMBL] Falha de rede ao buscar decoys: {erro}")
            break
        if resposta.status_code != 200:
            print(f"[AVISO][ChEMBL] HTTP {resposta.status_code} ao buscar decoys.")
            break

        dados = resposta.json()
        moleculas = dados.get("molecules", [])
        if not moleculas:
            break

        antes = len(decoys)
        for mol in moleculas:
            nome = mol.get("pref_name")
            estrutura = mol.get("molecule_structures") or {}
            smiles = estrutura.get("canonical_smiles")
            if not nome or not smiles:
                continue
            if nome.strip().lower() in nomes_excluir:
                continue
            decoys.append((nome, smiles))
            if len(decoys) >= n_decoys:
                break

        tentativas_sem_novidade = tentativas_sem_novidade + 1 if len(decoys) == antes else 0
        offset += limite_pagina

    return decoys[:n_decoys]


# =========================================================================
# Leave-one-out + ranking do conjunto misto (droga conhecida + decoys)
# =========================================================================

def calcular_percentual_topo(ranking: list[dict], nome_alvo: str) -> float | None:
    """Retorna o percentil de topo (100 = melhor colocado, ~100/N = pior) do nome_alvo no ranking."""
    n = len(ranking)
    for item in ranking:
        if item["Nome"] == nome_alvo:
            return round((1 - (item["Ranking"] - 1) / n) * 100, 2)
    return None


def rodar_iteracao_loo(
    nome_alvo: str,
    linhas_referencia_completa: list[dict],
    decoys_descritores: list[dict],
    colunas_farmaco_alvo: list[str],
) -> tuple[float | None, list[dict]]:
    """
    Uma rodada de leave-one-out: recalcula pesos sem 'nome_alvo', monta o
    pool (droga-alvo + decoys) e roda o ranking. Retorna (percentil, ranking_pool).
    """
    referencia_sem_alvo = [l for l in linhas_referencia_completa if l["Nome"] != nome_alvo]
    linha_alvo = next(l for l in linhas_referencia_completa if l["Nome"] == nome_alvo)

    pesos = calcular_pesos(referencia_sem_alvo, colunas_farmaco_alvo)
    if not pesos:
        return None, []

    pool = [linha_alvo] + decoys_descritores
    ranking = gerar_ranking(pool, pesos)

    percentil = calcular_percentual_topo(ranking, nome_alvo)
    return percentil, ranking


# =========================================================================
# Métricas de enriquecimento
# =========================================================================

def calcular_enrichment_factor(percentis_positivos: list[float], corte_pct: float) -> float:
    """EF no top X%: taxa de acerto observada / taxa esperada ao acaso (X/100)."""
    if not percentis_positivos:
        return float("nan")
    acertos = sum(1 for p in percentis_positivos if p >= (100 - corte_pct))
    taxa_observada = acertos / len(percentis_positivos)
    taxa_esperada = corte_pct / 100
    return round(taxa_observada / taxa_esperada, 3) if taxa_esperada > 0 else float("nan")


def calcular_auc(percentis: list[float], rotulos_positivo: list[bool]) -> float:
    """AUC via estatística de Mann-Whitney sobre os percentis pool-a-pool."""
    percentis = np.array(percentis)
    rotulos = np.array(rotulos_positivo, dtype=bool)
    n_pos, n_neg = rotulos.sum(), (~rotulos).sum()
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    ranks = rankdata(percentis)
    soma_ranks_pos = ranks[rotulos].sum()
    return round((soma_ranks_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg), 4)


# =========================================================================
# Main
# =========================================================================

def ler_referencia(caminho: str) -> list[tuple[str, str]]:
    delimitador = ";"
    with open(caminho, "r", encoding="utf-8-sig", newline="") as f:
        primeira = f.readline()
        f.seek(0)
        if "," in primeira and ";" not in primeira:
            delimitador = ","
        reader = csv.DictReader(f, delimiter=delimitador)
        return [(l["Nome"], l["SMILES"]) for l in reader if l.get("Nome") and l.get("SMILES")]


def main():
    parser = argparse.ArgumentParser(description="Validação por leave-one-out com decoys (enriquecimento) do ranking de drogas de referência.")
    parser.add_argument("--referencia", default=None, help="CSV de drogas de referência. Padrão: ranking.referencia.csv, ou derivado de 'ranking_output' do config.txt.")
    parser.add_argument("--n-decoys", type=int, default=None, help="Número de decoys (drogas de outras indicações) por rodada. Padrão: 50, ou 'n_decoys' do config.txt.")
    parser.add_argument("--output", default="validacao_enriquecimento.csv", help="CSV de saída com o detalhe por droga.")
    parser.add_argument("--config", default="config.txt", help="Caminho do config.txt (padrão: config.txt; opcional).")
    args = parser.parse_args()

    config = carregar_config(args.config)
    if args.referencia is None:
        ranking_output = get_str(config, "ranking_output", "ranking.csv")
        args.referencia = str(Path(ranking_output).with_suffix("")) + ".referencia.csv"
    if args.n_decoys is None:
        args.n_decoys = get_int(config, "n_decoys", 50)

    if not Path(args.referencia).is_file():
        print(f"[ERRO] Arquivo não encontrado: {args.referencia}", file=sys.stderr)
        sys.exit(1)

    nomes_smiles_referencia = ler_referencia(args.referencia)
    nomes_referencia = {n.strip().lower() for n, _s in nomes_smiles_referencia}
    print(f"[INFO] {len(nomes_smiles_referencia)} droga(s) de referência lida(s).")

    if len(nomes_smiles_referencia) < 3:
        print("[ERRO] São necessárias pelo menos 3 drogas de referência para leave-one-out ter sentido.", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Buscando {args.n_decoys} decoy(s) no ChEMBL (drogas aprovadas para outras indicações)...")
    decoys = buscar_decoys_chembl(nomes_referencia, args.n_decoys)
    print(f"[INFO] {len(decoys)} decoy(s) obtido(s).")
    if len(decoys) < args.n_decoys * 0.5:
        print(f"[AVISO] Menos decoys do que pedido ({len(decoys)}/{args.n_decoys}) -- resultado ainda válido, mas com menor poder estatístico.")

    print("[INFO] Calculando descritores de drogabilidade (RDKit) para referências e decoys...")
    linhas_referencia = processar_lista_druglikeness(nomes_smiles_referencia)
    linhas_decoys = processar_lista_druglikeness(decoys)
    colunas_alvo = [c for c in linhas_referencia[0].keys() if c not in ("Nome", "SMILES", "Valido")]

    print(f"[INFO] Rodando leave-one-out para {len(linhas_referencia)} droga(s) de referência...")
    resultados = []
    percentis_pool: list[float] = []
    rotulos_pool: list[bool] = []

    for i, linha in enumerate(linhas_referencia, start=1):
        nome = linha["Nome"]
        print(f"[INFO] ({i}/{len(linhas_referencia)}) Leave-one-out: '{nome}'...")
        percentil, ranking_pool = rodar_iteracao_loo(nome, linhas_referencia, linhas_decoys, colunas_alvo)
        if percentil is None:
            print(f"[AVISO]   Não foi possível calcular pesos/ranking sem '{nome}'. Ignorada.")
            continue
        resultados.append({"Nome": nome, "Percentil_Topo": percentil, "N_Decoys": len(linhas_decoys)})
        print(f"[INFO]   Percentil de topo: {percentil:.1f}%")

        # Acumula o percentil REAL desta rodada (não uma aproximação) tanto
        # da droga-alvo quanto de cada decoy, para a AUC agregada abaixo.
        percentis_pool.append(percentil)
        rotulos_pool.append(True)
        for item in ranking_pool:
            if item["Nome"] == nome:
                continue
            p_decoy = calcular_percentual_topo(ranking_pool, item["Nome"])
            if p_decoy is not None:
                percentis_pool.append(p_decoy)
                rotulos_pool.append(False)

    if not resultados:
        print("[ERRO] Nenhuma rodada de leave-one-out produziu resultado.", file=sys.stderr)
        sys.exit(1)

    with open(args.output, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["Nome", "Percentil_Topo", "N_Decoys"], delimiter=";")
        writer.writeheader()
        writer.writerows(resultados)
    print(f"[INFO] Detalhe por droga salvo em: {args.output}")

    percentis = [r["Percentil_Topo"] for r in resultados]
    n_positivos = len(percentis)

    print("\n" + "=" * 60)
    print("RESUMO DA VALIDAÇÃO POR ENRIQUECIMENTO")
    print("=" * 60)
    print(f"Percentil de topo médio das drogas conhecidas:   {np.mean(percentis):.1f}%")
    print(f"Percentil de topo mediano das drogas conhecidas: {np.median(percentis):.1f}%")

    print(f"\nRecuperação por corte de topo ({n_positivos} droga(s) retirada(s), uma por rodada):")
    print(f"  {'Corte':>8s}  {'Recuperadas':>12s}  {'% Recuperado':>13s}  {'EF (vs. acaso)':>15s}")
    tabela_cortes = []
    for corte in (1, 5, 10, 20):
        n_recuperadas = sum(1 for p in percentis if p >= (100 - corte))
        pct_recuperado = round(100 * n_recuperadas / n_positivos, 1)
        ef = calcular_enrichment_factor(percentis, corte)
        print(f"  top {corte:>3d}%   {n_recuperadas:>6d}/{n_positivos:<5d}  {pct_recuperado:>12.1f}%  {ef:>15}")
        tabela_cortes.append({
            "Corte_pct": corte, "N_Recuperadas": n_recuperadas, "N_Total": n_positivos,
            "Percentual_Recuperado": pct_recuperado, "Enrichment_Factor": ef,
        })

    auc = calcular_auc(percentis_pool, rotulos_pool)
    print(f"\nROC-AUC (droga conhecida vs. decoys): {auc}")
    print("=" * 60)
    print(
        "\nReferência de leitura: AUC=0.5 é equivalente a acaso; AUC>0.8 é considerado "
        "bom em triagem virtual retrospectiva (Mysinger et al., 2012, J. Med. Chem., "
        "metodologia de decoys DUD-E)."
    )

    caminho_cortes = str(Path(args.output).with_suffix("")) + ".cortes.csv"
    with open(caminho_cortes, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["Corte_pct", "N_Recuperadas", "N_Total", "Percentual_Recuperado", "Enrichment_Factor"], delimiter=";")
        writer.writeheader()
        writer.writerows(tabela_cortes)
    print(f"[INFO] Tabela de recuperação por corte salva em: {caminho_cortes}")


if __name__ == "__main__":
    main()
