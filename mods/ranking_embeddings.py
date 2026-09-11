#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ranking_embeddings.py

Ranking por embeddings de modelo pré-treinado em SMILES, seguindo os 5
passos propostos:

  1. Gera embeddings das drogas de referência.
  2. Representa cada candidata no mesmo espaço.
  3. Mede proximidade com cada referência (e com o grupo/centroide).
  4. Produz o ranking (por similaridade estrutural aprendida).
  5. Reordena usando ADMET (combina com a pontuação ADMET já calculada
     pelo Deep-PK, presente no farmaco.csv).

MODELO
------
Por padrão usa o ChemBERTa (seyonec/ChemBERTa-zinc-base-v1, RoBERTa
treinado em 100k SMILES do ZINC15, via HuggingFace/transformers) --
é o mais simples de carregar localmente (AutoModel/AutoTokenizer padrão,
sem infraestrutura extra). Outras opções mencionadas (MolFormer,
MegaMolBART) são compatíveis via --modelo, mas com ressalvas:
  - MolFormer (ex: "ibm/MoLFormer-XL-both-10pct") geralmente requer
    trust_remote_code=True.
  - MegaMolBART (NVIDIA) não é um checkpoint HuggingFace padrão -- é
    normalmente servido via NeMo/BioNeMo, exigindo mais infraestrutura
    para rodar localmente. Não é suportado diretamente por este script.

EMBEDDING
---------
Mean-pooling sobre os tokens (ponderado pela attention_mask), seguido de
normalização L2 -- a abordagem padrão para extrair um vetor de sentença
de um modelo tipo BERT sem um "pooler" de sentença dedicado.

REORDENAÇÃO POR ADMET (passo 5)
---------------------------------
Pontuação final = peso_embedding * similaridade_estrutural (0-100)
                 + peso_admet     * desejabilidade_admet (0-100)

A desejabilidade ADMET reaproveita a mesma lógica de pontuação por
similaridade à média das referências já usada no ranking_multibase.py
(pontuar_composto), mas RESTRITA às colunas ADMET dinâmicas (exclui os
descritores de Lipinski/Veber, que aqui já são cobertos indiretamente
pelo embedding estrutural).

USO
---
    python ranking_embeddings.py --farmaco farmaco.csv --referencia ranking.referencia.csv \
        --pesos ranking.pesos.json --output ranking_embeddings.csv

    # Sem reordenação por ADMET (só similaridade estrutural aprendida):
    python ranking_embeddings.py --farmaco farmaco.csv --referencia ranking.referencia.csv --peso-admet 0

REQUISITOS
----------
    pip install transformers torch numpy
    ranking_multibase.py precisa estar na mesma pasta.

LIMITAÇÃO CONHECIDA
--------------------
Este script não foi validado ponta-a-ponta com o modelo real (ambiente de
desenvolvimento sem acesso ao Hugging Face / sem espaço para o PyTorch).
A lógica de similaridade, ranking e reordenação por ADMET foi testada com
embeddings simulados; o carregamento e a inferência do modelo real
precisam ser verificados pelo usuário.
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

try:
    from ranking_multibase import (
        ler_farmaco_csv,
        pontuar_composto,
        COLUNAS_NAO_CARACTERISTICA,
        COLUNAS_NUMERICAS,
        COLUNAS_BOOLEANAS,
        COLUNAS_IGNORADAS,
    )
except ImportError:
    print("[ERRO] Não foi possível importar 'ranking_multibase.py'. Coloque este script na mesma pasta.", file=sys.stderr)
    sys.exit(1)

MODELO_PADRAO = "seyonec/ChemBERTa-zinc-base-v1"


# =========================================================================
# Modelo + geração de embeddings
# =========================================================================

def carregar_modelo(nome_modelo: str):
    """Carrega tokenizer + modelo do HuggingFace. Feito sob demanda (não no
    import do módulo) para que testes da lógica de similaridade/ranking não
    dependam de ter transformers/torch instalados."""
    try:
        import torch  # noqa: F401
        from transformers import AutoModel, AutoTokenizer
    except ImportError as erro:
        raise RuntimeError(
            "Este script precisa de 'transformers' e 'torch' instalados: "
            "pip install transformers torch"
        ) from erro

    print(f"[INFO] Carregando modelo '{nome_modelo}' (primeira vez pode baixar alguns GB)...")
    tokenizer = AutoTokenizer.from_pretrained(nome_modelo)
    modelo = AutoModel.from_pretrained(nome_modelo)
    modelo.eval()
    return tokenizer, modelo


def _mean_pooling(saida_modelo, attention_mask):
    import torch
    token_embeddings = saida_modelo[0]  # (batch, seq_len, hidden)
    mascara = attention_mask.unsqueeze(-1).expand(token_embeddings.size()).float()
    soma = torch.sum(token_embeddings * mascara, dim=1)
    contagem = torch.clamp(mascara.sum(dim=1), min=1e-9)
    return soma / contagem


def gerar_embeddings(lista_smiles: list[str], tokenizer, modelo, batch_size: int = 16) -> np.ndarray:
    """Gera embeddings normalizados (L2) para uma lista de SMILES. Retorna
    matriz (N x D). Requer transformers + torch (importado sob demanda)."""
    import torch

    todos_embeddings = []
    with torch.no_grad():
        for i in range(0, len(lista_smiles), batch_size):
            lote = lista_smiles[i:i + batch_size]
            entradas = tokenizer(lote, padding=True, truncation=True, return_tensors="pt")
            saida = modelo(**entradas)
            emb = _mean_pooling(saida, entradas["attention_mask"])
            emb = torch.nn.functional.normalize(emb, p=2, dim=1)
            todos_embeddings.append(emb.numpy())

    return np.vstack(todos_embeddings)


# =========================================================================
# Similaridade estrutural aprendida (passos 1-4)
# =========================================================================

def calcular_similaridades(embeddings_candidatos: np.ndarray, embeddings_referencia: np.ndarray) -> dict[str, np.ndarray]:
    """
    Retorna similaridades de cosseno (já normalizadas -> produto escalar
    direto) de cada candidata:
      - 'centroide': contra o vetor médio das referências (perfil do grupo).
      - 'max_individual': contra a referência individual mais próxima.
    """
    centroide = embeddings_referencia.mean(axis=0)
    centroide = centroide / (np.linalg.norm(centroide) + 1e-12)

    sim_centroide = embeddings_candidatos @ centroide  # (N,)
    sim_por_referencia = embeddings_candidatos @ embeddings_referencia.T  # (N, R)
    sim_max_individual = sim_por_referencia.max(axis=1)
    indice_mais_proxima = sim_por_referencia.argmax(axis=1)

    return {
        "centroide": sim_centroide,
        "max_individual": sim_max_individual,
        "indice_mais_proxima": indice_mais_proxima,
    }


# =========================================================================
# Reordenação por ADMET (passo 5)
# =========================================================================

def filtrar_pesos_admet(pesos: dict[str, dict]) -> dict[str, dict]:
    """Mantém só as colunas ADMET dinâmicas (exclui os descritores fixos de
    Lipinski/Veber, já cobertos pelo embedding estrutural)."""
    descritores_fixos = COLUNAS_NUMERICAS | COLUNAS_BOOLEANAS | COLUNAS_NAO_CARACTERISTICA | COLUNAS_IGNORADAS
    return {c: info for c, info in pesos.items() if c not in descritores_fixos}


def calcular_desejabilidade_admet(linhas_farmaco: list[dict], pesos_admet: dict[str, dict]) -> np.ndarray:
    if not pesos_admet:
        return np.zeros(len(linhas_farmaco))
    valores = []
    for linha in linhas_farmaco:
        pontuacao, _detalhes = pontuar_composto(linha, pesos_admet)
        valores.append(pontuacao)
    return np.array(valores)


def _normalizar_0_100(vetor: np.ndarray) -> np.ndarray:
    minimo, maximo = vetor.min(), vetor.max()
    if maximo - minimo < 1e-12:
        return np.full_like(vetor, 50.0)
    return (vetor - minimo) / (maximo - minimo) * 100


# =========================================================================
# Leitura auxiliar + main
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
    parser = argparse.ArgumentParser(description="Ranking por embeddings (ChemBERTa/MolFormer) + reordenação por ADMET.")
    parser.add_argument("--farmaco", default="farmaco.csv", help="CSV de compostos candidatos.")
    parser.add_argument("--referencia", default="ranking.referencia.csv", help="CSV de drogas de referência.")
    parser.add_argument("--pesos", default="ranking.pesos.json", help="JSON de pesos (para a parte ADMET do reordenamento).")
    parser.add_argument("--output", default="ranking_embeddings.csv", help="CSV de saída.")
    parser.add_argument("--modelo", default=MODELO_PADRAO, help=f"Checkpoint HuggingFace a usar. Padrão: {MODELO_PADRAO}")
    parser.add_argument("--similaridade", default="centroide", choices=["centroide", "max_individual"], help="Medir proximidade contra o centroide do grupo ou contra a referência individual mais próxima.")
    parser.add_argument("--peso-embedding", type=float, default=0.5, help="Peso da similaridade estrutural (embedding) na pontuação final. Padrão: 0.5")
    parser.add_argument("--peso-admet", type=float, default=0.5, help="Peso da desejabilidade ADMET na pontuação final. Padrão: 0.5")
    parser.add_argument("--batch-size", type=int, default=16, help="Tamanho do lote para inferência do modelo.")
    args = parser.parse_args()

    for caminho, nome in [(args.farmaco, "farmaco"), (args.referencia, "referência")]:
        if not Path(caminho).is_file():
            print(f"[ERRO] Arquivo de {nome} não encontrado: {caminho}", file=sys.stderr)
            sys.exit(1)

    _colunas, linhas_farmaco = ler_farmaco_csv(args.farmaco)
    nomes_smiles_referencia = ler_referencia(args.referencia)
    print(f"[INFO] {len(linhas_farmaco)} composto(s) candidato(s), {len(nomes_smiles_referencia)} droga(s) de referência.")

    # ---- Passos 1-2: embeddings ----
    tokenizer, modelo = carregar_modelo(args.modelo)
    smiles_referencia = [s for _n, s in nomes_smiles_referencia]
    smiles_candidatos = [l.get("SMILES", "") for l in linhas_farmaco]

    print("[INFO] Gerando embeddings das drogas de referência...")
    embeddings_referencia = gerar_embeddings(smiles_referencia, tokenizer, modelo, args.batch_size)
    print("[INFO] Gerando embeddings dos compostos candidatos...")
    embeddings_candidatos = gerar_embeddings(smiles_candidatos, tokenizer, modelo, args.batch_size)

    # ---- Passo 3: proximidade ----
    similaridades = calcular_similaridades(embeddings_candidatos, embeddings_referencia)
    sim_escolhida = similaridades[args.similaridade]
    sim_normalizada = _normalizar_0_100((sim_escolhida + 1) / 2)  # cosseno [-1,1] -> [0,1] -> [0,100]

    # ---- Passo 4: ranking por similaridade estrutural ----
    # ---- Passo 5: reordenação por ADMET ----
    pesos_admet = {}
    if args.peso_admet > 0:
        if not Path(args.pesos).is_file():
            print(f"[AVISO] --peso-admet > 0 mas arquivo de pesos não encontrado ({args.pesos}) -- ADMET será ignorado.")
        else:
            with open(args.pesos, "r", encoding="utf-8") as f:
                pesos_completos = json.load(f)
            pesos_admet = filtrar_pesos_admet(pesos_completos)
            print(f"[INFO] {len(pesos_admet)} coluna(s) ADMET usada(s) na reordenação.")

    desejabilidade_admet = calcular_desejabilidade_admet(linhas_farmaco, pesos_admet)
    desejabilidade_admet_norm = _normalizar_0_100(desejabilidade_admet) if pesos_admet else np.zeros(len(linhas_farmaco))

    soma_pesos = args.peso_embedding + args.peso_admet
    if soma_pesos <= 0:
        print("[ERRO] peso-embedding + peso-admet precisa ser > 0.", file=sys.stderr)
        sys.exit(1)
    pontuacao_final = (args.peso_embedding * sim_normalizada + args.peso_admet * desejabilidade_admet_norm) / soma_pesos

    resultado = []
    for i, linha in enumerate(linhas_farmaco):
        idx_mais_proxima = int(similaridades["indice_mais_proxima"][i])
        resultado.append({
            "Nome": linha.get("Nome", ""),
            "SMILES": linha.get("SMILES", ""),
            "Similaridade_Embedding": round(float(sim_normalizada[i]), 2),
            "Desejabilidade_ADMET": round(float(desejabilidade_admet_norm[i]), 2),
            "Pontuacao_Final": round(float(pontuacao_final[i]), 2),
            "Referencia_Mais_Proxima": nomes_smiles_referencia[idx_mais_proxima][0],
        })

    resultado.sort(key=lambda r: r["Pontuacao_Final"], reverse=True)
    for i, item in enumerate(resultado, start=1):
        item["Ranking"] = i

    colunas_saida = ["Ranking", "Nome", "SMILES", "Pontuacao_Final", "Similaridade_Embedding", "Desejabilidade_ADMET", "Referencia_Mais_Proxima"]
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=colunas_saida, delimiter=";")
        writer.writeheader()
        for item in resultado:
            writer.writerow({c: item[c] for c in colunas_saida})
    print(f"[INFO] Ranking salvo em: {args.output}")

    print("\n[INFO] Top 10:")
    for item in resultado[:10]:
        print(f"  {item['Ranking']:>3}. {item['Nome']:<30s} {item['Pontuacao_Final']:.2f}  (mais próxima: {item['Referencia_Mais_Proxima']})")


if __name__ == "__main__":
    main()
