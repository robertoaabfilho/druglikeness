#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
checar_alertas_adversos.py

Verifica, para cada composto candidato do ranking, se a bula dele (via
openFDA) menciona o(s) termo(s) de alerta configurados (ex: "neuropathy")
nas seções de ADVERTÊNCIAS / REAÇÕES ADVERSAS / TARJA PRETA -- nunca nas
indicações. Isso resolve um problema real encontrado numa execução: drogas
famosas por CAUSAR determinado efeito (ex: Ciprofloxacino, que tem alerta
de tarja preta da FDA para neuropatia periférica) podem aparecer no topo
de um ranking que busca drogas que TRATAM esse mesmo efeito, porque o
método de pontuação (`perfil`, em ranking_multibase.py) mede semelhança
físico-química/ADMET, não mecanismo de ação -- e as duas classes de droga
podem convergir nesse espaço sem nenhuma relação farmacológica real.

O QUE ESTE SCRIPT FAZ
-----------------------
  1. Lê o ranking.csv e, para cada composto, busca a bula no openFDA
     (por nome genérico, princípio ativo ou nome comercial).
  2. Procura o(s) termo(s) de alerta nos campos de advertência/reação
     adversa/tarja preta da bula.
  3. Aplica uma PENALIDADE MULTIPLICATIVA (configurável, padrão 0.5) na
     pontuação dos compostos sinalizados, e gera um ranking ajustado.

O QUE ISSO NÃO FAZ (de propósito)
------------------------------------
Não remove silenciosamente os compostos sinalizados -- eles continuam
visíveis no CSV de saída, com o trecho da bula que gerou o alerta, para
revisão humana. Busca em texto livre não entende negação (uma bula que
diz "não há relatos de neuropatia" ainda bateria no termo "neuropatia")
-- por isso a penalidade é PARCIAL por padrão, não uma exclusão total, e
o trecho encontrado é sempre mostrado para conferência.

USO
---
    python checar_alertas_adversos.py --ranking results/ranking.csv \
        --termos "neuropathy,neuritis" --output results/ranking.alertas.csv

REQUISITOS
----------
    pip install requests
"""

import argparse
import csv
import sys
import time
from pathlib import Path
from urllib.parse import quote

import requests

try:
    from configuracao import carregar_config, get_str, get_float
except ImportError:
    print("[ERRO] Não foi possível importar 'configuracao.py'. Coloque este script na mesma pasta.", file=sys.stderr)
    sys.exit(1)

OPENFDA_BASE = "https://api.fda.gov/drug/label.json"
HEADERS_HTTP = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) checar_alertas_adversos/1.0",
    "Accept": "application/json",
}

# Campos da bula onde um "alerta" tem significado real (efeito adverso/risco).
# NUNCA inclui 'indications_and_usage' -- lá a menção ao termo seria o
# oposto do que estamos procurando (a droga tratando a condição).
CAMPOS_ALERTA = ["boxed_warning", "warnings", "warnings_and_cautions", "adverse_reactions"]

CAMPOS_NOME = ["generic_name", "substance_name", "brand_name"]


def buscar_rotulo_openfda(nome_droga: str) -> dict | None:
    """Busca a bula (label) mais provável de uma droga pelo nome, tentando
    genérico -> princípio ativo -> nome comercial, nessa ordem."""
    for campo in CAMPOS_NOME:
        try:
            url = f"{OPENFDA_BASE}"
            params = {"search": f'openfda.{campo}:"{nome_droga}"', "limit": 1}
            resposta = requests.get(url, params=params, headers=HEADERS_HTTP, timeout=30)
        except requests.exceptions.RequestException:
            continue
        if resposta.status_code == 200:
            dados = resposta.json()
            resultados = dados.get("results", [])
            if resultados:
                return resultados[0]
        # 404 = não encontrado por esse campo; tenta o próximo campo de nome.
    return None


def extrair_trecho(texto: str, termo: str, janela: int = 100) -> str:
    """Extrai um trecho curto ao redor da primeira ocorrência do termo (case-insensitive)."""
    idx = texto.lower().find(termo.lower())
    if idx == -1:
        return texto[:janela] + "..."
    inicio = max(0, idx - janela)
    fim = min(len(texto), idx + len(termo) + janela)
    prefixo = "..." if inicio > 0 else ""
    sufixo = "..." if fim < len(texto) else ""
    return prefixo + texto[inicio:fim].strip() + sufixo


def checar_alerta(nome_droga: str, termos: list[str]) -> dict | None:
    """
    Retorna {'campo': ..., 'termo': ..., 'trecho': ...} se algum termo de
    alerta aparecer nos campos de advertência/reação adversa/tarja preta
    da bula da droga. Retorna None se não encontrar nada (ou não achar a
    bula).
    """
    rotulo = buscar_rotulo_openfda(nome_droga)
    if not rotulo:
        return None

    for campo in CAMPOS_ALERTA:
        valores = rotulo.get(campo, [])
        for texto in valores:
            texto_lower = texto.lower()
            for termo in termos:
                if termo.lower() in texto_lower:
                    return {"campo": campo, "termo": termo, "trecho": extrair_trecho(texto, termo)}
    return None


def ler_ranking(caminho: str) -> tuple[str, list[dict]]:
    with open(caminho, "r", encoding="utf-8-sig", newline="") as f:
        primeira = f.readline()
        delimitador = ";" if ";" in primeira else ","
        f.seek(0)
        reader = csv.DictReader(f, delimiter=delimitador)
        linhas = list(reader)
    return delimitador, linhas


def main():
    parser = argparse.ArgumentParser(
        description="Verifica alertas de efeito adverso (bula, openFDA) para os compostos de um ranking, e ajusta a pontuação."
    )
    parser.add_argument("--ranking", default=None, help="CSV do ranking a checar. Padrão: 'ranking_output' do config.txt.")
    parser.add_argument("--termos", default=None, help="Termo(s) de alerta, separados por vírgula (ex: 'neuropathy,neuritis'). Padrão: 'termos_alerta' do config.txt, ou 'neuropathy'.")
    parser.add_argument("--penalidade", type=float, default=None, help="Multiplicador aplicado à pontuação de compostos sinalizados (0=zera, 1=sem penalidade). Padrão: 0.5, ou 'penalidade_alerta' do config.txt.")
    parser.add_argument("--output", default=None, help="CSV de saída. Padrão: <ranking>.alertas.csv")
    parser.add_argument("--config", default="config.txt", help="Caminho do config.txt (opcional).")
    args = parser.parse_args()

    config = carregar_config(args.config)
    if args.ranking is None:
        args.ranking = get_str(config, "ranking_output", "ranking.csv")
    if args.termos is None:
        args.termos = get_str(config, "termos_alerta", "neuropathy")
    termos = [t.strip() for t in args.termos.split(",") if t.strip()]
    if args.penalidade is None:
        args.penalidade = get_float(config, "penalidade_alerta", 0.5)
    if args.output is None:
        args.output = str(Path(args.ranking).with_suffix("")) + ".alertas.csv"

    if not Path(args.ranking).is_file():
        print(f"[ERRO] Arquivo não encontrado: {args.ranking}", file=sys.stderr)
        sys.exit(1)
    if not (0.0 <= args.penalidade <= 1.0):
        print(f"[ERRO] --penalidade precisa estar entre 0 e 1 (recebido: {args.penalidade}).", file=sys.stderr)
        sys.exit(1)

    delimitador, linhas = ler_ranking(args.ranking)
    print(f"[INFO] {len(linhas)} composto(s) lido(s) de '{args.ranking}'.")
    print(f"[INFO] Termo(s) de alerta: {termos}. Penalidade: {args.penalidade} (multiplicador).")

    n_sinalizados = 0
    for i, linha in enumerate(linhas, start=1):
        nome = linha.get("Nome", "")
        print(f"[INFO] ({i}/{len(linhas)}) Checando bula de '{nome}'...")
        alerta = checar_alerta(nome, termos)

        pontuacao_original = float(linha.get("Pontuacao_Final", 0) or 0)
        if alerta:
            n_sinalizados += 1
            print(f"[AVISO]   ALERTA: '{nome}' menciona '{alerta['termo']}' em {alerta['campo']}.")
            linha["Alerta_Adverso"] = "SIM"
            linha["Campo_Alerta"] = alerta["campo"]
            linha["Trecho_Alerta"] = alerta["trecho"]
            linha["Pontuacao_Ajustada"] = round(pontuacao_original * args.penalidade, 2)
        else:
            linha["Alerta_Adverso"] = "nao"
            linha["Campo_Alerta"] = ""
            linha["Trecho_Alerta"] = ""
            linha["Pontuacao_Ajustada"] = pontuacao_original

        time.sleep(0.2)  # gentileza com a API pública

    linhas.sort(key=lambda l: l["Pontuacao_Ajustada"], reverse=True)
    for i, linha in enumerate(linhas, start=1):
        linha["Ranking_Ajustado"] = i

    colunas_saida = ["Ranking_Ajustado", "Ranking", "Nome", "SMILES", "Pontuacao_Final",
                      "Pontuacao_Ajustada", "Alerta_Adverso", "Campo_Alerta", "Trecho_Alerta"]
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=colunas_saida, delimiter=delimitador)
        writer.writeheader()
        for linha in linhas:
            writer.writerow({c: linha.get(c, "") for c in colunas_saida})

    print(f"\n[INFO] {n_sinalizados} de {len(linhas)} composto(s) sinalizado(s) com alerta adverso.")
    print(f"[INFO] Ranking ajustado salvo em: {args.output}")

    print("\n[INFO] Top 10 do ranking ajustado:")
    for linha in linhas[:10]:
        marca = " ⚠" if linha["Alerta_Adverso"] == "SIM" else ""
        print(f"  {linha['Ranking_Ajustado']:>3}. {linha['Nome']:<30s} {linha['Pontuacao_Ajustada']:.2f}{marca}")


if __name__ == "__main__":
    main()
