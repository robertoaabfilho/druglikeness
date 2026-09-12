#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calibrar_deeppk.py

Automatiza a ETAPA DE COLETA para uma calibração/validação das predições
ADMET do Deep-PK: para cada droga de referência (as mesmas usadas pelo
ranking_multibase.py), busca:
  1. As predições do Deep-PK (as mesmas que entram no cálculo dos pesos).
  2. O texto farmacocinético EXPERIMENTAL disponível no PubChem (seções
     "Absorption, Distribution and Excretion", "Biological Half-Life",
     "Metabolism/Metabolites", "Pharmacokinetics", agregadas de fontes como
     DrugBank, HSDB, bulas — com a fonte original citada).

O QUE ESTE SCRIPT NÃO FAZ (de propósito)
------------------------------------------
Não tenta converter automaticamente o texto experimental em números
comparáveis direto contra cada predição do Deep-PK. Fazer isso por regex
daria uma falsa sensação de precisão -- textos de farmacocinética são
heterogêneos (unidades diferentes, populações diferentes, valores em
faixas). Em vez disso, o script gera um RELATÓRIO organizado (Markdown)
com as predições e os textos lado a lado, para que a comparação/anotação
final seja feita por um humano (ou por uma segunda passada, mais tarde,
já com uma lista clara de quais endpoints realmente têm dado experimental
disponível).

Cobertura esperada: a maioria dos medicamentos aprovados tem alguma seção
de farmacocinética no PubChem (agregada de DrugBank/HSDB/bula). Drogas
menos conhecidas, candidatos pré-clínicos ou compostos de referência mais
obscuros costumam não ter nada — o script reporta isso claramente por
droga, em vez de silenciar a ausência.

USO
---
    python calibrar_deeppk.py --referencia ranking.referencia.csv --output calibracao_deeppk.md

REQUISITOS
----------
    pip install requests
    farmaco_completo.py precisa estar na mesma pasta (reaproveita a
    consulta ao Deep-PK de lá).
"""

import argparse
import csv
import json
import sys
import time
from pathlib import Path
from urllib.parse import quote

import requests

try:
    from farmaco_completo import (
        submeter_job,
        aguardar_resultado,
        resumir_molecula,
        extrair_colunas,
        construir_colunas_deeppk,
        COLUNA_OUTRAS,
    )
except ImportError:
    print(
        "[ERRO] Não foi possível importar 'farmaco_completo.py'. "
        "Coloque farmaco_completo.py na mesma pasta que este script.",
        file=sys.stderr,
    )
    sys.exit(1)

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_VIEW_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"
HEADERS_HTTP = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) calibrar_deeppk/1.0",
    "Accept": "application/json",
}

# Seções do PubChem PUG View que carregam informação farmacocinética
# experimental (agregada de fontes como DrugBank, HSDB, bulas).
SECOES_ALVO = {
    "absorption, distribution and excretion",
    "biological half-life",
    "metabolism/metabolites",
    "pharmacokinetics",
    "protein binding",
}


# =========================================================================
# PubChem: nome -> CID -> texto experimental de farmacocinética
# =========================================================================

def buscar_cid_pubchem(nome_droga: str) -> int | None:
    try:
        url = f"{PUBCHEM_BASE}/compound/name/{quote(nome_droga)}/cids/JSON"
        resposta = requests.get(url, headers=HEADERS_HTTP, timeout=30)
        if resposta.status_code != 200:
            return None
        dados = resposta.json()
        cids = dados.get("IdentifierList", {}).get("CID", [])
        return cids[0] if cids else None
    except (requests.exceptions.RequestException, ValueError):
        return None


def _caminhar_secoes(secoes: list, encontrados: list) -> None:
    """Percorre recursivamente a árvore de seções do PUG View, coletando o
    texto de qualquer seção cujo título bata com SECOES_ALVO."""
    for secao in secoes or []:
        titulo = str(secao.get("TOCHeading", "")).strip().lower()
        if titulo in SECOES_ALVO:
            for info in secao.get("Information", []):
                valores = info.get("Value", {}).get("StringWithMarkup", [])
                fontes = info.get("Reference", [])
                nome_fonte = fontes[0].get("SourceName", "PubChem") if fontes else "PubChem"
                for v in valores:
                    texto = (v.get("String") or "").strip()
                    if texto:
                        encontrados.append({
                            "secao": secao.get("TOCHeading", ""),
                            "fonte": nome_fonte,
                            "texto": texto,
                        })
        _caminhar_secoes(secao.get("Section"), encontrados)


def buscar_texto_farmacocinetica_pubchem(cid: int) -> list[dict]:
    try:
        url = f"{PUBCHEM_VIEW_BASE}/data/compound/{cid}/JSON"
        resposta = requests.get(url, headers=HEADERS_HTTP, timeout=60)
        if resposta.status_code != 200:
            return []
        dados = resposta.json()
    except (requests.exceptions.RequestException, ValueError):
        return []

    encontrados: list[dict] = []
    _caminhar_secoes(dados.get("Record", {}).get("Section", []), encontrados)
    return encontrados


# =========================================================================
# Deep-PK: predições para as drogas de referência (reaproveita farmaco_completo.py)
# =========================================================================

def consultar_deeppk_lote(nomes_smiles: list[tuple[str, str]], pred_type: str = "admet", batch_size: int = 5) -> dict[str, dict]:
    """Consulta o Deep-PK em lotes pequenos e retorna {nome: {coluna: valor}}."""
    resultado: dict[str, dict] = {}
    colunas_parametros: list[str] | None = None
    total = len(nomes_smiles)
    indice = 0

    while indice < total:
        lote = nomes_smiles[indice: indice + batch_size]
        smiles_lote = [s for _, s in lote]
        print(f"[INFO][Deep-PK] Consultando {indice + 1}-{indice + len(lote)} de {total}...")
        try:
            job_id = submeter_job(smiles_lote, pred_type, None)
            resultados_job = aguardar_resultado(job_id)
        except Exception as erro:
            print(f"[AVISO][Deep-PK] Falha neste lote: {erro}. Lote ignorado.")
            indice += len(lote)
            continue

        if colunas_parametros is None:
            colunas_parametros = extrair_colunas(resultados_job)

        for idx_local_str, dados_molecula in resultados_job.items():
            if not idx_local_str.isdigit() or not isinstance(dados_molecula, dict):
                continue
            idx_local = int(idx_local_str)
            if idx_local >= len(lote):
                continue
            nome_droga = lote[idx_local][0]
            resultado[nome_droga] = resumir_molecula(dados_molecula, colunas_parametros)

        indice += len(lote)

    return resultado


# =========================================================================
# Leitura da referência + geração do relatório
# =========================================================================

def ler_referencia(caminho: str) -> list[tuple[str, str]]:
    delimitador = ";"
    with open(caminho, "r", encoding="utf-8-sig", newline="") as f:
        amostra = f.read(2048)
        f.seek(0)
        if "," in amostra.splitlines()[0] and ";" not in amostra.splitlines()[0]:
            delimitador = ","
        reader = csv.DictReader(f, delimiter=delimitador)
        return [(linha["Nome"], linha["SMILES"]) for linha in reader if linha.get("Nome") and linha.get("SMILES")]


def gerar_relatorio(
    nomes_smiles: list[tuple[str, str]],
    admet_deeppk: dict[str, dict],
    textos_pubchem: dict[str, list[dict]],
    caminho_saida: str,
) -> None:
    linhas_md = [
        "# Calibração Deep-PK vs. dados experimentais (PubChem)\n",
        (
            "Este relatório junta, lado a lado, as predições ADMET do Deep-PK e o texto "
            "farmacocinético experimental disponível no PubChem para cada droga de referência. "
            "A comparação numérica endpoint-a-endpoint precisa ser feita manualmente — "
            "os textos abaixo vêm de fontes heterogêneas (DrugBank, HSDB, bulas), com unidades "
            "e populações de estudo que nem sempre correspondem 1:1 aos endpoints do Deep-PK.\n"
        ),
    ]

    n_com_texto = 0
    n_sem_texto = 0

    for nome, _smiles in nomes_smiles:
        linhas_md.append(f"\n## {nome}\n")

        linhas_md.append("**Predições Deep-PK:**\n")
        predicoes = admet_deeppk.get(nome)
        if predicoes:
            for coluna, valor in predicoes.items():
                if coluna == COLUNA_OUTRAS or not str(valor).strip():
                    continue
                linhas_md.append(f"- `{coluna}`: {valor}")
        else:
            linhas_md.append("- *(sem predição Deep-PK disponível para esta droga)*")

        linhas_md.append("\n**Texto experimental (PubChem):**\n")
        textos = textos_pubchem.get(nome, [])
        if textos:
            n_com_texto += 1
            for item in textos:
                linhas_md.append(f"- *{item['secao']}* (fonte: {item['fonte']}): {item['texto']}")
        else:
            n_sem_texto += 1
            linhas_md.append("- *(nenhum texto farmacocinético encontrado no PubChem para esta droga)*")

    resumo = (
        f"\n\n---\n\n**Resumo de cobertura:** {n_com_texto} de {len(nomes_smiles)} droga(s) de referência "
        f"tinham texto experimental de farmacocinética disponível no PubChem; "
        f"{n_sem_texto} não tinham nada utilizável.\n"
    )
    linhas_md.append(resumo)

    with open(caminho_saida, "w", encoding="utf-8") as f:
        f.write("\n".join(linhas_md))
    print(f"[INFO] Relatório salvo em: {caminho_saida}")
    print(f"[INFO] Cobertura: {n_com_texto}/{len(nomes_smiles)} drogas com texto experimental encontrado.")


def main():
    parser = argparse.ArgumentParser(
        description="Coleta predições Deep-PK e texto farmacocinético experimental (PubChem) para as drogas de referência, lado a lado."
    )
    parser.add_argument("--referencia", default="ranking.referencia.csv", help="CSV de drogas de referência (saída do ranking_multibase.py). Padrão: ranking.referencia.csv")
    parser.add_argument("--output", default="calibracao_deeppk.md", help="Relatório Markdown de saída.")
    parser.add_argument("--pred-type", default="admet", choices=["absorption", "distribution", "metabolism", "excretion", "admet"], help="Tipo de predição do Deep-PK.")
    parser.add_argument("--batch-size", type=int, default=5, help="Tamanho do lote enviado ao Deep-PK.")
    args = parser.parse_args()

    if not Path(args.referencia).is_file():
        print(f"[ERRO] Arquivo não encontrado: {args.referencia}", file=sys.stderr)
        sys.exit(1)

    nomes_smiles = ler_referencia(args.referencia)
    if not nomes_smiles:
        print(f"[ERRO] Nenhuma droga de referência lida de {args.referencia}.", file=sys.stderr)
        sys.exit(1)
    print(f"[INFO] {len(nomes_smiles)} droga(s) de referência lida(s) de '{args.referencia}'.")

    # ---- Deep-PK ----
    admet_deeppk = consultar_deeppk_lote(nomes_smiles, args.pred_type, args.batch_size)

    # ---- PubChem: texto experimental ----
    textos_pubchem: dict[str, list[dict]] = {}
    for i, (nome, _smiles) in enumerate(nomes_smiles, start=1):
        print(f"[INFO][PubChem] ({i}/{len(nomes_smiles)}) Buscando texto farmacocinético de '{nome}'...")
        cid = buscar_cid_pubchem(nome)
        if not cid:
            print(f"[AVISO][PubChem]   CID não encontrado para '{nome}'.")
            textos_pubchem[nome] = []
            continue
        textos = buscar_texto_farmacocinetica_pubchem(cid)
        textos_pubchem[nome] = textos
        print(f"[INFO][PubChem]   {len(textos)} trecho(s) de texto farmacocinético encontrado(s) (CID {cid}).")
        time.sleep(0.2)

    gerar_relatorio(nomes_smiles, admet_deeppk, textos_pubchem, args.output)


if __name__ == "__main__":
    main()
