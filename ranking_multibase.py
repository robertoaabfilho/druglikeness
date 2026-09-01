#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ranking_multibase.py

Versão sem DrugBank: consulta Open Targets, ChEMBL, DrugCentral e openFDA
SIMULTANEAMENTE para achar drogas indicadas para uma doença, mescla os
resultados (deduplicados por nome), deriva pesos por característica a
partir delas, e ranqueia os compostos de um 'farmaco.csv'.

FONTES (todas gratuitas, sem chave de API)
-------------------------------------------
1. Open Targets Platform (GraphQL) — https://api.platform.opentargets.org
   Busca a doença por nome -> disease.drugAndClinicalCandidates (drogas
   conhecidas + fase clínica máxima). Retorna o ID ChEMBL de cada droga.
2. ChEMBL (REST) — https://www.ebi.ac.uk/chembl/api/data
   Endpoint drug_indication, filtrando por termo da doença (EFO/MeSH).
   Também usado para obter o SMILES de qualquer droga via seu ID ChEMBL.
3. DrugCentral — via conexão direta ao Postgres público de só-leitura
   (unmtid-dbs.net:5433), já que não há uma API REST oficial estável.
   Requer 'psycopg2'. Se não estiver instalado, esta fonte é ignorada
   automaticamente (as outras continuam funcionando).
4. openFDA (REST) — https://api.fda.gov/drug/label.json
   Busca no texto de bula ('indications_and_usage') por menção à doença.
   Não retorna SMILES nem ID ChEMBL — usa o fallback via PubChem.

Uma droga encontrada em mais de uma fonte é listada com suas fontes de
origem (Open Targets / ChEMBL / DrugCentral / openFDA) — útil como um
indicador informal de confiabilidade da referência.

USO
---
    python ranking_multibase.py --doenca "Alzheimer's disease" \
        --farmaco farmaco.csv --output ranking.csv

    # Desativar alguma fonte:
    python ranking_multibase.py --doenca Alzheimer --no-drugcentral --no-openfda

    # Sem consultar ADMET das drogas de referência (mais rápido):
    python ranking_multibase.py --doenca Alzheimer --no-admet

REQUISITOS
----------
    pip install rdkit requests
    pip install psycopg2-binary   # opcional, só para a fonte DrugCentral

DEPENDÊNCIA DE ARQUIVO
-----------------------
Este script IMPORTA `farmaco_completo.py` (mesma pasta) para calcular os
descritores de drogabilidade das drogas de referência — a mesma função
(`processar_lista_druglikeness`) usada para gerar o próprio farmaco.csv,
evitando duplicar essa lógica em dois lugares.
"""

import argparse
import csv
import json
import math
import os
import re
import statistics
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path
from urllib.parse import quote

import requests
from rdkit import Chem
from rdkit.Chem import Descriptors

try:
    import psycopg2
    PSYCOPG2_DISPONIVEL = True
except ImportError:
    PSYCOPG2_DISPONIVEL = False

# A Etapa 3 (descritores de drogabilidade locais) reaproveita diretamente o
# farmaco_completo.py — mesma lógica, mesmos nomes de coluna, sem duplicar
# código. Este arquivo precisa estar na mesma pasta (ou no PYTHONPATH).
try:
    from farmaco_completo import processar_lista_druglikeness
except ImportError:
    print(
        "[ERRO] Não foi possível importar 'farmaco_completo.py'. "
        "Coloque farmaco_completo.py na mesma pasta que este script "
        "(ou no PYTHONPATH) e tente novamente.",
        file=sys.stderr,
    )
    sys.exit(1)

# =========================================================================
# CONFIGURAÇÃO GERAL
# =========================================================================

OPENTARGETS_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"
CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
OPENFDA_BASE = "https://api.fda.gov/drug/label.json"

DRUGCENTRAL_DBHOST = "unmtid-dbs.net"
DRUGCENTRAL_DBPORT = 5433
DRUGCENTRAL_DBNAME = "drugcentral"
DRUGCENTRAL_DBUSER = "drugman"
DRUGCENTRAL_DBPASS = "dosage"

DEEPPK_API_URL = "https://biosig.lab.uq.edu.au/deeppk/api/predict"
DEEPPK_POLL_INTERVAL = 15
DEEPPK_MAX_WAIT = 30 * 60
DEEPPK_BATCH_SIZE = 5

CSV_DELIMITADOR_PADRAO = ";"

# Alguns servidores (EBI/Open Targets) tratam requisições sem um User-Agent
# "normal" como tráfego suspeito. Isso não custa nada e evita esse tipo de
# bloqueio silencioso.
HEADERS_HTTP = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) ranking_multibase/1.0",
    "Accept": "application/json",
}

COLUNAS_NAO_CARACTERISTICA = {"Nome", "SMILES", "Valido"}
COLUNAS_NUMERICAS = {
    "PesoMolecular", "DoadoresHB", "AceptoresHB", "TPSA",
    "LigacoesRotacionaveis", "RefratividadeMolar", "AtomosPesados",
}
COLUNAS_BOOLEANAS = {"Passa_Ro5_parcial", "Passa_Veber"}
COLUNAS_IGNORADAS = {"Violacoes_Ro5_parcial", "Outras_Propriedades"}


# =========================================================================
# PARTE 1a — Open Targets: doença -> drogas conhecidas (ChEMBL IDs)
# =========================================================================

def _graphql(query: str, variables: dict) -> dict:
    resposta = requests.post(
        OPENTARGETS_GRAPHQL,
        json={"query": query, "variables": variables},
        headers=HEADERS_HTTP,
        timeout=60,
    )
    # Servidores GraphQL (este usa Sangria) costumam responder HTTP 400 quando
    # a PRÓPRIA QUERY tem um erro de sintaxe/validação, e o corpo da resposta
    # traz a mensagem exata do que está errado. Sem mostrar esse corpo, um 400
    # vira uma caixa-preta — por isso capturamos e exibimos o texto aqui.
    if resposta.status_code >= 400:
        raise RuntimeError(
            f"HTTP {resposta.status_code} do Open Targets. Corpo da resposta: {resposta.text[:1000]}"
        )
    payload = resposta.json()
    if "errors" in payload:
        raise RuntimeError(f"Erro da API do Open Targets: {payload['errors']}")
    return payload["data"]


def _open_targets_buscar_efo_id(termo: str) -> str | None:
    query = """
    query BuscarDoenca($q: String!) {
      search(queryString: $q, entityNames: ["disease"], page: {index: 0, size: 5}) {
        hits { id name entity }
      }
    }
    """
    dados = _graphql(query, {"q": termo})
    hits = dados.get("search", {}).get("hits", [])
    if not hits:
        return None
    return hits[0]["id"]  # melhor correspondência


def _ordinal_fase_clinica(valor) -> float:
    """
    Converte uma fase clinica (string do Open Targets tipo 'Phase III'/
    'Approved', ou numero do ChEMBL 0-4) em um valor comparavel -- quanto
    maior, mais avancada/relevante a indicacao. Usado para ordenar e
    escolher as top N mais relevantes de cada fonte.
    """
    if valor is None:
        return -1.0
    if isinstance(valor, (int, float)):
        return float(valor)
    texto = str(valor).strip().lower()
    if not texto or texto in ("none", "-", "n/a", "unknown", "null"):
        return -1.0
    if "approved" in texto or "aprovad" in texto:
        return 5.0
    m = re.search(r"phase\s+(iv|iii|ii|i|\d+(\.\d+)?)", texto)
    if m:
        valor_fase = m.group(1)
        romanos = {"iv": 4.0, "iii": 3.0, "ii": 2.0, "i": 1.0}
        if valor_fase in romanos:
            return romanos[valor_fase]
        return float(valor_fase)
    if "preclinical" in texto or "pre-clinic" in texto:
        return 0.0
    m = re.search(r"(\d+(\.\d+)?)", texto)
    if m:
        return float(m.group(1))
    return -1.0


def _open_targets_drogas_para_efo(efo_id: str) -> dict[str, dict]:
    # 'drugAndClinicalCandidates' NAO aceita argumento de paginacao (confirmado
    # pelo proprio erro da API: "Unknown argument 'page' on field
    # 'drugAndClinicalCandidates'") -- retorna a lista inteira, sem corte
    # nosso: pegamos tudo e so depois decidimos o que e mais relevante.
    query = """
    query DrogasDaDoenca($efoId: String!) {
      disease(efoId: $efoId) {
        name
        drugAndClinicalCandidates {
          count
          rows {
            maxClinicalStage
            drug { id name }
          }
        }
      }
    }
    """
    dados = _graphql(query, {"efoId": efo_id})
    doenca = dados.get("disease")
    if not doenca:
        return {}

    encontrados = {}
    for linha in doenca.get("drugAndClinicalCandidates", {}).get("rows", []):
        drug = linha.get("drug")
        if not drug or not drug.get("name"):
            continue
        chave = drug["name"].strip().lower()
        encontrados[chave] = {
            "nome": drug["name"],
            "chembl_id": drug.get("id"),
            "fase_clinica": linha.get("maxClinicalStage"),
        }
    return encontrados


def buscar_open_targets(keywords: list[str], top_n: int) -> tuple[dict[str, dict], list[str]]:
    """
    Busca TODAS as drogas indicadas para a doenca no Open Targets (sem corte),
    depois ordena por fase clinica maxima (Aprovado > Fase IV > ... >
    Pre-clinico) e mantem so as 'top_n' mais relevantes.
    Retorna (encontrados, efo_ids_resolvidos).
    """
    todas: dict[str, dict] = {}
    efo_ids_resolvidos: list[str] = []
    for termo in keywords:
        print(f"[INFO][OpenTargets] Buscando doenca por: '{termo}'...")
        try:
            efo_id = _open_targets_buscar_efo_id(termo)
            if not efo_id:
                print(f"[AVISO][OpenTargets] Nenhuma doenca encontrada para '{termo}'.")
                continue
            efo_ids_resolvidos.append(efo_id)
            drogas = _open_targets_drogas_para_efo(efo_id)
            print(f"[INFO][OpenTargets]   {len(drogas)} droga(s) encontrada(s) no total para '{termo}' ({efo_id}).")
            todas.update(drogas)
        except Exception as erro:
            print(f"[AVISO][OpenTargets] Falha ao consultar '{termo}': {erro}")

    ordenadas = sorted(todas.items(), key=lambda kv: _ordinal_fase_clinica(kv[1].get("fase_clinica")), reverse=True)
    top = dict(ordenadas[:top_n])
    print(f"[INFO][OpenTargets] Mantendo as {len(top)} mais relevantes (maior fase clinica) de {len(todas)} encontradas.")
    return top, efo_ids_resolvidos


# =========================================================================
# PARTE 1b — ChEMBL: drug_indication filtrando por termo da doença
# =========================================================================

def _chembl_buscar_por_campo(campo_filtro: str, termo: str, encontrados: dict, limite_seguranca: int = 3000) -> None:
    """
    Consulta /drug_indication filtrando por um campo (efo_id, efo_term ou
    mesh_heading) com '__icontains', paginando ATE O FIM (sem cortar cedo --
    a escolha de quais manter e feita depois, por relevancia clinica).
    'limite_seguranca' e so uma trava para nunca rodar indefinidamente.
    Imprime o motivo em caso de falha, em vez de falhar silenciosamente.
    """
    offset = 0
    limite_pagina = 100
    while offset < limite_seguranca:
        params = {campo_filtro: termo, "limit": limite_pagina, "offset": offset, "format": "json"}
        resposta = requests.get(f"{CHEMBL_BASE}/drug_indication", params=params, headers=HEADERS_HTTP, timeout=60)
        if resposta.status_code != 200:
            print(
                f"[AVISO][ChEMBL] HTTP {resposta.status_code} ao buscar por {campo_filtro}='{termo}'. "
                f"Corpo: {resposta.text[:300]}"
            )
            return
        dados = resposta.json()
        itens = dados.get("drug_indications", [])
        if not itens:
            if offset == 0:
                print(f"[INFO][ChEMBL]   0 resultados para {campo_filtro}='{termo}'. Resposta bruta (trecho): {resposta.text[:200]}")
            return
        for item in itens:
            molecule_chembl_id = item.get("molecule_chembl_id")
            if not molecule_chembl_id:
                continue
            registro = encontrados.setdefault(
                molecule_chembl_id, {"chembl_id": molecule_chembl_id, "nome": None, "fase_clinica": None}
            )
            fase = item.get("max_phase_for_ind")
            if fase is not None and (registro["fase_clinica"] is None or fase > registro["fase_clinica"]):
                registro["fase_clinica"] = fase
        offset += limite_pagina
        if len(itens) < limite_pagina:
            return


def buscar_chembl_indicacao(keywords: list[str], top_n: int, efo_ids: list[str] | None = None) -> dict[str, dict]:
    """
    Busca TODAS as indicacoes no ChEMBL para a doenca (sem corte durante a
    coleta), combinando 3 estrategias (mais robusto que depender de uma so):
      1. efo_id (a partir dos IDs EFO ja resolvidos pelo Open Targets) --
         correspondencia por IDENTIFICADOR, nao por texto, entao nao depende
         de acertar a grafia exata do termo em ingles.
      2. efo_term__icontains -- texto livre no termo EFO.
      3. mesh_heading__icontains -- o MeSH e descrito pelo proprio ChEMBL
         como o identificador PRIMARIO de indicacao (o EFO nem sempre esta
         mapeado), entao usar so efo_term pode perder resultados.
    Depois ordena por 'max_phase_for_ind' (fase clinica maxima da indicacao;
    4 = aprovado) e mantem so as 'top_n' mais relevantes.
    Retorna {chembl_id: {nome, chembl_id, fase_clinica}}.
    """
    todas: dict[str, dict] = {}

    for efo_id in efo_ids or []:
        # IDs EFO costumam vir como "EFO_0000249"; o campo efo_id do ChEMBL
        # pode usar prefixo/separador diferente (ex: "EFO:0000249"), entao
        # buscamos so pelo numero, por seguranca.
        numero = re.sub(r"\D", "", efo_id)
        if not numero:
            continue
        print(f"[INFO][ChEMBL] Buscando indicacoes pelo ID EFO: '{efo_id}' (n. {numero})...")
        try:
            antes = len(todas)
            _chembl_buscar_por_campo("efo_id__icontains", numero, todas)
            print(f"[INFO][ChEMBL]   +{len(todas) - antes} ID(s) ChEMBL via efo_id. Total ate agora: {len(todas)}.")
        except Exception as erro:
            print(f"[AVISO][ChEMBL] Falha ao consultar por efo_id '{efo_id}': {erro}")

    for termo in keywords:
        print(f"[INFO][ChEMBL] Buscando indicacoes por texto para: '{termo}'...")
        try:
            antes = len(todas)
            _chembl_buscar_por_campo("efo_term__icontains", termo, todas)
            _chembl_buscar_por_campo("mesh_heading__icontains", termo, todas)
            print(f"[INFO][ChEMBL]   +{len(todas) - antes} ID(s) ChEMBL via texto. Total: {len(todas)}.")
        except Exception as erro:
            print(f"[AVISO][ChEMBL] Falha ao consultar '{termo}': {erro}")

    if not todas:
        return {}

    # So ordena e resolve nome das top_n candidatas -- nao vale a pena gastar
    # uma chamada de API por nome para milhares de moleculas que vao ser
    # descartadas de qualquer forma.
    ordenadas = sorted(todas.items(), key=lambda kv: _ordinal_fase_clinica(kv[1].get("fase_clinica")), reverse=True)
    top = dict(ordenadas[:top_n])
    print(f"[INFO][ChEMBL] Mantendo as {len(top)} mais relevantes (maior fase clinica) de {len(todas)} encontradas.")

    for chembl_id, info in top.items():
        try:
            resposta = requests.get(f"{CHEMBL_BASE}/molecule/{chembl_id}", params={"format": "json"}, headers=HEADERS_HTTP, timeout=30)
            if resposta.status_code == 200:
                nome = resposta.json().get("pref_name")
                info["nome"] = nome or chembl_id
            else:
                info["nome"] = chembl_id
        except Exception:
            info["nome"] = chembl_id

    return top


def obter_smiles_via_chembl(chembl_id: str) -> str | None:
    try:
        resposta = requests.get(f"{CHEMBL_BASE}/molecule/{chembl_id}", params={"format": "json"}, headers=HEADERS_HTTP, timeout=30)
        if resposta.status_code != 200:
            return None
        dados = resposta.json()
        estrutura = dados.get("molecule_structures") or {}
        return estrutura.get("canonical_smiles")
    except (requests.exceptions.RequestException, ValueError):
        return None


# =========================================================================
# PARTE 1c — DrugCentral: conexão direta ao Postgres público (opcional)
# =========================================================================

def _drugcentral_consultar_termo(cursor, termo: str) -> list[tuple]:
    """
    Tenta a consulta 'enriquecida' (com contagem de produtos aprovados via
    tabela 'product', se ela existir nesse dump do DrugCentral). Se falhar
    (nome de tabela/coluna diferente do esperado), cai para a consulta
    simples -- sem quebrar o script por causa de uma suposicao de schema
    que nao se confirmou.
    Retorna linhas (nome, smiles, match_exato, n_produtos).
    """
    try:
        cursor.execute(
            """
            SELECT s.name, s.smiles,
                   bool_or(LOWER(r.concept_name) = LOWER(%(termo)s)) AS match_exato,
                   COUNT(DISTINCT p.ndc_product_code) AS n_produtos
            FROM structures s
            JOIN omop_relationship r ON r.struct_id = s.id
            LEFT JOIN product p ON p.struct_id = s.id
            WHERE r.relationship_name = 'indication'
              AND r.concept_name ILIKE %(termo_like)s
              AND s.smiles IS NOT NULL
            GROUP BY s.name, s.smiles
            """,
            {"termo": termo, "termo_like": f"%{termo}%"},
        )
        return cursor.fetchall()
    except Exception:
        cursor.connection.rollback()  # limpa a transacao abortada antes de tentar de novo
        cursor.execute(
            """
            SELECT s.name, s.smiles,
                   bool_or(LOWER(r.concept_name) = LOWER(%(termo)s)) AS match_exato,
                   0 AS n_produtos
            FROM structures s
            JOIN omop_relationship r ON r.struct_id = s.id
            WHERE r.relationship_name = 'indication'
              AND r.concept_name ILIKE %(termo_like)s
              AND s.smiles IS NOT NULL
            GROUP BY s.name, s.smiles
            """,
            {"termo": termo, "termo_like": f"%{termo}%"},
        )
        return cursor.fetchall()


def buscar_drugcentral(keywords: list[str], top_n: int) -> dict[str, dict]:
    """
    Retorna {nome_normalizado: {nome, smiles}} consultando o Postgres público
    do DrugCentral. Como o DrugCentral só armazena indicações já aprovadas e
    curadas, não há "fase clínica" para ordenar -- em vez disso, prioriza:
      1. Correspondência EXATA do nome da condição (não só substring) --
         ex: droga indicada para exatamente "Alzheimer's Disease" vem antes
         de uma indicada para "Alzheimer's Disease, Late Onset, ...".
      2. Número de produtos farmacêuticos aprovados com aquele ingrediente
         (proxy de quão estabelecida/mainstream é a droga), quando essa
         informação está disponível no dump.
    Busca TODAS as correspondências primeiro, corta as 'top_n' só no final.
    """
    if not PSYCOPG2_DISPONIVEL:
        print("[AVISO][DrugCentral] Pacote 'psycopg2' não instalado — fonte DrugCentral ignorada. "
              "Instale com: pip install psycopg2-binary")
        return {}

    todas: dict[str, dict] = {}
    try:
        conexao = psycopg2.connect(
            host=DRUGCENTRAL_DBHOST, port=DRUGCENTRAL_DBPORT, dbname=DRUGCENTRAL_DBNAME,
            user=DRUGCENTRAL_DBUSER, password=DRUGCENTRAL_DBPASS, connect_timeout=20,
        )
    except Exception as erro:
        print(f"[AVISO][DrugCentral] Não foi possível conectar ao banco público: {erro}. Fonte ignorada.")
        return {}

    try:
        with conexao.cursor() as cursor:
            for termo in keywords:
                print(f"[INFO][DrugCentral] Buscando indicações para: '{termo}'...")
                linhas = _drugcentral_consultar_termo(cursor, termo)
                for nome, smiles, match_exato, n_produtos in linhas:
                    if not nome or not smiles:
                        continue
                    chave = nome.strip().lower()
                    registro = todas.setdefault(
                        chave, {"nome": nome, "smiles": smiles, "match_exato": False, "n_produtos": 0}
                    )
                    registro["match_exato"] = registro["match_exato"] or bool(match_exato)
                    registro["n_produtos"] = max(registro["n_produtos"], n_produtos or 0)
                print(f"[INFO][DrugCentral]   {len(linhas)} droga(s) encontrada(s) no total para '{termo}'.")
    except Exception as erro:
        print(f"[AVISO][DrugCentral] Erro durante a consulta: {erro}")
    finally:
        conexao.close()

    ordenadas = sorted(
        todas.items(),
        key=lambda kv: (kv[1]["match_exato"], kv[1]["n_produtos"], kv[0]),
        reverse=True,
    )
    top = dict(ordenadas[:top_n])
    print(f"[INFO][DrugCentral] Mantendo as {len(top)} mais relevantes (correspondência exata + nº de produtos aprovados) de {len(todas)} encontradas.")
    return top


# =========================================================================
# PARTE 1d — openFDA: texto de bula (indications_and_usage)
# =========================================================================

def buscar_openfda(keywords: list[str], top_n: int) -> dict[str, dict]:
    """
    Busca no texto de bula (campo 'indications_and_usage') das bulas da FDA
    por menção à doença. Gratuito, sem chave (limite de 40 requisições/min).

    Diferente das outras fontes, o openFDA não tem "fase clínica" nem
    correspondência estruturada a uma doença -- é busca textual livre em
    bula. Como proxy de relevância, agrupa por princípio ativo
    (openfda.substance_name / generic_name) e conta em quantos RÓTULOS DE
    BULA diferentes (produtos comerciais distintos) aquele princípio ativo
    aparece mencionando a doença: mais rótulos = droga mais estabelecida/
    mainstream para essa indicação (parecido com o critério de "nº de
    produtos aprovados" usado no DrugCentral).
    Retorna {nome_normalizado: {nome, n_rotulos}}.
    """
    todas: dict[str, dict] = {}

    for termo in keywords:
        termo_busca = termo.replace('"', "'")
        query = f'indications_and_usage:"{termo_busca}"'
        print(f"[INFO][openFDA] Buscando bulas com indicação para: '{termo}'...")

        skip = 0
        limite_pagina = 100
        limite_seguranca = 2000  # trava de segurança, nunca roda indefinidamente
        total_rotulos_termo = 0

        while skip < limite_seguranca:
            params = {"search": query, "limit": limite_pagina, "skip": skip}
            try:
                resposta = requests.get(OPENFDA_BASE, params=params, headers=HEADERS_HTTP, timeout=60)
            except requests.exceptions.RequestException as erro:
                print(f"[AVISO][openFDA] Falha de rede ao consultar '{termo}': {erro}")
                break

            if resposta.status_code == 404:
                # openFDA usa 404 para "nenhum resultado encontrado" (comportamento documentado da API)
                if skip == 0:
                    print(f"[INFO][openFDA]   0 resultados para '{termo}'.")
                break
            if resposta.status_code != 200:
                print(f"[AVISO][openFDA] HTTP {resposta.status_code} ao consultar '{termo}'. Corpo: {resposta.text[:300]}")
                break

            dados = resposta.json()
            resultados = dados.get("results", [])
            if not resultados:
                break

            for item in resultados:
                info_openfda = item.get("openfda", {}) or {}
                substancias = info_openfda.get("substance_name") or info_openfda.get("generic_name") or []
                for substancia in substancias:
                    chave = (substancia or "").strip().lower()
                    if not chave:
                        continue
                    registro = todas.setdefault(chave, {"nome": substancia.strip(), "n_rotulos": 0})
                    registro["n_rotulos"] += 1

            total_rotulos_termo += len(resultados)
            skip += limite_pagina
            if len(resultados) < limite_pagina:
                break

        print(f"[INFO][openFDA]   {total_rotulos_termo} rótulo(s) de bula encontrados para '{termo}'.")

    if not todas:
        return {}

    ordenadas = sorted(todas.items(), key=lambda kv: kv[1]["n_rotulos"], reverse=True)
    top = dict(ordenadas[:top_n])
    print(f"[INFO][openFDA] Mantendo as {len(top)} mais mencionadas (em mais rótulos de bula) de {len(todas)} substância(s) única(s).")
    return top


# =========================================================================
# PARTE 1e — Fusão das quatro fontes + fallback de SMILES via PubChem
# =========================================================================

def obter_smiles_via_pubchem(nome_droga: str) -> str | None:
    """
    Busca o SMILES de uma droga pelo nome no PubChem.

    IMPORTANTE: o PubChem passou a devolver a propriedade sob a chave
    'ConnectivitySMILES' (ou às vezes 'IsomericSMILES') no JSON, MESMO
    quando se pede 'CanonicalSMILES' na URL — pedir só uma propriedade e
    procurar só por essa mesma chave no retorno faz o valor ser descartado
    silenciosamente mesmo quando a molécula existe. Por isso pedimos as
    três e aceitamos qualquer uma que vier preenchida.
    """
    try:
        propriedades_pedidas = "IsomericSMILES,CanonicalSMILES,ConnectivitySMILES"
        url = f"{PUBCHEM_BASE}/compound/name/{quote(nome_droga)}/property/{propriedades_pedidas}/JSON"
        resposta = requests.get(url, headers=HEADERS_HTTP, timeout=30)
        if resposta.status_code != 200:
            if resposta.status_code != 404:  # 404 = "não encontrado", não é erro de verdade
                print(f"[AVISO][PubChem] HTTP {resposta.status_code} ao buscar '{nome_droga}'. Corpo: {resposta.text[:200]}")
            return None
        dados = resposta.json()
        propriedades = dados.get("PropertyTable", {}).get("Properties", [])
        if not propriedades:
            return None
        registro = propriedades[0]
        for chave in ("IsomericSMILES", "CanonicalSMILES", "ConnectivitySMILES", "SMILES"):
            if registro.get(chave):
                return registro[chave]
    except (requests.exceptions.RequestException, ValueError) as erro:
        print(f"[AVISO][PubChem] Falha ao buscar '{nome_droga}': {erro}")
        return None
    return None


def mesclar_fontes(
    resultado_opentargets: dict[str, dict],
    resultado_chembl: dict[str, dict],
    resultado_drugcentral: dict[str, dict],
    resultado_openfda: dict[str, dict],
) -> dict[str, dict]:
    """
    Mescla as quatro fontes. A chave de deduplicação é o ID ChEMBL quando
    disponível (Open Targets e ChEMBL sempre trazem um) — é uma identidade
    muito mais confiável que comparar nomes em texto, já que "Donepezil"
    (Open Targets) e "DONEPEZIL HYDROCHLORIDE" (ChEMBL) são o mesmo
    medicamento mas não bateriam por string. Só cai para nome normalizado
    quando a entrada não tem ID ChEMBL (caso do DrugCentral e do openFDA).
    Retorna {chave: {nome, chembl_id, smiles, fontes}}.
    """
    mesclado: dict[str, dict] = {}

    for origem_nome, dicionario in (
        ("opentargets", resultado_opentargets),
        ("chembl", resultado_chembl),
        ("drugcentral", resultado_drugcentral),
        ("openfda", resultado_openfda),
    ):
        for info in dicionario.values():
            chembl_id = (info.get("chembl_id") or "").strip().upper()
            nome = (info.get("nome") or "").strip()
            chave = chembl_id or nome.lower()
            if not chave:
                continue

            item = mesclado.setdefault(chave, {"nome": nome or chave, "fontes": set()})
            item["fontes"].add(origem_nome)
            if chembl_id:
                item.setdefault("chembl_id", chembl_id)
            if info.get("smiles"):
                item.setdefault("smiles", info["smiles"])
            # Prefere nomes "bonitos" (Open Targets/DrugCentral tendem a vir
            # com capitalização normal; ChEMBL às vezes só tem o ID como nome).
            if nome and nome.upper() != nome and item["nome"].upper() == item["nome"]:
                item["nome"] = nome

    # Sem corte aqui: cada fonte ja trouxe so as suas top_n mais relevantes,
    # entao o total final e naturalmente pequeno (no maximo a soma das 3).
    # Ainda assim ordena por numero de fontes, para o CSV de referencia
    # mostrar primeiro as drogas confirmadas por mais de uma base.
    mesclado = dict(sorted(mesclado.items(), key=lambda kv: len(kv[1]["fontes"]), reverse=True))
    return mesclado


def obter_smiles_das_drogas(drogas: dict[str, dict]) -> list[tuple[str, str, set]]:
    """
    Garante o SMILES de cada droga mesclada: usa o já obtido (DrugCentral),
    senão tenta via ChEMBL (se houver chembl_id), senão via PubChem (nome).
    Retorna lista de (nome, smiles, fontes) só para as que têm SMILES válido.
    """
    resultado = []
    total = len(drogas)
    for i, (_chave, info) in enumerate(drogas.items(), start=1):
        nome = info["nome"]
        print(f"[INFO] ({i}/{total}) Garantindo SMILES de '{nome}' (fontes: {', '.join(sorted(info['fontes']))})...")

        smiles = info.get("smiles")
        if not smiles and info.get("chembl_id"):
            smiles = obter_smiles_via_chembl(info["chembl_id"])
        if not smiles:
            smiles = obter_smiles_via_pubchem(nome)

        if not smiles:
            print(f"[AVISO]   SMILES não encontrado para '{nome}'. Droga ignorada na referência.")
            continue
        if Chem.MolFromSmiles(smiles) is None:
            print(f"[AVISO]   SMILES inválido para '{nome}' (RDKit não interpretou). Ignorada.")
            continue

        resultado.append((nome, smiles, info["fontes"]))
        time.sleep(0.1)

    print(f"[INFO] SMILES obtido com sucesso para {len(resultado)} de {total} droga(s) de referência.")
    return resultado


# =========================================================================
# PARTE 2 — ADMET (Deep-PK) para as drogas de referência (opcional)
# =========================================================================

REGEX_CHAVE = re.compile(r"^\[(?P<categoria>[^/\]]+)/(?P<propriedade>.+)\]\s+(?P<tipo>Predictions|Probability|Interpretation)$")

CATEGORIA_PREFIXOS = {
    "General Properties": "Geral",
    "Absorption": "Absorcao",
    "Distribution": "Distribuicao",
    "Metabolism": "Metabolismo",
    "Excretion": "Excrecao",
    "Toxicity": "Toxicidade",
}


def slugify(texto: str) -> str:
    texto = re.sub(r"[^\w]+", "_", texto, flags=re.UNICODE)
    return texto.strip("_")


def limpar_texto(texto) -> str:
    texto = str(texto)
    texto = re.sub(r"<br\s*/?>", " ", texto, flags=re.IGNORECASE)
    texto = texto.replace("&nbsp;", " ")
    return re.sub(r"\s+", " ", texto).strip()


def _decodificar_json(texto: str):
    dado = json.loads(texto)
    if isinstance(dado, str):
        dado = json.loads(dado)
    return dado


def submeter_job_deeppk(smiles_lote: list[str], pred_type: str = "admet") -> str:
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False, encoding="utf-8") as tmp:
        tmp.write("\n".join(smiles_lote))
        caminho_tmp = tmp.name
    try:
        with open(caminho_tmp, "rb") as f:
            resposta = requests.post(
                DEEPPK_API_URL, files={"smiles_file": f}, data={"pred_type": pred_type}, timeout=120
            )
    finally:
        os.remove(caminho_tmp)
    resposta.raise_for_status()
    payload = _decodificar_json(resposta.text)
    job_id = payload.get("job_id") if isinstance(payload, dict) else None
    if not job_id:
        raise RuntimeError(f"Resposta inesperada do Deep-PK ao submeter job: {payload}")
    return job_id


def aguardar_resultado_deeppk(job_id: str) -> dict:
    inicio = time.time()
    while True:
        resposta = requests.get(DEEPPK_API_URL, data={"job_id": job_id}, timeout=60)
        resposta.raise_for_status()
        payload = _decodificar_json(resposta.text)
        if isinstance(payload, dict) and payload.get("status") == "running":
            if time.time() - inicio > DEEPPK_MAX_WAIT:
                raise TimeoutError(f"Tempo máximo excedido aguardando job {job_id} do Deep-PK.")
            time.sleep(DEEPPK_POLL_INTERVAL)
            continue
        if not isinstance(payload, dict):
            raise RuntimeError(f"Resposta inesperada do Deep-PK para o job {job_id}.")
        return payload


def resumir_molecula_deeppk(dados_molecula: dict) -> dict:
    valores: dict[str, str] = {}
    for chave, valor in dados_molecula.items():
        m = REGEX_CHAVE.match(chave)
        if not m or m.group("tipo") != "Predictions":
            continue
        categoria = m.group("categoria").strip()
        propriedade = m.group("propriedade").strip()
        prefixo = CATEGORIA_PREFIXOS.get(categoria)
        if not prefixo:
            continue
        nome_coluna = f"{prefixo}_{slugify(propriedade)}"

        chave_interp = f"[{categoria}/{propriedade}] Interpretation"
        interpretacao = dados_molecula.get(chave_interp, "-")
        tem_interp = interpretacao is not None and str(interpretacao).strip() not in ("-", "None", "")
        if tem_interp:
            interpretacao = limpar_texto(interpretacao)
            valores[nome_coluna] = f"{valor} ({interpretacao})"
        else:
            valores[nome_coluna] = f"{valor}"
    return valores


def consultar_admet_referencia(nomes_smiles: list[tuple[str, str]], pred_type: str = "admet") -> dict[str, dict]:
    resultado_por_nome: dict[str, dict] = {}
    total = len(nomes_smiles)
    indice = 0
    while indice < total:
        lote = nomes_smiles[indice: indice + DEEPPK_BATCH_SIZE]
        smiles_lote = [s for _, s in lote]
        print(f"[INFO] Consultando ADMET (Deep-PK) para drogas de referência {indice + 1}-{indice + len(lote)} de {total}...")
        try:
            job_id = submeter_job_deeppk(smiles_lote, pred_type)
            resultados = aguardar_resultado_deeppk(job_id)
        except Exception as erro:
            print(f"[AVISO]   Falha ao consultar Deep-PK para este lote: {erro}. Lote ignorado.")
            indice += len(lote)
            continue

        for indice_local_str, dados_molecula in resultados.items():
            if not indice_local_str.isdigit() or not isinstance(dados_molecula, dict):
                continue
            idx_local = int(indice_local_str)
            if idx_local >= len(lote):
                continue
            nome_droga = lote[idx_local][0]
            resultado_por_nome[nome_droga] = resumir_molecula_deeppk(dados_molecula)

        indice += len(lote)

    return resultado_por_nome


# =========================================================================
# PARTE 4 — Leitura do farmaco.csv
# =========================================================================

def detectar_delimitador(caminho: str) -> str:
    with open(caminho, "r", encoding="utf-8-sig", newline="") as f:
        amostra = f.read(4096)
    try:
        return csv.Sniffer().sniff(amostra, delimiters=";,\t").delimiter
    except csv.Error:
        return CSV_DELIMITADOR_PADRAO


def ler_farmaco_csv(caminho: str) -> tuple[list[str], list[dict]]:
    delimitador = detectar_delimitador(caminho)
    with open(caminho, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter=delimitador)
        colunas = reader.fieldnames or []
        linhas = list(reader)
    return colunas, linhas


# =========================================================================
# PARTE 5 — Cálculo dos pesos por característica
# =========================================================================

def _para_float(valor) -> float | None:
    if valor is None or valor == "":
        return None
    try:
        return float(valor)
    except (TypeError, ValueError):
        return None


def _para_bool(valor) -> bool | None:
    if isinstance(valor, bool):
        return valor
    if valor is None:
        return None
    texto = str(valor).strip().lower()
    if texto in ("true", "1", "sim", "yes"):
        return True
    if texto in ("false", "0", "nao", "não", "no"):
        return False
    return None


def _categoria_de(valor) -> str | None:
    if valor is None:
        return None
    texto = str(valor).strip()
    if not texto:
        return None
    return texto.split(" (")[0].strip()


def calcular_pesos(linhas_referencia: list[dict], colunas_alvo: list[str]) -> dict[str, dict]:
    pesos: dict[str, dict] = {}

    for coluna in colunas_alvo:
        if coluna in COLUNAS_NAO_CARACTERISTICA or coluna in COLUNAS_IGNORADAS:
            continue

        if coluna in COLUNAS_NUMERICAS:
            valores = [_para_float(l.get(coluna)) for l in linhas_referencia]
            valores = [v for v in valores if v is not None]
            if len(valores) < 2:
                continue
            media = statistics.mean(valores)
            desvio = statistics.pstdev(valores)
            cv = (desvio / abs(media)) if media not in (0, None) else 1.0
            peso = 1.0 / (1.0 + cv)
            pesos[coluna] = {
                "tipo": "numerico", "peso": round(peso, 4),
                "media": round(media, 3), "desvio": round(desvio, 3), "n": len(valores),
            }

        elif coluna in COLUNAS_BOOLEANAS:
            valores = [_para_bool(l.get(coluna)) for l in linhas_referencia]
            valores = [v for v in valores if v is not None]
            if not valores:
                continue
            contagem = Counter(valores)
            valor_moda, freq_moda = contagem.most_common(1)[0]
            peso = freq_moda / len(valores)
            pesos[coluna] = {"tipo": "booleano", "peso": round(peso, 4), "moda": valor_moda, "n": len(valores)}

        else:
            valores = [_categoria_de(l.get(coluna)) for l in linhas_referencia]
            valores = [v for v in valores if v]
            if not valores:
                continue
            contagem = Counter(valores)
            valor_moda, freq_moda = contagem.most_common(1)[0]
            peso = freq_moda / len(valores)
            pesos[coluna] = {"tipo": "categorico", "peso": round(peso, 4), "moda": valor_moda, "n": len(valores)}

    return pesos


# =========================================================================
# PARTE 6 — Pontuação e ranking dos compostos do farmaco.csv
# =========================================================================

def pontuar_composto(linha: dict, pesos: dict[str, dict]) -> tuple[float, dict]:
    soma_ponderada = 0.0
    soma_pesos_usados = 0.0
    detalhes = {}

    for coluna, info in pesos.items():
        peso = info["peso"]
        valor_bruto = linha.get(coluna)

        if info["tipo"] == "numerico":
            valor = _para_float(valor_bruto)
            if valor is None:
                continue
            media, desvio = info["media"], info["desvio"]
            if desvio and desvio > 0:
                z = (valor - media) / desvio
                similaridade = math.exp(-0.5 * (z ** 2))
            else:
                similaridade = 1.0 if valor == media else 0.0

        elif info["tipo"] == "booleano":
            valor = _para_bool(valor_bruto)
            if valor is None:
                continue
            similaridade = 1.0 if valor == info["moda"] else 0.0

        else:
            valor = _categoria_de(valor_bruto)
            if not valor:
                continue
            similaridade = 1.0 if valor == info["moda"] else 0.0

        soma_ponderada += peso * similaridade
        soma_pesos_usados += peso
        detalhes[coluna] = round(similaridade, 3)

    if soma_pesos_usados == 0:
        return 0.0, detalhes

    pontuacao = (soma_ponderada / soma_pesos_usados) * 100
    return round(pontuacao, 2), detalhes


def gerar_ranking(linhas_farmaco: list[dict], pesos: dict[str, dict]) -> list[dict]:
    ranking = []
    for linha in linhas_farmaco:
        pontuacao, _detalhes = pontuar_composto(linha, pesos)
        ranking.append({
            "Nome": linha.get("Nome", ""),
            "SMILES": linha.get("SMILES", ""),
            "Pontuacao_Final": pontuacao,
        })
    ranking.sort(key=lambda r: r["Pontuacao_Final"], reverse=True)
    for i, item in enumerate(ranking, start=1):
        item["Ranking"] = i
    return ranking


def salvar_ranking(ranking: list[dict], caminho_saida: str, delimitador: str = CSV_DELIMITADOR_PADRAO) -> None:
    colunas = ["Ranking", "Nome", "SMILES", "Pontuacao_Final"]
    with open(caminho_saida, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=colunas, delimiter=delimitador)
        writer.writeheader()
        for item in ranking:
            writer.writerow({c: item.get(c) for c in colunas})
    print(f"[INFO] Ranking salvo em: {caminho_saida}")


def salvar_pesos(pesos: dict[str, dict], caminho_saida: str) -> None:
    with open(caminho_saida, "w", encoding="utf-8") as f:
        json.dump(pesos, f, ensure_ascii=False, indent=2, default=str)
    print(f"[INFO] Pesos das características salvos em: {caminho_saida}")


def salvar_referencia(drogas: list[tuple[str, str, set]], caminho_saida: str) -> None:
    with open(caminho_saida, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter=CSV_DELIMITADOR_PADRAO)
        writer.writerow(["Nome", "SMILES", "Fontes"])
        for nome, smiles, fontes in drogas:
            writer.writerow([nome, smiles, "+".join(sorted(fontes))])
    print(f"[INFO] Drogas de referência (com fontes) salvas em: {caminho_saida}")


# =========================================================================
# MAIN
# =========================================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Consulta Open Targets + ChEMBL + DrugCentral simultaneamente para achar drogas "
            "indicadas para uma doença, deriva pesos por característica, e ranqueia um farmaco.csv."
        )
    )
    parser.add_argument("--doenca", required=True, nargs="+", help="Palavra(s)-chave da doença.")
    parser.add_argument("--farmaco", default="farmaco.csv", help="Caminho do farmaco.csv de entrada (padrão: farmaco.csv).")
    parser.add_argument("--output", default="ranking.csv", help="Caminho do CSV de ranking de saída (padrão: ranking.csv).")
    parser.add_argument("--weights-output", default=None, help="Caminho para salvar os pesos calculados em JSON.")
    parser.add_argument("--referencia-output", default=None, help="Caminho para salvar a lista de drogas de referência usadas.")
    parser.add_argument("--top-por-fonte", type=int, default=10, help="Quantas drogas manter de CADA fonte (Open Targets, ChEMBL, DrugCentral), priorizando maior fase clínica/aprovação. Cada fonte é consultada por inteiro, sem limite, e só depois cortada. Padrão: 10.")
    parser.add_argument("--no-opentargets", action="store_true", help="Não consulta o Open Targets.")
    parser.add_argument("--no-chembl", action="store_true", help="Não consulta o ChEMBL.")
    parser.add_argument("--no-drugcentral", action="store_true", help="Não consulta o DrugCentral.")
    parser.add_argument("--no-openfda", action="store_true", help="Não consulta o openFDA.")
    parser.add_argument("--no-admet", action="store_true", help="Não consulta o Deep-PK para as drogas de referência.")
    parser.add_argument("--pred-type", default="admet", choices=["absorption", "distribution", "metabolism", "excretion", "admet"], help="Tipo de predição do Deep-PK.")
    args = parser.parse_args()

    if args.no_opentargets and args.no_chembl and args.no_drugcentral and args.no_openfda:
        print("[ERRO] Todas as fontes foram desativadas — não há como buscar drogas de referência.", file=sys.stderr)
        sys.exit(1)

    if not Path(args.farmaco).is_file():
        print(f"[ERRO] Arquivo não encontrado: {args.farmaco}", file=sys.stderr)
        sys.exit(1)

    # ---- Etapa 1: buscar drogas de referência nas 4 fontes ----
    resultado_ot, efo_ids = ({}, []) if args.no_opentargets else buscar_open_targets(args.doenca, args.top_por_fonte)
    resultado_chembl = {} if args.no_chembl else buscar_chembl_indicacao(args.doenca, args.top_por_fonte, efo_ids=efo_ids)
    resultado_dc = {} if args.no_drugcentral else buscar_drugcentral(args.doenca, args.top_por_fonte)
    resultado_fda = {} if args.no_openfda else buscar_openfda(args.doenca, args.top_por_fonte)

    drogas_mescladas = mesclar_fontes(resultado_ot, resultado_chembl, resultado_dc, resultado_fda)
    print(f"[INFO] Total de drogas únicas mescladas das 4 fontes: {len(drogas_mescladas)}")
    if not drogas_mescladas:
        print("[ERRO] Nenhuma droga de referência encontrada em nenhuma fonte. Encerrando.", file=sys.stderr)
        sys.exit(1)

    # ---- Etapa 2: garantir SMILES de cada droga mesclada ----
    nomes_smiles_fontes = obter_smiles_das_drogas(drogas_mescladas)
    if not nomes_smiles_fontes:
        print("[ERRO] Não foi possível obter SMILES de nenhuma droga de referência. Encerrando.", file=sys.stderr)
        sys.exit(1)

    caminho_referencia = args.referencia_output or (str(Path(args.output).with_suffix("")) + ".referencia.csv")
    salvar_referencia(nomes_smiles_fontes, caminho_referencia)

    # ---- Etapa 3: descritores locais das drogas de referência ----
    # Reaproveita a mesma função do farmaco_completo.py usada para gerar o
    # próprio farmaco.csv, garantindo colunas e regras idênticas.
    linhas_referencia = processar_lista_druglikeness(
        [(nome, smiles) for nome, smiles, _fontes in nomes_smiles_fontes]
    )
    print(f"[INFO] Descritores de drogabilidade calculados para {len(linhas_referencia)} droga(s) de referência.")

    # ---- Etapa 4: colunas do farmaco.csv ----
    colunas_farmaco, linhas_farmaco = ler_farmaco_csv(args.farmaco)
    print(f"[INFO] '{args.farmaco}' lido: {len(linhas_farmaco)} composto(s), {len(colunas_farmaco)} coluna(s).")

    colunas_admet_no_farmaco = [
        c for c in colunas_farmaco
        if c not in COLUNAS_NAO_CARACTERISTICA and c not in COLUNAS_NUMERICAS
        and c not in COLUNAS_BOOLEANAS and c not in COLUNAS_IGNORADAS
    ]

    # ---- Etapa 5 (opcional): ADMET das drogas de referência via Deep-PK ----
    if not args.no_admet and colunas_admet_no_farmaco:
        nomes_smiles = [(nome, smiles) for nome, smiles, _f in nomes_smiles_fontes]
        admet_referencia = consultar_admet_referencia(nomes_smiles, args.pred_type)
        for linha in linhas_referencia:
            valores_admet = admet_referencia.get(linha["Nome"], {})
            for coluna in colunas_admet_no_farmaco:
                linha[coluna] = valores_admet.get(coluna, "")
    elif colunas_admet_no_farmaco:
        print("[INFO] --no-admet ativado: colunas ADMET do farmaco.csv não serão pesadas.")

    # ---- Etapa 6: cálculo dos pesos ----
    pesos = calcular_pesos(linhas_referencia, colunas_farmaco)
    if not pesos:
        print("[ERRO] Não foi possível calcular nenhum peso (dados de referência insuficientes). Encerrando.", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Pesos calculados para {len(pesos)} característica(s):")
    for coluna, info in sorted(pesos.items(), key=lambda kv: kv[1]["peso"], reverse=True):
        print(f"         {coluna:35s} peso={info['peso']:.3f}  (tipo={info['tipo']})")

    caminho_pesos = args.weights_output or (str(Path(args.output).with_suffix("")) + ".pesos.json")
    salvar_pesos(pesos, caminho_pesos)

    # ---- Etapa 7: ranking do farmaco.csv com base nos pesos ----
    ranking = gerar_ranking(linhas_farmaco, pesos)
    salvar_ranking(ranking, args.output)

    print("\n[INFO] Top 10 do ranking:")
    for item in ranking[:10]:
        print(f"  {item['Ranking']:>3}. {item['Nome']:<30s} {item['Pontuacao_Final']:.2f}")


if __name__ == "__main__":
    main()
