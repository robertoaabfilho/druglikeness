#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
configuracao.py

Leitor simples de 'config.txt': um arquivo de texto com linhas no formato
    chave = valor
(uma por linha; linhas vazias ou iniciadas com '#' são ignoradas).

Usado por main.py e pelos demais scripts do pipeline para permitir
configuração central num só lugar, em vez de repetir os mesmos parâmetros
como argumento de linha de comando em cada script. Argumentos passados
explicitamente na linha de comando sempre têm prioridade sobre o config.txt
-- o config.txt só preenche o que não foi passado.

Exemplo de config.txt:

    # Sintoma(s) a pesquisar nas bases de referência
    sintomas = neuropathic pain, peripheral neuropathy

    # farmaco_completo.py
    farmaco_input = input.txt
    farmaco_output = farmaco.csv
    pred_type = admet
    batch_size = 5

    # ranking_multibase.py
    ranking_output = ranking.csv
    top_por_fonte = 10
    no_admet = false

    # validar_enriquecimento.py
    n_decoys = 50
"""

from pathlib import Path

CAMINHO_PADRAO = "config.txt"


def carregar_config(caminho: str = CAMINHO_PADRAO) -> dict[str, str]:
    """
    Lê o config.txt e retorna um dicionário {chave: valor}, ambos como
    string (a conversão de tipo é feita pelos getters abaixo). Se o
    arquivo não existir, retorna um dicionário vazio silenciosamente --
    config.txt é opcional, todos os scripts continuam funcionando só com
    argumentos de linha de comando.
    """
    config: dict[str, str] = {}
    caminho_path = Path(caminho)
    if not caminho_path.is_file():
        return config

    with open(caminho_path, "r", encoding="utf-8-sig") as f:
        for numero_linha, linha in enumerate(f, start=1):
            linha = linha.strip()
            if not linha or linha.startswith("#"):
                continue
            if "=" not in linha:
                print(f"[AVISO][config] Linha {numero_linha} de '{caminho}' ignorada (sem '='): {linha!r}")
                continue
            chave, _, valor = linha.partition("=")
            config[chave.strip().lower()] = valor.strip()

    return config


def get_str(config: dict, chave: str, default: str | None = None) -> str | None:
    valor = config.get(chave)
    return valor if valor not in (None, "") else default


def get_int(config: dict, chave: str, default: int | None = None) -> int | None:
    valor = config.get(chave)
    if valor in (None, ""):
        return default
    try:
        return int(valor)
    except ValueError:
        print(f"[AVISO][config] Valor inválido para '{chave}': {valor!r} (esperado inteiro). Usando padrão: {default}")
        return default


def get_float(config: dict, chave: str, default: float | None = None) -> float | None:
    valor = config.get(chave)
    if valor in (None, ""):
        return default
    try:
        return float(valor)
    except ValueError:
        print(f"[AVISO][config] Valor inválido para '{chave}': {valor!r} (esperado número). Usando padrão: {default}")
        return default


def get_bool(config: dict, chave: str, default: bool = False) -> bool:
    valor = config.get(chave)
    if valor in (None, ""):
        return default
    return valor.strip().lower() in ("1", "true", "sim", "yes", "on", "verdadeiro")


def get_list(config: dict, chave: str, default: list | None = None) -> list[str]:
    """Lê um valor separado por vírgulas (ex: 'sintomas = a, b, c') como lista de strings."""
    valor = config.get(chave)
    if valor in (None, ""):
        return list(default) if default else []
    return [item.strip() for item in valor.split(",") if item.strip()]
