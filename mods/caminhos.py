# -*- coding: utf-8 -*-
"""
caminhos.py - onde fica cada coisa na pasta do projeto

    <projeto>/
        main.py       orquestrador das tres etapas
        painel.py     interface: operar e analisar
        input/        o que VOCE fornece: ligantes.tsv, alvos.txt,
                      sintomas.txt, evidencias_chembl.csv, procedencia.json
        mods/         codigo: as tres etapas, as CLIs e a infraestrutura
        modelo/       pesos .pt e target_labels.tsv (fora do git)
        json/         config.json e cache_rede/<versao-dos-servicos>/
        tsv/          o que a MAQUINA gerou, uma pasta por etapa:
                      etapa1/, portao/, ranking/ -- cada tabela com o seu
                      manifesto.json ao lado -- e .execucao/, o estado da
                      execucao que o painel disparou
        dados/        dado de referencia versionado, com licenca propria
        exemplos/     fixtures de validacao
        tests/        a suite; nenhum teste usa rede

A regra que desempata `input/` de `tsv/`, porque ligantes e alvos tambem sao
TSV: **procedencia, nao extensao**. O que voce forneceu fica em `input/`; o que
uma etapa escreveu fica em `tsv/`.

Sem efeitos colaterais: este modulo so define caminhos. Nenhuma pasta e criada
no import -- quem cria e a etapa que vai escrever, no momento de publicar.
"""
from pathlib import Path

MODS = Path(__file__).resolve().parent
RAIZ = MODS.parent

INPUT = RAIZ / "input"
MODELO = RAIZ / "modelo"
JSON = RAIZ / "json"
TSV = RAIZ / "tsv"
DADOS = RAIZ / "dados"
EXEMPLOS = RAIZ / "exemplos"
TESTS = RAIZ / "tests"

# Estado da execucao disparada pelo painel: PID, log, ritmo e o resultado da
# ultima. Fica em tsv/ porque e coisa que a MAQUINA escreveu, e o ponto no
# nome mantem a pasta fora do caminho de quem procura tabela.
EXECUCAO = TSV / ".execucao"

CACHE = JSON / "cache_rede"
CONFIG = JSON / "config.json"
PROCEDENCIA = INPUT / "procedencia.json"

# As tres etapas, na ordem em que rodam. O nome e o da pasta em tsv/ e o `id`
# da tabela ETAPAS do main.py: um nome so para a etapa em todo o projeto.
ETAPAS = ("etapa1", "portao", "ranking")


def etapa(nome: str) -> Path:
    """Pasta de saida de uma etapa: tsv/<nome>. Nome desconhecido e erro, nunca
    uma pasta inventada em silencio -- a mesma regra das flags das CLIs."""
    if nome not in ETAPAS:
        raise ValueError(
            f"etapa desconhecida: {nome!r}; as etapas sao {', '.join(ETAPAS)}")
    return TSV / nome


def tsv(*partes: str) -> Path:
    """Caminho dentro de tsv/: tsv('etapa1', 'evidencias.tsv')."""
    return TSV.joinpath(*partes)


def cache_rede(versao: str) -> Path:
    """Cache de uma versao dos servicos: json/cache_rede/<versao>. Versao nova
    nunca le o cache velho, e por isso a versao esta no caminho."""
    return CACHE / versao
