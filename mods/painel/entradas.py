# -*- coding: utf-8 -*-
"""
entradas.py - receber o que o pesquisador fornece, e conferir antes de rodar.

Os tres arquivos de entrada sao os dele: a lista de ligantes, a lista de
alvos e a lista de sintomas. O pipeline nao deduz nenhuma das tres.

**A conferencia e a MESMA da etapa 1.** Quem le e valida aqui e
`mods.evidencias.entrada`, os mesmos leitores que `coletar_evidencias.py`
chama. Uma validacao propria do painel divergiria da do pipeline, e o painel
aprovaria entrada que a etapa recusa -- ou pior, o contrario.

**O upload diz o que substituiu.** Trocar `input/ligantes.tsv` em silencio
deixaria o pesquisador sem saber qual base gerou o resultado que ele tem na
tela, e o manifesto registra o MD5 da base, nao o nome do arquivo que ele
subiu.

**A gravacao e atomica.** Upload interrompido no meio nao pode destruir a
base anterior: ela e entrada de uma execucao de horas.
"""
import hashlib
from pathlib import Path

import streamlit as st

from mods import caminhos
from mods.evidencias import atomico, entrada

# Nome canonico de cada papel em `input/`. O upload pode ter qualquer nome; o
# que a configuracao aponta e este.
PAPEIS = {
    "ligantes": "ligantes.tsv",
    "alvos": "alvos.txt",
    "sintomas": "sintomas.txt",
}

# Quantas mensagens de problema a tela mostra. Dezessete mil linhas sem nome
# produzem dezessete mil problemas; o que informa e as primeiras mais a
# contagem.
TETO_MENSAGENS = 10


@st.cache_data(max_entries=12, show_spinner=False)
def _md5_de(caminho_str: str, mtime: float, tamanho: int) -> str:
    soma = hashlib.md5()
    with open(caminho_str, "rb") as arq:
        for bloco in iter(lambda: arq.read(1 << 20), b""):
            soma.update(bloco)
    return soma.hexdigest()


def _md5(caminho: Path) -> str | None:
    """MD5 do arquivo, cacheado por `(caminho, mtime, tamanho)`.

    `input/ligantes.tsv` tem 79 MB nesta arvore, e `estado()` e chamada a cada
    redesenho da aba Operar -- qualquer widget, qualquer clique. Sem o cache
    eram ~290 ms por rerun, com o disco quente, contra o que o docstring de
    `estado()` promete."""
    if not caminho.is_file():
        return None
    info = caminho.stat()
    return _md5_de(str(caminho), info.st_mtime, info.st_size)


def caminho_do_papel(papel: str) -> Path:
    if papel not in PAPEIS:
        raise ValueError(
            f"papel desconhecido: {papel!r}; os papeis sao "
            f"{', '.join(PAPEIS)}")
    return caminhos.INPUT / PAPEIS[papel]


def receber(envio, papel: str) -> dict:
    """Grava o arquivo enviado em `input/` e relata o que mudou.

    Devolve `destino`, `substituiu`, `md5_anterior`, `md5`, `bytes` e a
    `conferencia` do arquivo novo -- a mesma que a etapa 1 faria."""
    destino = caminho_do_papel(papel)
    anterior = _md5(destino)
    destino.parent.mkdir(parents=True, exist_ok=True)
    dados = envio.getvalue()
    tmp = atomico.temporario(destino)
    try:
        tmp.write_bytes(dados)
        atomico.publicar(tmp, destino)
    finally:
        tmp.unlink(missing_ok=True)  # so o temporario, nunca o arquivo do usuario
    return {"destino": destino, "substituiu": anterior is not None,
            "md5_anterior": anterior, "md5": _md5(destino),
            "bytes": len(dados), "enviado_como": envio.name,
            "conferencia": conferir(destino, papel)}


def conferir(caminho, papel: str) -> dict:
    """Le a entrada com os leitores da etapa 1 e resume o que achou.

    Nunca levanta por conteudo: arquivo ausente, coluna que falta e molecula
    que o RDKit nao le saem como contagem e mensagem. Levantar aqui tiraria
    da tela justamente a informacao que ela existe para dar."""
    caminho = Path(caminho)
    if papel not in PAPEIS:
        raise ValueError(
            f"papel desconhecido: {papel!r}; os papeis sao "
            f"{', '.join(PAPEIS)}")
    if not caminho.is_file():
        return {"aceitos": 0, "erros": 1, "avisos": 0, "problemas": 1,
                "mensagens": [f"arquivo nao encontrado: {caminho}"]}
    try:
        if papel == "sintomas":
            itens, problemas = entrada.ler_sintomas(caminho), []
        elif papel == "alvos":
            itens, problemas = entrada.ler_alvos(caminho)
        else:
            itens, problemas = entrada.ler_ligantes(caminho)
    except OSError as erro:
        return {"aceitos": 0, "erros": 1, "avisos": 0, "problemas": 1,
                "mensagens": [str(erro)]}
    erros = [p for p in problemas if p.gravidade == "erro"]
    return {
        "aceitos": len(itens),
        "erros": len(erros),
        "avisos": len(problemas) - len(erros),
        "problemas": len(problemas),
        "mensagens": [f"linha {p.linha}, {p.campo}: {p.mensagem}"
                      for p in problemas[:TETO_MENSAGENS]],
    }


def estado() -> dict:
    """O que existe em `input/` hoje, por papel: `existe`, `caminho`, `md5`,
    `bytes` e `mtime`. Sem conferir -- conferir le o arquivo inteiro, e a
    tela nao faz isso a cada redesenho."""
    saida = {}
    for papel in PAPEIS:
        caminho = caminho_do_papel(papel)
        existe = caminho.is_file()
        saida[papel] = {
            "existe": existe, "caminho": caminho,
            "md5": _md5(caminho) if existe else None,
            "bytes": caminho.stat().st_size if existe else None,
            "mtime": caminho.stat().st_mtime if existe else None,
        }
    return saida
