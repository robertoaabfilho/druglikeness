# -*- coding: utf-8 -*-
"""
leitura.py - as tabelas publicadas, lidas do disco, com cache.

O painel nao recebe resultado de ninguem: ele le o sistema de arquivos. E
isso que faz as duas portas conviverem -- execucao disparada no terminal
aparece no painel, execucao disparada no painel e lida por um script -- e o
sistema de arquivos e o UNICO acoplamento entre elas.

**O cache tem chave `(caminho, mtime)`.** Nao basta o caminho: a etapa
republica o arquivo no mesmo lugar, e um cache por caminho so entregaria a
tabela da execucao anterior como se fosse a nova. O `mtime` e lido do disco a
cada chamada; so a leitura e cacheada.

**Quem le e `ler_tsv`, nunca `pandas.read_csv`.** `EscritorTSV` escapa
tabulacao, nova linha e barra invertida dentro do campo; `ler_tsv` e o leitor
que reverte isso. `read_csv`, mesmo com `QUOTE_NONE`, devolveria `\\t` literal
dentro do nome -- e nomes com aspas e tabulacao existem de verdade na base.

**Toda leitura tem teto.** `ranking.tsv` da execucao real tem milhoes de
linhas; carregar tudo num DataFrame trava a aba. Quem chama passa `limite`, e
a tela diz quando a tabela esta truncada. O TSV continua inteiro no disco.
"""
import json
from pathlib import Path

import pandas as pd
import streamlit as st

from mods import caminhos
from mods.portao.esquema import ler_tsv

# Teto de linhas que uma tela carrega por padrao. Acima disso a tela avisa.
# 50.000 linhas de ranking cabem em memoria e sao mais do que alguem le; o
# numero existe para que nenhuma tela possa, por descuido, pedir a tabela
# inteira de 23 milhoes de pares.
LIMITE_PADRAO = 50_000


class TabelaInvalida(Exception):
    """O arquivo existe e nao e um TSV do projeto: cabecalho `#` ausente,
    linha com largura diferente do cabecalho, ou zero bytes.

    E um erro tipado e nao um traceback porque quem le e uma tela: o
    `ValueError` do `ler_tsv` subia como stack trace, e o Review Focus 5 diz
    que nenhum estado de tabela pode fazer isso. Tambem nao vira "tabela
    vazia": arquivo truncado apresentado como "rodou e nada qualificou"
    esconderia o truncamento."""


def _etapas() -> dict:
    """A tabela ETAPAS do orquestrador, por id. Importada aqui dentro para
    que `import mods.painel.leitura` nao arraste o main -- e para que a
    sequencia das etapas tenha UMA fonte, que e o main.py."""
    import main
    return {e["id"]: e for e in main.ETAPAS}


@st.cache_data(show_spinner="Lendo a tabela...", max_entries=12)
def _ler(caminho_str: str, mtime: float, limite: int | None) -> pd.DataFrame:
    if Path(caminho_str).stat().st_size == 0:
        raise TabelaInvalida(f"{caminho_str}: arquivo vazio (zero bytes), "
                             "nem o cabecalho foi escrito")
    linhas = []
    colunas: list[str] = []
    try:
        for i, linha in enumerate(ler_tsv(caminho_str)):
            if not colunas:
                colunas = list(linha)
            if limite is not None and i >= limite:
                break
            linhas.append(linha)
    except ValueError as erro:
        raise TabelaInvalida(str(erro)) from erro
    if not colunas:  # so cabecalho: pega os nomes sem gerar nenhuma linha
        colunas = _colunas(caminho_str)
    return pd.DataFrame(linhas, columns=colunas, dtype=str)


def _colunas(caminho_str: str) -> list[str]:
    """O cabecalho de um TSV do projeto, sem ler o resto."""
    with open(caminho_str, encoding="utf-8", newline="\n") as arq:
        primeira = arq.readline().rstrip("\n")
    return primeira[1:].split("\t") if primeira.startswith("#") else []


def ler_tabela(caminho, limite: int | None = LIMITE_PADRAO):
    """DataFrame de um TSV do projeto; `None` se o arquivo nao existe.

    Tabela vazia nao e ausencia: um arquivo com cabecalho e nada mais e
    resultado legitimo (nenhum ligante ranqueou), e vem como quadro vazio COM
    as colunas, para a tela poder dizer qual tabela esta vazia."""
    caminho = Path(caminho)
    if not caminho.is_file():
        return None
    return _ler(str(caminho), caminho.stat().st_mtime, limite)


@st.cache_data(max_entries=12)
def _contar(caminho_str: str, mtime: float) -> int:
    with open(caminho_str, "rb") as arq:
        return max(0, sum(1 for _ in arq) - 1)  # menos o cabecalho


def contar_linhas(caminho) -> int | None:
    """Registros da tabela, sem o cabecalho; `None` se o arquivo nao existe.

    Contar bytes e barato onde carregar o DataFrame nao e: e assim que a tela
    diz "107.123 ligantes" sem ler 107.123 linhas."""
    caminho = Path(caminho)
    if not caminho.is_file():
        return None
    return _contar(str(caminho), caminho.stat().st_mtime)


@st.cache_data(max_entries=12)
def _ler_json(caminho_str: str, mtime: float):
    try:
        return json.loads(Path(caminho_str).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def ler_manifesto(caminho):
    """O manifesto de uma etapa; `None` se falta ou esta quebrado.

    Manifesto ilegivel recebe o mesmo tratamento de manifesto ausente, e por
    um motivo: nos dois casos a tabela nao e citavel, e e isso que a tela tem
    de dizer. Derrubar a aba nao informaria nada melhor."""
    caminho = Path(caminho)
    if not caminho.is_file():
        return None
    return _ler_json(str(caminho), caminho.stat().st_mtime)


def caminhos_da_etapa(nome: str) -> tuple[Path, Path]:
    """`(tabela, manifesto)` daquela etapa, pelos nomes que o `main.py` usa.

    Os nomes vem da tabela ETAPAS e a pasta de `caminhos.etapa`: assim
    renomear uma saida e uma mudanca num lugar so."""
    etapa = _etapas().get(nome)
    if etapa is None:
        raise ValueError(
            f"etapa desconhecida: {nome!r}; as etapas sao "
            f"{', '.join(caminhos.ETAPAS)}")
    pasta = caminhos.TSV / nome
    return pasta / Path(etapa["saida"]).name, pasta / Path(etapa["manifesto"]).name


def temporario_da_etapa(nome: str) -> Path | None:
    """O `.tmp` que a etapa esta escrevendo para a TABELA daquela etapa.

    `atomico.temporario` nomeia o temporario `<destino>.<aleatorio>.tmp` no
    mesmo diretorio do destino, e e por isso que encontra-lo e um glob e nao
    um palpite. O temporario do manifesto fica de fora: ele nao tem linhas
    para mostrar, e confundi-lo com a tabela poria um parcial de uma linha na
    tela. Se houver mais de um (execucao interrompida deixa o anterior), vale
    o mais recente."""
    tabela, _ = caminhos_da_etapa(nome)
    # O `.tmp` pode virar o destino entre o glob e o stat -- e exatamente o
    # instante da publicacao, quando o pesquisador esta olhando a tela.
    datados = []
    for candidato in tabela.parent.glob(tabela.name + ".*.tmp"):
        try:
            datados.append((candidato.stat().st_mtime, candidato))
        except OSError:
            continue  # publicado no meio da varredura
    return max(datados)[1] if datados else None


def copiar_parcial(nome: str, destino, limite: int | None = LIMITE_PADRAO):
    """Copia o parcial da etapa para `destino` e devolve quantos registros a
    copia tem; `None` se nao ha execucao escrevendo aquela tabela.

    **A ultima linha de dados e descartada sempre**, completa ou nao. A etapa
    esta escrevendo: a linha pode estar cortada no meio de um campo (e
    `ler_tsv` levanta quando o numero de campos nao bate), ou o `\\n` pode ter
    sido gravado com o resto ainda no buffer. Julgar qual dos dois e adivinhar;
    uma linha a menos num resultado parcial nao custa nada.

    O cabecalho e preservado: ele nao e dado, e sem ele a copia nao e um TSV
    do projeto."""
    tmp = temporario_da_etapa(nome)
    if tmp is None:
        return None
    teto = None if limite is None else limite + 2  # cabecalho + limite + a descartada
    linhas: list[bytes] = []
    truncada = False
    try:
        with open(tmp, "rb") as arq:
            for bruta in arq:
                if not bruta.endswith(b"\n"):
                    truncada = True  # fim do arquivo no meio de uma linha
                    break
                linhas.append(bruta)
                if teto is not None and len(linhas) >= teto:
                    break
    except OSError:
        return None  # o temporario virou o destino durante a copia: publicou
    cabecalho, dados = (linhas[:1], linhas[1:]) if linhas else ([], [])
    if not truncada and dados:
        dados.pop()
    if limite is not None:
        dados = dados[:limite]
    destino = Path(destino)
    destino.parent.mkdir(parents=True, exist_ok=True)
    with open(destino, "wb") as copia:
        copia.writelines(cabecalho + dados)
    return len(dados)


def estado_da_etapa(nome: str) -> dict:
    """O que existe no disco para aquela etapa.

    `estado` e um de tres:

    * `"ausente"`  -- nada no disco para esta etapa.
    * `"parcial"`  -- ha um temporario sendo escrito. Vence a tabela
      publicada: e a execucao em curso que o pesquisador esta acompanhando, e
      a tabela publicada, se existe, e a ANTERIOR -- continua em `tabela`, e
      continua citavel pelo manifesto dela.
    * `"pronta"`   -- tabela publicada, nenhuma execucao escrevendo.

    Tabela presente com manifesto ausente ainda e `"pronta"`: o estado
    descreve a TABELA. Quem avisa que ela nao e citavel e a tela -- esconder
    a tabela que o pesquisador acabou de gerar seria pior."""
    tabela, manifesto = caminhos_da_etapa(nome)
    existe = tabela.is_file()
    parcial = temporario_da_etapa(nome) is not None
    if not existe and not parcial:
        return {"estado": "ausente", "tabela": None, "manifesto": None,
                "linhas": None, "mtime": None}
    return {"estado": "parcial" if parcial else "pronta",
            "tabela": tabela if existe else None,
            "manifesto": ler_manifesto(manifesto) if existe else None,
            "linhas": contar_linhas(tabela) if existe else None,
            "mtime": tabela.stat().st_mtime if existe else None}
