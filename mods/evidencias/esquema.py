# -*- coding: utf-8 -*-
"""
esquema.py

Esquema de `evidencias.tsv`: cada linha e UMA evidencia sobre um par
ligante/alvo, medida (ChEMBL) ou prevista (modelo).

Regras do formato
-----------------
* TSV, nunca CSV: 42,8% dos nomes de produto natural tem virgula.
* Separador decimal e o ponto em toda coluna legivel por maquina. A virgula
  so pode aparecer no texto de exibicao de `atividade_relatada`.
* Escrita em fluxo: cada linha vai direto ao disco, nada fica em memoria
  (a base real tem milhoes de linhas).
* Procedencia nunca se mistura: `pchembl` so recebe medida; probabilidade
  de modelo nao tem coluna aqui.

Escape (reversivel)
-------------------
Barra invertida como caractere de escape, como no formato TSV do MySQL:
    \\ -> \\\\     tabulacao -> \\t     quebra de linha -> \\n
    retorno de carro -> \\r
Na leitura a sequencia e revertida em uma unica passada, entao o texto
original volta exato, inclusive se ele ja continha "\\t" literal.
"""

import re
from pathlib import Path

COLUNAS_EVIDENCIA = [
    "ligante", "ligante_inchikey_parental", "ligante_chembl_id",
    "alvo_chembl_id", "alvo", "organismo", "tipo_alvo", "origem_alvo",
    "tipo_evidencia", "fonte", "atividade_relatada", "medida", "relacao",
    "valor", "unidade", "pchembl", "tipo_ensaio", "ensaio_descricao",
    "triagem_hts", "confianca_ensaio", "n_moleculas_treino", "documento",
    "pmid", "ano", "procedencia",
]

_ESCAPE = {"\\": "\\\\", "\t": "\\t", "\n": "\\n", "\r": "\\r"}
_RE_ESCAPAR = re.compile(r"[\\\t\n\r]")
_DESESCAPE = {"\\": "\\", "t": "\t", "n": "\n", "r": "\r"}
_RE_DESESCAPAR = re.compile(r"\\(.)", re.DOTALL)


def escapar(valor) -> str:
    """Normaliza um campo para TSV; None vira vazio. Reversivel."""
    if valor is None:
        return ""
    return _RE_ESCAPAR.sub(lambda m: _ESCAPE[m.group()], str(valor))


def _desescapar(texto: str) -> str:
    # Sequencia desconhecida (arquivo escrito por outra ferramenta) fica como esta.
    return _RE_DESESCAPAR.sub(
        lambda m: _DESESCAPE.get(m.group(1), m.group()), texto)


def formatar_atividade(medida: str, relacao: str, valor, unidade: str) -> str:
    """Texto de exibicao, ex. `IC50 = 120 nM`. O chamador deve passar `valor`
    como numero (nao texto ja formatado); assim o decimal sai com ponto."""
    if isinstance(valor, float) and valor.is_integer():
        valor = int(valor)
    partes = [medida, relacao, valor, unidade]
    return " ".join(str(p) for p in partes if p not in (None, ""))


class EscritorEvidencias:
    """Gerenciador de contexto que grava linhas uma a uma, sem acumular."""

    def __init__(self, caminho):
        self.caminho = Path(caminho)
        self._arq = None

    def __enter__(self):
        self._arq = open(self.caminho, "w", encoding="utf-8", newline="\n")
        self._arq.write("#" + "\t".join(COLUNAS_EVIDENCIA) + "\n")
        return self

    def __exit__(self, *exc):
        self._arq.close()
        self._arq = None
        return False

    def escrever(self, linha: dict) -> None:
        desconhecidas = set(linha) - set(COLUNAS_EVIDENCIA)
        if desconhecidas:
            raise ValueError(f"colunas fora do esquema: {sorted(desconhecidas)}")
        self._arq.write(
            "\t".join(escapar(linha.get(c)) for c in COLUNAS_EVIDENCIA) + "\n")


def ler_evidencias(caminho):
    """Gera um dict por linha de dados, lendo em fluxo."""
    with open(caminho, encoding="utf-8", newline="\n") as arq:
        colunas = None
        for numero, bruta in enumerate(arq, start=1):
            texto = bruta.rstrip("\n")
            if colunas is None:
                if not texto.startswith("#"):
                    raise ValueError("cabecalho '#...' ausente em " + str(caminho))
                colunas = texto[1:].split("\t")
                continue
            if not texto:
                continue
            campos = texto.split("\t")
            if len(campos) != len(colunas):
                raise ValueError(
                    f"linha {numero} de {caminho}: {len(campos)} campos "
                    f"encontrados, {len(colunas)} esperados")
            yield {c: _desescapar(v) for c, v in zip(colunas, campos)}
