# -*- coding: utf-8 -*-
"""
esquema.py

Esquema das quatro saidas da etapa 2. Mesmas regras de formato da etapa 1
(`mods/evidencias/esquema.py`), reaproveitadas, nao reescritas: TSV, cabecalho
com um '#', escape reversivel, ponto decimal.

Por que tabelas diferentes
--------------------------
* `propriedades`: uma linha por LIGANTE. ADMET e descritor nao sao afirmacao
  sobre um par ligante/alvo, e enfia-los no `evidencias.tsv` exigiria um alvo
  vazio em cada linha -- justamente o que torna `pchembl` a unica coluna
  comparavel.
* `portao`: uma linha por LIGANTE, so o veredito. Separada de `propriedades`
  porque e a tabela que a fase 3 le, e ela tem de caber na tela.
* `alertas_bula`: uma linha por ACHADO (ligante, secao, termo). O rastro e o
  produto: sem o trecho e a secao, "tem alerta" nao se confere.
* `alvos_anotados`: uma linha por ALVO. Expressao e propriedade do alvo.

Colunas de ADMET nao estao fixadas aqui: sao as 104 que o ADMET-AI publica,
com os nomes dele, descobertas em tempo de execucao. Fixar a lista aqui faria
o repositorio mentir sobre a versao instalada.
"""

from pathlib import Path

from mods.evidencias.esquema import desescapar, escapar

COLUNAS_PROPRIEDADES = [
    "ligante", "ligante_inchikey_parental", "ligante_chembl_id", "smiles",
    "peso_molecular", "doadores_hb", "aceptores_hb", "tpsa",
    "ligacoes_rotacionaveis", "refratividade_molar", "atomos_pesados", "qed",
    "alerta_pains", "alerta_brenk", "alerta_nih", "procedencia",
]

COLUNAS_ALERTAS_BULA = [
    "ligante", "ligante_inchikey_parental", "termo", "secao", "forca",
    "trecho", "n_bulas_consultadas", "campo_do_nome", "procedencia",
]

COLUNAS_ALVOS_EXTRA = [
    "expressao_maxima", "tecido_maximo", "expressao_nos_tecidos_pedidos",
    "tecidos_pedidos_casados", "unidade_expressao", "fonte_expressao",
    "escore_rna_expression", "expressao_veredito", "expressao_observacao",
]


def colunas_portao(criterios) -> list[str]:
    """Cabecalho de `portao.tsv`: identidade, veredito e, por criterio, o
    veredito dele e o valor que o decidiu. A lista de criterios vem da
    politica, nao de uma constante: promover um criterio nao deve exigir
    editar o esquema."""
    colunas = ["ligante", "ligante_inchikey_parental", "aprovado", "motivos"]
    for c in criterios:
        colunas += [f"{c}_veredito", f"{c}_valor"]
    return colunas + ["politica"]


class EscritorTSV:
    """Gravacao em fluxo, uma linha por vez, com cabecalho fixo.

    A base real tem 108.935 ligantes; `propriedades.tsv` com 120 colunas nao
    cabe em memoria como lista de dicts, e nao precisa caber."""

    def __init__(self, caminho, colunas: list[str]):
        self.caminho = Path(caminho)
        self.colunas = list(colunas)
        self._arq = None
        self.linhas = 0

    def __enter__(self):
        self._arq = open(self.caminho, "w", encoding="utf-8", newline="\n")
        self._arq.write("#" + "\t".join(self.colunas) + "\n")
        return self

    def __exit__(self, *exc):
        self._arq.close()
        self._arq = None
        return False

    def escrever(self, linha: dict) -> None:
        desconhecidas = set(linha) - set(self.colunas)
        if desconhecidas:
            raise ValueError(f"colunas fora do esquema: {sorted(desconhecidas)}")
        self._arq.write(
            "\t".join(formatar(linha.get(c)) for c in self.colunas) + "\n")
        self.linhas += 1


def formatar(valor) -> str:
    """Celula pronta para o TSV. Booleano vira 'sim'/'nao' (a convencao do
    `alvos_enriquecidos.tsv` da etapa 1) e float sai com ponto decimal e sem
    notacao cientifica curta demais para reconstruir o numero."""
    if isinstance(valor, bool):
        return "sim" if valor else "nao"
    if isinstance(valor, float):
        if valor != valor or valor in (float("inf"), float("-inf")):
            return ""  # NaN/inf do ADMET-AI: celula vazia diz 'sem valor'
        return escapar(repr(round(valor, 6)) if abs(valor) < 1e15 else repr(valor))
    return escapar(valor)


def ler_tsv(caminho):
    """Gera um dict por linha de dados. Mesmo contrato de `ler_evidencias`."""
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
            yield {c: desescapar(v) for c, v in zip(colunas, campos)}
