# -*- coding: utf-8 -*-
"""
vistas.py - recortes das tabelas para a tela. Funcoes puras.

Nada aqui decide: as faixas foram atribuidas pela fase 3, o veredito pelo
portao, a expressao pela etapa 2. Estas funcoes agrupam, filtram e nomeiam --
e o lugar onde errar seria mentir sobre o que a tabela diz.

Duas regras do pipeline que esta camada tem de honrar, e que os testes prendem:

**`indeterminado` nunca vira aprovacao.** Criterio nao pedido sai
`indeterminado`, nunca `passa`. Traduzir isso para "sem alerta" na tela
desfaria a distincao que o portao existe para manter: ausencia de documento
nao e ausencia de risco.

**A ordem das faixas e a da evidencia, nao a da contagem.** 1a antes de 1b
antes de 2a: medida antes de prevista, alvo do usuario antes de alvo
enriquecido. Ordenar pela contagem poria a faixa mais populosa no topo e
perderia a hierarquia que ordena o ranking inteiro.
"""
import pandas as pd

from mods.ranque.faixa import FAIXAS

# O que cada faixa significa, na linguagem de quem vai citar o resultado.
# Vem do desenho da fase 3, nao de rotulo inventado na tela.
LEGENDAS = {
    "1a": "evidencia MEDIDA (ChEMBL) sobre um alvo que VOCE declarou; "
          "ordena por pchembl",
    "1b": "evidencia MEDIDA sobre um alvo que o enriquecimento acrescentou "
          "(Open Targets ou Reactome); ordena por pchembl",
    "2a": "evidencia PREVISTA pelo modelo sobre um alvo que VOCE declarou; "
          "ordena por probabilidade",
    "2b": "evidencia PREVISTA sobre um alvo que o enriquecimento acrescentou; "
          "ordena por probabilidade",
    "3": "sem evidencia em alvo afirmativo: NAO RANQUEIA, e aparece na tabela "
         "de nao ranqueados com o motivo",
}

# Rotulo de tela para cada veredito do portao. `indeterminado` fica como esta:
# nao ha palavra mais curta que nao minta.
SITUACOES = {
    "sim": "aprovado",
    "nao": "reprovado",
    "indeterminado": "indeterminado",
    "passa": "passa",
    "marca": "marca",
    "reprova": "reprova",
}

ORDEM_DA_FORCA = {"forte": 0, "moderada": 1, "fraca": 2}


def legenda_da_faixa(faixa: str) -> str:
    """O que a faixa significa; string vazia para faixa que nao existe."""
    return LEGENDAS.get(str(faixa), "")


def resumo_das_faixas(df: pd.DataFrame) -> pd.DataFrame:
    """Contagem por faixa, na ordem da evidencia, com a legenda ao lado.

    Faixa sem nenhuma linha nao aparece: uma linha com zero diria que a faixa
    foi avaliada e ficou vazia, e nem sempre foi."""
    if df is None or "faixa" not in getattr(df, "columns", []) or not len(df):
        return pd.DataFrame(columns=["faixa", "ligantes", "significado"])
    contagem = df["faixa"].value_counts()
    linhas = [{"faixa": f, "ligantes": int(contagem[f]),
               "significado": legenda_da_faixa(f)}
              for f in FAIXAS if f in contagem.index]
    return pd.DataFrame(linhas, columns=["faixa", "ligantes", "significado"])


def filtrar_faixa(df: pd.DataFrame, faixas) -> pd.DataFrame:
    """Nenhuma faixa marcada significa TUDO, nao nada: uma tela que comeca
    vazia parece resultado vazio."""
    if df is None or not len(df) or not faixas:
        return df
    return df[df["faixa"].isin(list(faixas))]


def filtrar_alvo(df: pd.DataFrame, alvos) -> pd.DataFrame:
    if df is None or not len(df) or not alvos or "alvo" not in df.columns:
        return df
    return df[df["alvo"].isin(list(alvos))]


def portao_para_tela(df: pd.DataFrame) -> pd.DataFrame:
    """Acrescenta `situacao`, legivel, SEM tocar em `aprovado`.

    A coluna original fica: ela e a que se cita, e quem exportar a tela tem de
    exportar o token do pipeline, nao a traducao."""
    if df is None or not len(df):
        return df
    copia = df.copy()
    coluna = "aprovado" if "aprovado" in copia.columns else "portao_aprovado"
    copia["situacao"] = copia[coluna].map(
        lambda v: SITUACOES.get(str(v).strip().lower(), str(v)))
    return copia


def motivos_dos_nao_ranqueados(df: pd.DataFrame) -> pd.DataFrame:
    """Quantos ligantes por motivo de nao ranquear.

    Sem coluna `posicao`: quem nao ranqueia nao tem posicao, e uma coluna
    vazia convidaria a ordenar esta tabela como se fosse ranking."""
    if df is None or "motivo" not in getattr(df, "columns", []) or not len(df):
        return pd.DataFrame(columns=["motivo", "ligantes"])
    contagem = df["motivo"].value_counts()
    return pd.DataFrame({"motivo": contagem.index,
                         "ligantes": contagem.to_numpy(dtype=int)})


def alertas_por_forca(df: pd.DataFrame) -> pd.DataFrame:
    """Os achados de bula, mais fortes primeiro, com o trecho literal intacto.

    O rastro e o produto: sem o trecho e a secao, "tem alerta" nao se confere,
    e a negacao em texto livre ("no evidence of neuropathy") continua nao
    resolvida -- e por isso que o trecho fica exposto para leitura humana, em
    vez de virar um booleano."""
    if df is None or not len(df):
        return df
    copia = df.copy()
    if "forca" in copia.columns:
        copia = copia.assign(_ordem=copia["forca"].map(
            lambda v: ORDEM_DA_FORCA.get(str(v).strip().lower(), 9)))
        copia = copia.sort_values("_ordem", kind="stable").drop(columns="_ordem")
    return copia


def alvos_com_expressao(df: pd.DataFrame) -> pd.DataFrame:
    """A lista de alvos anotados, com as colunas de expressao a frente.

    Expressao e propriedade do ALVO: um alvo de expressao baixa nao reprova
    ligante nenhum, porque o ligante pode ter outros alvos. A tela mostra o
    veredito do alvo, e nao o soma ao veredito do ligante."""
    if df is None or not len(df):
        return df
    frente = [c for c in ("alvo", "alvo_chembl_id", "origem_alvo",
                          "expressao_veredito", "expressao_maxima",
                          "tecido_maximo", "no_modelo")
              if c in df.columns]
    return df[frente + [c for c in df.columns if c not in frente]]
