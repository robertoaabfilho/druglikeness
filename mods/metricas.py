# -*- coding: utf-8 -*-
"""
metricas.py - a AUC que decide o desenho do ranqueamento.

AUC pela estatistica de Mann-Whitney: a fracao de pares (positivo, negativo)
em que o positivo vem antes, **com empate valendo meio**. E a medida que o
projeto usa para perguntar "a ordem poe quem TRATA antes de quem CAUSA?".

Este calculo vivia dentro de um script de validacao que importava sua
configuracao de forma plana e, por isso, nao era importavel sob o pytest --
o que levou quatro helpers a reimplementa-lo dentro de um arquivo de teste.
Aqui ele e um modulo, com uma unica definicao.

**`nan`, nao zero, quando falta um dos lados.** Zero afirmaria "ordem
invertida"; ausencia de dado nao afirma nada. E a mesma regra que faz
`indeterminado` nunca penalizar um composto no ranqueamento.

A granularidade importa na leitura: com 16 positivos e 3 negativos sao 48
pares, logo o menor passo e 1/48 ~ 0,021. Um valor de AUC sobre tres negativos
nao sustenta conclusao -- so direcao.
"""
import numpy as np
from scipy.stats import rankdata


def calcular_auc(positivos: list[float], negativos: list[float]) -> float:
    """Fracao de pares em que o positivo supera o negativo, empate = meio.

    `positivos` e a classe que se quer no topo. Devolve `nan` quando um dos
    lados esta vazio, e arredonda a 4 casas, que e como o projeto cita.
    """
    n_pos, n_neg = len(positivos), len(negativos)
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    valores = np.asarray(list(positivos) + list(negativos), dtype=float)
    ranks = rankdata(valores)
    soma_ranks_pos = ranks[:n_pos].sum()
    return round(
        (soma_ranks_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg), 4)
