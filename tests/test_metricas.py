# -*- coding: utf-8 -*-
"""A AUC que decide o desenho do ranqueamento.

Este calculo existia duplicado em quatro helpers dentro de um teste, porque o
modulo que o continha nao era importavel sob o pytest. Agora e um modulo, e os
numeros que o projeto publica (0,750 / 0,806 / 0,583) saem daqui.

A convencao: o positivo e a classe que se quer no topo (`trata`), e empate
vale meio ponto.
"""
import math

from mods.metricas import calcular_auc


def test_separacao_perfeita_e_um():
    assert calcular_auc([2.0, 3.0], [0.5, 1.0]) == 1.0


def test_ordem_invertida_e_zero():
    assert calcular_auc([0.5], [2.0, 3.0]) == 0.0


def test_empate_vale_meio():
    """A regra que faz 12 farmacos empatados darem 0,5 em vez de 0 ou 1."""
    assert calcular_auc([1.0], [1.0]) == 0.5
    assert calcular_auc([1.0, 1.0], [1.0, 1.0]) == 0.5


def test_um_acerto_e_um_erro_em_quatro_pares_da_meio():
    assert calcular_auc([3.0, 1.0], [2.0, 0.0]) == 0.75


def test_sem_positivo_ou_sem_negativo_e_indefinido_nao_zero():
    """Zero afirmaria 'ordem invertida'. Ausencia de dado nao afirma nada --
    a mesma regra que faz `indeterminado` nunca penalizar no ranqueamento."""
    assert math.isnan(calcular_auc([1.0, 2.0], []))
    assert math.isnan(calcular_auc([], [1.0, 2.0]))
    assert math.isnan(calcular_auc([], []))


def test_arredondado_a_quatro_casas():
    """As AUC do projeto sao citadas com quatro casas (0,6146, 0,5145, 0,3438);
    o arredondamento mora aqui, nao em cada chamador."""
    valor = calcular_auc([1.0, 2.0, 3.0], [0.0, 1.5, 2.5, 4.0])
    assert valor == round(valor, 4)


def test_a_granularidade_pequena_do_conjunto_real_e_representavel():
    """16 `trata` x 3 `causa` = 48 pares, granularidade de 1/48. O calculo tem
    de devolver esse valor, nao arredondar para 0,02."""
    assert calcular_auc([1.0] * 16, [0.0, 0.0, 2.0]) == round(32 / 48, 4)
