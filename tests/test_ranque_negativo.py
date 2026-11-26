# -*- coding: utf-8 -*-
"""O canal negativo deduzido da fase 3.

Rotulos sinteticos em `tests/dados/rotulos_minimo.tsv`: nenhum teste depende do
`target_labels.tsv` real, que nao esta no repositorio.
"""

import numpy as np

from mods.prever_alvos import carregar_rotulos
from mods.ranque.negativo import agregar, mascara

ROTULOS = carregar_rotulos("tests/dados/rotulos_minimo.tsv")

IDX_HUMANO = 0              # CHEMBL228 SERT humano, 900 ativos
IDX_NAO_HUMANO_FORTE = 1    # CHEMBL360 M. tuberculosis, 11.049 ativos
IDX_NAO_HUMANO_FRACO = 2    # CHEMBL5555 DHFR de P. falciparum, 12 ativos
IDX_NET_HUMANO = 3          # CHEMBL222 NET humano
IDX_NET_RATO = 4            # CHEMBL304 NET de rato -- ortologo do pedido


def test_alvo_humano_nunca_entra_no_canal():
    assert not mascara(ROTULOS, [])[IDX_HUMANO]


def test_alvo_pedido_nao_entra_no_canal():
    assert not mascara(ROTULOS, [IDX_NAO_HUMANO_FORTE])[IDX_NAO_HUMANO_FORTE]


def test_alvo_pouco_treinado_nao_acusa():
    # requisito 1: esparso e ausencia de evidencia, nao evidencia de atividade
    assert not mascara(ROTULOS, [])[IDX_NAO_HUMANO_FRACO]


def test_ortologo_de_alvo_pedido_fica_fora():
    # CHEMBL304 (NET de rato) e ortologo do CHEMBL222 pedido
    assert not mascara(ROTULOS, [IDX_NET_HUMANO])[IDX_NET_RATO]


def test_ortologo_entra_quando_o_alvo_humano_NAO_foi_pedido():
    # A exclusao e de ortologo de alvo PEDIDO, nao de todo homonimo
    assert mascara(ROTULOS, [])[IDX_NET_RATO]


def test_nao_humano_treinado_e_nao_pedido_entra():
    assert mascara(ROTULOS, [])[IDX_NAO_HUMANO_FORTE]


def test_minimo_treino_e_configuravel():
    assert mascara(ROTULOS, [], minimo_treino=10)[IDX_NAO_HUMANO_FRACO]


def test_agregar_e_maximo():
    assert agregar(np.array([0.1, 0.9, 0.4]),
                   np.array([True, True, False])) == 0.9


def test_agregar_ignora_o_que_esta_fora_da_mascara():
    assert agregar(np.array([0.1, 0.99]), np.array([True, False])) == 0.1


def test_mascara_vazia_devolve_none_nao_zero():
    # zero e "nenhuma atividade indesejada", que e afirmacao, nao ausencia
    assert agregar(np.array([0.9]), np.array([False])) is None
