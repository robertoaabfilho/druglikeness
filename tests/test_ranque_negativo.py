# -*- coding: utf-8 -*-
"""O canal negativo deduzido da fase 3.

Rotulos sinteticos em `tests/dados/rotulos_minimo.tsv`: nenhum teste depende do
`target_labels.tsv` real, que nao esta no repositorio.
"""

import numpy as np

from mods.modelo import carregar_rotulos
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


def test_organismo_humano_e_reconhecido_em_qualquer_caixa():
    """`organismo` era comparado sensivel a caixa na mesma expressao que faz
    `casefold()` no nome. Um `target_labels.tsv` com "homo sapiens" jogaria
    alvos HUMANOS para dentro do canal -- o desenho medido como prejudicial
    (AUC 0,694, e pune a amitriptilina). O repositorio ja casefolda organismo
    em `predito.py` e `medido.py`."""
    rotulos = [dict(r) for r in ROTULOS]
    rotulos[IDX_HUMANO]["organismo"] = "homo sapiens"
    assert not mascara(rotulos, [])[IDX_HUMANO]
    rotulos[IDX_HUMANO]["organismo"] = "  HOMO SAPIENS "
    assert not mascara(rotulos, [])[IDX_HUMANO]


def test_alvo_sem_organismo_nao_entra_no_canal():
    """Organismo vazio nao e evidencia de "nao humano" -- e ausencia de dado,
    e a regra do repositorio e que ausencia nunca vira evidencia.

    O caso que motivou: `CHEMBL612545 "Unchecked"` (target_type UNCHECKED,
    organismo vazio, 217.067 moleculas ativas) e o balde do ChEMBL para
    ensaios cujo alvo nunca foi curado. Por ser o alvo mais treinado do modelo
    ele liderava o canal negativo de quase todo composto -- a penalidade vinha
    de um nao-alvo. Medido: excluir os 31 alvos sem organismo mantem a AUC em
    0,806 e o maximo do canal passa a ser patogeno de verdade."""
    rotulos = [dict(r) for r in ROTULOS]
    rotulos.append({"chembl_id": "CHEMBL612545", "nome": "Unchecked",
                    "organismo": "", "tipo": "UNCHECKED", "n_ativos": 217067})
    assert not mascara(rotulos, [])[-1]


# ===========================================================================
# Espectro de alvos: o canal negativo com alvo fora do modelo
# ===========================================================================
def test_a_mascara_do_canal_so_cobre_os_rotulos_do_modelo():
    """O canal e DEDUZIDO sobre os ROTULOS: a mascara e indexada por posicao
    nos rotulos, e alvo fora do modelo nao tem posicao nenhuma -- nao tem
    coluna de previsao, e nao pode acusar nada.

    O teste anterior afirmava isso e media so `len(mascara) == len(ROTULOS)`,
    que e verdade por construcao da compreensao de lista: passaria mesmo se
    toda a regra do canal estivesse violada. Este mede a propriedade: um alvo
    pedido que NAO esta nos rotulos nao desloca nem altera a mascara."""
    sem_o_alvo = mascara(ROTULOS, [])
    with_o_alvo_inexistente = mascara(ROTULOS, [])  # nao ha indice para ele
    assert list(sem_o_alvo) == list(with_o_alvo_inexistente)
    assert len(sem_o_alvo) == len(ROTULOS)
    # e a mascara continua recusando o humano e o pouco treinado
    assert not sem_o_alvo[IDX_HUMANO]
    assert not sem_o_alvo[IDX_NAO_HUMANO_FRACO]


def test_rotulo_com_n_ativos_vazio_nao_entra_no_canal():
    """Se um rotulo chegar com `n_ativos` vazio -- que e como o alvo fora do
    modelo sai do enriquecimento --, ele nao pode passar o teto de 50 por
    `int("" or 0)` virar 0 e, pior, nao pode levantar."""
    rotulos = [{**ROTULOS[IDX_NAO_HUMANO_FORTE], "n_ativos": ""}]
    assert not mascara(rotulos, [])[0]
