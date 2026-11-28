# -*- coding: utf-8 -*-
"""As duas AUC do adendo, como teste.

`tests/dados/controle_previsoes.tsv` traz, por farmaco do controle, a
probabilidade em cada um dos 9 alvos pedidos, o maximo do canal negativo e a
contagem de promiscuidade -- gravados do `ecfp4_fold_4.pt` real.

O que este arquivo fixa e a ARITMETICA, nao o modelo: se o checkpoint mudar, os
numeros da fixture continuam valendo como entrada. A alternativa seria exigir o
checkpoint real na suite, e nenhum teste deste repositorio depende de arquivo
fora dele. A procedencia esta em `tests/dados/LICENSE`.
"""


import pytest

from mods.portao.esquema import ler_tsv
from mods.metricas import calcular_auc
from mods.ranque.faixa import chave, classificar

FIXTURE = "tests/dados/controle_previsoes.tsv"
PEDIDOS = ["CHEMBL1919", "CHEMBL3896", "CHEMBL228", "CHEMBL222", "CHEMBL1845",
           "CHEMBL4794", "CHEMBL233", "CHEMBL2015", "CHEMBL1904"]


def _previsoes():
    return list(ler_tsv(FIXTURE))


def _linhas_de_evidencia(previsao):
    """As probabilidades da fixture em forma de linhas de `evidencias.tsv`.

    A fase 3 classifica linhas de evidencia, e nao probabilidades cruas; montar
    as linhas aqui evita uma segunda porta de entrada em `classificar` que so o
    teste usaria."""
    return [{"ligante": previsao["nome"],
             "ligante_inchikey_parental": previsao["nome"],
             "alvo_chembl_id": alvo, "alvo": alvo,
             "organismo": "Homo sapiens", "origem_alvo": "usuario",
             "tipo_evidencia": "predito", "valor": previsao[alvo],
             "pchembl": ""}
            for alvo in PEDIDOS]


def _auc(escore_por_nome, previsoes) -> float:
    """Pares concordantes, empate valendo meio. `trata` e o positivo.

    A conta mora em `mods/metricas.py`: aqui so se separa quem trata de quem
    causa. Antes ela era reimplementada em cada arquivo de teste."""
    trata = [escore_por_nome[p["nome"]] for p in previsoes
             if p["grupo_esperado"] == "trata"]
    causa = [escore_por_nome[p["nome"]] for p in previsoes
             if p["grupo_esperado"] == "causa"]
    return calcular_auc(trata, causa)


def _escores(previsoes, *, com_canal: bool):
    escores = {}
    for previsao in previsoes:
        p_negativo = float(previsao["p_negativo"]) if com_canal else None
        resultado = classificar(_linhas_de_evidencia(previsao),
                                p_negativo=p_negativo)
        escores[previsao["nome"]] = resultado["valor"]
    return escores


def test_a_fixture_tem_os_12_farmacos_do_controle():
    previsoes = _previsoes()
    assert len(previsoes) == 12
    assert sum(1 for p in previsoes if p["grupo_esperado"] == "trata") == 6


def test_sem_canal_negativo_reproduz_a_auc_de_0_750():
    previsoes = _previsoes()
    assert _auc(_escores(previsoes, com_canal=False), previsoes) == \
        pytest.approx(0.750, abs=0.001)


def test_com_canal_negativo_reproduz_a_auc_de_0_806():
    previsoes = _previsoes()
    assert _auc(_escores(previsoes, com_canal=True), previsoes) == \
        pytest.approx(0.806, abs=0.001)


def test_o_canal_negativo_melhora_a_ordem_neste_conjunto():
    previsoes = _previsoes()
    sem = _auc(_escores(previsoes, com_canal=False), previsoes)
    com = _auc(_escores(previsoes, com_canal=True), previsoes)
    assert com > sem


def test_contar_promiscuidade_seria_pior_e_zeraria_a_amitriptilina():
    # Desenho DESCARTADO, fixado aqui para nao voltar. A contagem nao existe em
    # `mods/`: o teste a calcula da coluna `n_promiscuo` da fixture, para que
    # nenhum codigo de producao precise suportar um desenho rejeitado.
    previsoes = _previsoes()
    pior = max(int(p["n_promiscuo"]) for p in previsoes)
    afirmativo = _escores(previsoes, com_canal=False)
    escores = {p["nome"]: afirmativo[p["nome"]]
               * (1 - int(p["n_promiscuo"]) / pior) for p in previsoes}
    assert _auc(escores, previsoes) == pytest.approx(0.583, abs=0.001)
    assert escores["amitriptyline"] == 0.0


def test_somar_daria_a_mesma_auc_neste_conjunto():
    """O controle NAO discrimina `max` de soma: as duas dao 0,750.

    Medido ao escrever este teste. Importa registrar porque o desenho
    justificava `max` dizendo que soma premiaria o composto promiscuo -- com
    n = 12 isso nao se sustenta, e `max` fica como escolha DECLARADA (nao
    somar grandezas de alvos diferentes), nao como escolha medida. O que a
    medicao realmente rejeita e a CONTAGEM de promiscuidade: 0,583.
    """
    previsoes = _previsoes()
    soma = {p["nome"]: sum(float(p[a]) for a in PEDIDOS) for p in previsoes}
    assert _auc(soma, previsoes) == pytest.approx(0.750, abs=0.001)


def test_a_ordem_de_max_e_de_soma_difere_mesmo_com_a_auc_igual():
    # AUC igual nao e ordem igual: a duloxetina e a amitriptilina trocam de
    # distancia. Prova que o teste acima nao e vacuo.
    previsoes = _previsoes()
    soma = {p["nome"]: sum(float(p[a]) for a in PEDIDOS) for p in previsoes}
    maximo = _escores(previsoes, com_canal=False)
    assert soma["duloxetine"] != maximo["duloxetine"]


def test_a_ordem_do_ranking_poe_a_duloxetina_antes_do_ciprofloxacino():
    previsoes = _previsoes()
    ordem = sorted(
        previsoes,
        key=lambda p: chave(classificar(_linhas_de_evidencia(p),
                                        p_negativo=float(p["p_negativo"])),
                            p["nome"]))
    nomes = [p["nome"] for p in ordem]
    assert nomes.index("duloxetine") < nomes.index("ciprofloxacin")
