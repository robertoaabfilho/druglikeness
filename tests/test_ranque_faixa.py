# -*- coding: utf-8 -*-
"""Faixa de evidencia e chave de ordem da fase 3.

Desenho: docs/superpowers/specs/2026-11-26-fase3-ranqueamento-design.md
"""

import pytest

from mods.ranque.faixa import FAIXAS, chave, classificar


def _linha(tipo, origem, **extra):
    """Uma linha de `evidencias.tsv`, com os campos que a faixa le."""
    linha = {"#ligante": "x", "ligante_inchikey_parental": "K" * 27,
             "alvo_chembl_id": "CHEMBL228", "alvo": "serotonin transporter",
             "organismo": "Homo sapiens", "origem_alvo": origem,
             "tipo_evidencia": tipo, "valor": "", "unidade": "",
             "pchembl": "", "n_moleculas_treino": "900", "no_modelo": "sim"}
    linha.update(extra)
    return linha


def _med(origem="usuario", pchembl="8.0", **extra):
    return _linha("medido", origem, pchembl=pchembl, unidade="nM", **extra)


def _prev(origem="usuario", valor="0.5", **extra):
    return _linha("previsto", origem, valor=valor, **extra)


def _uma_de_cada():
    """Listas de evidencia de 1a, 1b, 2a, 2b e 3, nessa ordem."""
    return [[_med(origem="usuario", pchembl="5.0")],
            [_med(origem="open_targets", pchembl="9.9")],
            [_prev(origem="usuario", valor="0.10")],
            [_prev(origem="reactome", valor="0.99")],
            []]


def test_medido_em_alvo_do_usuario_e_faixa_1a():
    r = classificar([_med(origem="usuario", pchembl="8.2")], p_negativo=0.9)
    assert r["faixa"] == "1a" and r["valor"] == 8.2


def test_medido_em_alvo_do_enriquecimento_e_faixa_1b():
    for origem in ("open_targets", "reactome"):
        assert classificar([_med(origem=origem, pchembl="7.0")],
                           p_negativo=None)["faixa"] == "1b"


def test_canal_negativo_nao_entra_no_valor_da_faixa_1():
    # pchembl e medido, p_negativo e previsto: multiplicar misturaria os dois
    a = classificar([_med(origem="usuario", pchembl="8.2")], p_negativo=0.0)
    b = classificar([_med(origem="usuario", pchembl="8.2")], p_negativo=0.99)
    assert a["valor"] == b["valor"] == 8.2


def test_previsto_multiplica_pelo_canal_negativo():
    r = classificar([_prev(origem="usuario", valor="0.90")], p_negativo=0.5)
    assert r["faixa"] == "2a" and r["valor"] == pytest.approx(0.45)


def test_previsto_sem_canal_negativo_ordena_pelo_p_sozinho():
    r = classificar([_prev(origem="usuario", valor="0.90")], p_negativo=None)
    assert r["valor"] == pytest.approx(0.90)


def test_sem_evidencia_em_alvo_afirmativo_e_faixa_3_sem_valor():
    r = classificar([], p_negativo=0.1)
    assert r["faixa"] == "3" and r["valor"] is None


def test_agrega_varios_alvos_por_maximo_nao_por_soma():
    r = classificar([_med(origem="usuario", pchembl="6.0"),
                     _med(origem="usuario", pchembl="7.5")], p_negativo=None)
    assert r["valor"] == 7.5


def test_medida_sem_pchembl_nao_ordena_e_nao_rebaixa_para_faixa_3():
    # 1 das 49 linhas medidas do evidencias_controle.tsv e assim
    r = classificar([_med(origem="usuario", pchembl=""),
                     _prev(origem="usuario", valor="0.80")], p_negativo=None)
    assert r["faixa"] == "2a"
    assert r["medida_sem_pchembl"] == "sim"
    assert r["valor"] == pytest.approx(0.80)


def test_medida_sem_pchembl_sozinha_nao_e_faixa_1():
    r = classificar([_med(origem="usuario", pchembl="")], p_negativo=None)
    assert r["faixa"] == "3" and r["medida_sem_pchembl"] == "sim"


def test_outra_linha_medida_com_pchembl_mantem_na_faixa_1():
    r = classificar([_med(origem="usuario", pchembl=""),
                     _med(origem="usuario", pchembl="9.1")], p_negativo=None)
    assert r["faixa"] == "1a" and r["valor"] == 9.1


def test_alvo_afirmativo_fora_do_modelo_continua_faixa_1():
    # sem previsao nenhuma: medido basta
    r = classificar([_med(origem="usuario", pchembl="8.0", no_modelo="nao")],
                    p_negativo=None)
    assert r["faixa"] == "1a"


def test_chave_ordena_as_faixas_na_ordem_do_desenho():
    # Chave ASCENDENTE: `sorted` sem reverse poe 1a primeiro e a faixa 3 por
    # ultimo. A lista de `_uma_de_cada` ja esta na ordem do desenho.
    ordem = [chave(classificar(l, p_negativo=None), "K") for l in _uma_de_cada()]
    assert ordem == sorted(ordem)


def test_alvo_que_decidiu_o_valor_sai_nomeado():
    # Formato "ID nome": um escore sem o alvo responsavel nao e auditavel, e o
    # ID sozinho nao se le. Quando nao ha ID, sai so o nome.
    r = classificar([_med(origem="usuario", pchembl="6.0",
                          alvo_chembl_id="CHEMBL228", alvo="SERT"),
                     _med(origem="usuario", pchembl="7.5",
                          alvo_chembl_id="CHEMBL222", alvo="NET")],
                    p_negativo=None)
    assert r["alvo"] == "CHEMBL222 NET"


def test_alvo_sem_id_sai_so_com_o_nome():
    r = classificar([_med(alvo_chembl_id="", alvo="cst-II")], p_negativo=None)
    assert r["alvo"] == "cst-II"


def test_faixa_3_nunca_vem_antes_de_uma_faixa_com_evidencia():
    vazia = chave(classificar([], p_negativo=None), "AAA")
    fraca = chave(classificar([_prev(origem="reactome", valor="0.001")],
                              p_negativo=0.99), "ZZZ")
    assert fraca < vazia


def test_empate_de_valor_desempata_por_inchikey_crescente():
    a = chave(classificar([_med(pchembl="8.0")], p_negativo=None), "AAA")
    z = chave(classificar([_med(pchembl="8.0")], p_negativo=None), "ZZZ")
    assert a < z


def test_faixas_sao_as_cinco_do_desenho():
    assert FAIXAS == ("1a", "1b", "2a", "2b", "3")
