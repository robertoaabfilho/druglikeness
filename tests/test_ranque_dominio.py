# -*- coding: utf-8 -*-
"""Guardas da zona de confiabilidade da previsao.

O que estes testes fixam e o resultado do adendo 8: dentro da faixa prevista,
a previsao na zona medida com habilidade ordena acima da previsao na zona
medida no azar, por maior que seja a probabilidade desta.

Cada guarda e provada contra o caso que ela TEM de pegar.
"""
import pytest

from mods.evidencias.esquema import TIPO_MEDIDO, TIPO_PREDITO
from mods.ranque import dominio
from mods.ranque.faixa import chave, classificar


def _prev(valor, origem="usuario", alvo="Receptor opioide mu"):
    return {"tipo_evidencia": TIPO_PREDITO, "origem_alvo": origem,
            "valor": str(valor), "alvo_pref_name": alvo, "alvo": alvo,
            "alvo_chembl_id": "CHEMBL233"}


def _med(pchembl, origem="usuario", alvo="Receptor opioide mu"):
    return {"tipo_evidencia": TIPO_MEDIDO, "origem_alvo": origem,
            "pchembl": str(pchembl), "alvo_pref_name": alvo, "alvo": alvo,
            "alvo_chembl_id": "CHEMBL233"}


# ---------------------------------------------------------------- zona
@pytest.mark.parametrize("sim, esperado", [
    (0.30, dominio.ANALOGO),
    (0.25, dominio.ANALOGO),      # limite fechado embaixo
    (0.2499, dominio.BANDA),
    (0.20, dominio.BANDA),        # limite fechado embaixo
    (0.1999, dominio.FORA),
    (0.0, dominio.FORA),
])
def test_zona_nos_limites_medidos(sim, esperado):
    assert dominio.zona(sim, n_referencias=5000) == esperado


def test_poucas_referencias_viram_indeterminado_e_nao_fora():
    """O caso que a guarda tem de pegar: alvo com bibliografia estreita.

    Medido: CHEMBL3896 tem 19 referencias e 99,9% da biblioteca abaixo de
    0,20. Chamar isso de `fora` rebaixaria o alvo por falta de bibliografia.
    """
    assert dominio.zona(0.05, n_referencias=19) == dominio.INDETERMINADO
    assert dominio.zona(0.05, n_referencias=dominio.MIN_REFERENCIAS - 1) == \
        dominio.INDETERMINADO
    assert dominio.zona(0.05, n_referencias=dominio.MIN_REFERENCIAS) == \
        dominio.FORA


def test_sem_similaridade_e_indeterminado_nunca_fora():
    assert dominio.zona(None, n_referencias=5000) == dominio.INDETERMINADO


def test_ordem_poe_o_indeterminado_entre_medido_e_azar():
    assert dominio.ordem(dominio.ANALOGO) == dominio.ordem(dominio.BANDA)
    assert dominio.ordem(dominio.BANDA) < dominio.ordem(dominio.INDETERMINADO)
    assert dominio.ordem(dominio.INDETERMINADO) < dominio.ordem(dominio.FORA)


def test_zona_ausente_nao_e_punida_como_azar():
    assert dominio.ordem(None) == dominio.ordem(dominio.INDETERMINADO)
    assert dominio.ordem("") == dominio.ordem(dominio.INDETERMINADO)


# ---------------------------------------------------------------- sim_top5
def test_sim_top5_e_media_das_cinco_maiores():
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
    fp = gen.GetFingerprint(Chem.MolFromSmiles("CCO"))
    iguais = [gen.GetFingerprint(Chem.MolFromSmiles("CCO")) for _ in range(5)]
    assert dominio.sim_top5(fp, iguais) == pytest.approx(1.0)
    # com menos de k referencias, a media e sobre o que existe
    assert dominio.sim_top5(fp, iguais[:2]) == pytest.approx(1.0)
    assert dominio.sim_top5(fp, []) is None
    assert dominio.sim_top5(None, iguais) is None


def test_sim_top5_cai_com_referencia_distante():
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
    fp = gen.GetFingerprint(Chem.MolFromSmiles("CCO"))
    longe = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in
             ("c1ccc2c(c1)nc(n2CCN3CCOCC3)c4ccc(cc4)S(=O)(=O)N",)]
    assert dominio.sim_top5(fp, longe) < 0.2


# ------------------------------------------------- efeito no ranqueamento
def test_banda_vence_probabilidade_maior_fora_do_dominio():
    """O ganho de performance que o adendo 8 mediu, em forma de teste.

    p = 0,99 na zona onde o modelo ficou em AUROC 0,4987 nao deve ordenar
    acima de p = 0,55 na zona onde ele deu 0,6819.
    """
    zonas_fora = {"CHEMBL233": dominio.FORA}
    zonas_banda = {"CHEMBL233": dominio.BANDA}
    forte_sem_dominio = chave(
        classificar([_prev(0.99)], p_negativo=None,
                    zona_por_alvo=zonas_fora), "AAA")
    fraca_na_banda = chave(
        classificar([_prev(0.55)], p_negativo=None,
                    zona_por_alvo=zonas_banda), "ZZZ")
    assert fraca_na_banda < forte_sem_dominio


def test_dentro_da_mesma_zona_a_probabilidade_manda():
    z = {"CHEMBL233": dominio.BANDA}
    forte = chave(classificar([_prev(0.9)], p_negativo=None,
                              zona_por_alvo=z), "AAA")
    fraca = chave(classificar([_prev(0.5)], p_negativo=None,
                              zona_por_alvo=z), "AAA")
    assert forte < fraca


def test_evidencia_medida_nao_depende_do_dominio_do_modelo():
    """`pchembl` e medido; a zona e propriedade do modelo. Nao se cruzam."""
    fora = chave(classificar([_med(8.0)], p_negativo=None,
                             zona_por_alvo={"CHEMBL233":
                                            dominio.FORA}), "AAA")
    banda = chave(classificar([_med(8.0)], p_negativo=None,
                              zona_por_alvo={"CHEMBL233":
                                             dominio.BANDA}), "AAA")
    assert fora == banda
    assert classificar([_med(8.0)], p_negativo=None)["zona_dominio"] == ""


def test_sem_medida_de_dominio_o_ranking_nao_muda():
    """Ausencia de medida nao vira punicao: so a ordem relativa importa."""
    linhas = [[_prev(v)] for v in (0.9, 0.5, 0.1)]
    com = [chave(classificar(l, p_negativo=None), "K") for l in linhas]
    assert com == sorted(com)


def test_faixa_prevista_registra_a_zona():
    r = classificar([_prev(0.7)], p_negativo=None,
                    zona_por_alvo={"CHEMBL233": dominio.ANALOGO})
    assert r["faixa"] == "2a" and r["zona_dominio"] == dominio.ANALOGO


def test_parametros_vao_para_o_manifesto():
    p = dominio.parametros()
    assert p["limite_banda"] == 0.20 and p["limite_analogo"] == 0.25
    assert p["minimo_referencias"] == dominio.MIN_REFERENCIAS
