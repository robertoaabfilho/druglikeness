# -*- coding: utf-8 -*-
import inspect

from mods.evidencias.predito import linhas_preditas

LIGANTES = [
    {"nome": "aspirina", "smiles": "CC(=O)Oc1ccccc1C(=O)O",
     "inchikey_parental": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"},
    {"nome": "cafeina", "smiles": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
     "inchikey_parental": "RYYVLZVUVIJVGH-UHFFFAOYSA-N"},
]
ALVOS = [{"alvo_chembl_id": "CHEMBL100"}, {"alvo_chembl_id": "CHEMBL101"},
         {"alvo_chembl_id": "CHEMBL103"}]


def test_respeita_o_corte(modelo_falso, rotulos_falsos):
    linhas = list(linhas_preditas(LIGANTES, ALVOS, checkpoint=modelo_falso,
                                  rotulos=rotulos_falsos, corte=0.99))
    assert all(float(l["pchembl"] or 0) == 0 for l in linhas)
    assert all(l["tipo_evidencia"] == "predito" for l in linhas)


def test_corte_filtra_de_verdade(modelo_falso, rotulos_falsos):
    todas = list(linhas_preditas(LIGANTES, ALVOS, checkpoint=modelo_falso,
                                 rotulos=rotulos_falsos, corte=0.0))
    assert len(todas) == len(LIGANTES) * len(ALVOS)
    # pesos aleatorios dao p perto de 0.5; 0.505 separa os pares sem depender
    # de valor exato
    corte = 0.505
    acima = list(linhas_preditas(LIGANTES, ALVOS, checkpoint=modelo_falso,
                                 rotulos=rotulos_falsos, corte=corte))
    assert 0 < len(acima) < len(todas)
    assert all(float(l["atividade_relatada"].split()[2]) >= 0.5 for l in acima)


def test_alvo_fora_do_modelo_nao_gera_linha(modelo_falso, rotulos_falsos):
    alvos = ALVOS + [{"alvo_chembl_id": "CHEMBL_INEXISTENTE", "no_modelo": False}]
    linhas = list(linhas_preditas(LIGANTES, alvos, checkpoint=modelo_falso,
                                  rotulos=rotulos_falsos, corte=0.0))
    assert not any(l["alvo_chembl_id"] == "CHEMBL_INEXISTENTE" for l in linhas)


def test_nunca_preenche_pchembl(modelo_falso, rotulos_falsos):
    linhas = list(linhas_preditas(LIGANTES, ALVOS, checkpoint=modelo_falso,
                                  rotulos=rotulos_falsos, corte=0.0))
    assert all(l["pchembl"] == "" for l in linhas)
    assert all(l["n_moleculas_treino"] != "" for l in linhas)


def test_formato_das_linhas_e_geracao_preguicosa(modelo_falso, rotulos_falsos):
    gen = linhas_preditas(LIGANTES, ALVOS, checkpoint=modelo_falso,
                          rotulos=rotulos_falsos, corte=0.0)
    assert inspect.isgenerator(gen)
    l = next(gen)
    assert l["atividade_relatada"].startswith("p = 0.")
    assert l["atividade_relatada"].endswith(" (ECFP4)")
    assert l["procedencia"] == "predito (PP4CKD ECFP4 ecfp4_falso.pt)"
    assert l["n_moleculas_treino"] == 4000
    assert l["ligante"] == "aspirina"
