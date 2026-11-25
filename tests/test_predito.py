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
    assert linhas == []
    todas = list(linhas_preditas(LIGANTES, ALVOS, checkpoint=modelo_falso,
                                 rotulos=rotulos_falsos, corte=0.0))
    # CHEMBL100, 101 e 103 existem no modelo falso
    assert len(todas) == len(LIGANTES) * len(ALVOS)
    assert all(l["tipo_evidencia"] == "predito" for l in todas)


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
    assert all(float(l["valor"]) >= corte for l in acima)
    assert all(float(l["valor"]) < corte for l in todas if l not in acima)


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


def test_valor_e_numero_legivel_por_maquina(modelo_falso, rotulos_falsos,
                                            monkeypatch):
    import math
    import torch
    from mods.evidencias import predito

    # rede de 1 alvo que sempre devolve p = 0.9341
    rede = torch.nn.Sequential(torch.nn.Linear(4096, 1))
    with torch.no_grad():
        rede[0].weight.zero_()
        rede[0].bias.fill_(math.log(0.9341 / (1 - 0.9341)))
    monkeypatch.setattr(predito, "carregar_modelo",
                        lambda *a, **k: (rede, {"raio": 2, "bits": 4096}, "cpu", "x.pt"))
    linhas = list(linhas_preditas(LIGANTES, [{"alvo_chembl_id": "CHEMBL100"}],
                                  checkpoint="x.pt", rotulos=rotulos_falsos,
                                  corte=0.0))
    assert len(linhas) == 2
    for l in linhas:
        assert l["medida"] == "probabilidade" and l["relacao"] == "="
        assert l["unidade"] == ""
        assert abs(float(l["valor"]) - 0.9341) < 1e-5
        assert float(l["valor"]) != 0.93
        assert l["atividade_relatada"] == "p = 0.93 (ECFP4)"


def test_valor_sobrevive_a_gravacao(modelo_falso, rotulos_falsos, tmp_path):
    from mods.evidencias.esquema import EscritorEvidencias, ler_evidencias
    destino = tmp_path / "e.tsv"
    with EscritorEvidencias(destino) as w:
        for l in linhas_preditas(LIGANTES, ALVOS, checkpoint=modelo_falso,
                                 rotulos=rotulos_falsos, corte=0.0):
            w.escrever(l)
    lidas = list(ler_evidencias(destino))
    assert lidas and all("," not in l["valor"] and 0 < float(l["valor"]) < 1
                         for l in lidas)


def test_smiles_ilegivel_e_contado(modelo_falso, rotulos_falsos):
    ruins = LIGANTES + [{"nome": "lixo", "smiles": "isto nao e smiles"}]
    contagens = {}
    linhas = list(linhas_preditas(ruins, ALVOS, checkpoint=modelo_falso,
                                  rotulos=rotulos_falsos, corte=0.0,
                                  contagens=contagens))
    assert contagens["smiles_ilegiveis"] == 1
    assert len(linhas) == len(LIGANTES) * len(ALVOS)
