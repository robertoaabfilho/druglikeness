# -*- coding: utf-8 -*-
"""Conjunto-controle versionado (`dados/evidencias_controle.tsv`).

Existe para que os criterios de aceite do desenho sejam verificaveis SEM rede
e sem o `evidencias_chembl.csv` de 60 MB. Ate aqui nenhum teste o lia: o
arquivo estava no repositorio e o criterio de aceite 2 nao era afirmado em
lugar nenhum.
"""
from pathlib import Path

from mods.evidencias.esquema import COLUNAS_EVIDENCIA, ler_evidencias

CONTROLE = Path(__file__).resolve().parent.parent / "dados" / "evidencias_controle.tsv"


def _linhas():
    return list(ler_evidencias(CONTROLE))


def test_criterio_2_amitriptilina_tem_medido_para_o_transportador_de_noradrenalina():
    """Criterio de aceite 2: o conjunto-controle produz evidencia em que a
    amitriptilina tem linhas `medido` para o transportador de noradrenalina
    (NET, CHEMBL222, humano)."""
    net = [l for l in _linhas()
           if l["ligante"].casefold() == "amitriptyline"
           and l["alvo_chembl_id"] == "CHEMBL222"
           and l["tipo_evidencia"] == "medido"]
    assert net, "nenhuma linha medida de amitriptilina sobre CHEMBL222"
    assert all(l["fonte"] == "ChEMBL" for l in net)
    assert all(l["organismo"] == "Homo sapiens" for l in net)
    # Rastro: sem documento nem ensaio a linha nao serve para publicacao.
    assert all(l["documento"] and l["ensaio_descricao"] for l in net)
    assert any(l["pchembl"] for l in net)


def test_criterio_3_nenhuma_linha_mistura_medido_com_predito():
    for l in _linhas():
        if l["tipo_evidencia"] == "medido":
            assert l["medida"] != "probabilidade" and l["fonte"] == "ChEMBL"
        else:
            assert l["pchembl"] == ""


def test_o_arquivo_de_controle_segue_o_esquema_e_traz_a_procedencia_do_alvo():
    linhas = _linhas()
    assert linhas
    for l in linhas:
        assert list(l) == COLUNAS_EVIDENCIA
        assert l["origem_alvo"], "origem_alvo vazio: a procedencia do alvo se perdeu"
        assert l["procedencia"]
