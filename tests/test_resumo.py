# -*- coding: utf-8 -*-
"""Resumos que o painel mostra no upload; testados sem subir o Streamlit."""
from mods.evidencias.entrada import Problema
from mods.evidencias.resumo import resumo_alvos, resumo_ligantes


def test_resumo_de_alvos_conta_os_problemas(rotulos_falsos):
    # O roteiro original usava CHEMBL222, que os rotulos falsos (CHEMBL100..105)
    # nao tem; aqui CHEMBL100 e o alvo presente.
    alvos = [{"alvo_chembl_id": "CHEMBL100", "organismo": "Homo sapiens"},
             {"alvo": "algum alvo", "organismo": ""},
             {"alvo_chembl_id": "CHEMBL_X", "organismo": "Homo sapiens"}]
    r = resumo_alvos(alvos, rotulos_falsos)
    assert r["sem_organismo"] == 1
    assert r["fora_do_modelo"] == 1
    assert r["total"] == 3 and r["no_modelo"] == 1
    assert r["pouco_treino"] == 0  # CHEMBL100 tem 4000 ativos


def test_resumo_de_alvos_aponta_pouco_treino_e_nome_com_organismo(rotulos_falsos):
    alvos = [{"alvo_chembl_id": "CHEMBL101"},  # 9 ativos
             {"alvo": "ALVO FALSO 0", "organismo": "homo sapiens"}]  # nome casa
    r = resumo_alvos(alvos, rotulos_falsos)
    assert r["pouco_treino"] == 1 and r["fora_do_modelo"] == 0
    assert [a["alvo_chembl_id"] for a in r["lista_pouco_treino"]] == ["CHEMBL101"]


def test_resumo_de_ligantes_separa_validos_de_invalidos():
    r = resumo_ligantes([{"smiles": "CCO"}, {"smiles": "xyz!"}], [])
    assert r["validos"] == 1 and r["invalidos"] == 1


def test_resumo_de_ligantes_soma_linhas_ja_descartadas_pela_leitura():
    problemas = [Problema(3, "smiles", "x: smiles que o RDKit nao le: ?", "erro"),
                 Problema(0, "nome", "duplicado", "aviso")]
    r = resumo_ligantes([{"smiles": "CCO"}], problemas)
    assert r["validos"] == 1 and r["invalidos"] == 1
    assert r["erros"] == 1 and r["avisos"] == 1


def test_resumo_de_alvos_conta_o_que_a_leitura_ja_descartou(rotulos_falsos):
    problemas = [Problema(2, "organismo", "x: organismo obrigatorio", "erro")]
    r = resumo_alvos([{"alvo_chembl_id": "CHEMBL100"}], rotulos_falsos,
                     problemas=problemas)
    assert r["sem_organismo"] == 1 and r["total"] == 2 and r["no_modelo"] == 1


def test_comando_de_coleta_nao_enriquece_e_o_de_enriquecer_nao_preve():
    from mods.evidencias.resumo import montar_comando
    base = dict(python="py", ligantes="l", alvos="a", sintomas="s", saida="o",
                rotulos="r")
    c = montar_comando("coletar", checkpoint="m.pt", **base)
    assert "--checkpoint" in c and "--alvos-enriquecidos" not in c
    e = montar_comando("enriquecer", alvos_enriquecidos="x.tsv",
                       efos={"dor": "EFO_1"}, checkpoint="m.pt", **base)
    assert "--checkpoint" not in e and "dor=EFO_1" in e
    import pytest
    with pytest.raises(ValueError):
        montar_comando("enriquecer", **base)


def test_chave_de_enriquecimento_muda_com_efo_e_parametros():
    from mods.evidencias.resumo import chave_enriquecimento as k
    base = k("abc", {"dor": "EFO_1"}, 0.4, 150, None)
    assert base == k("abc", {"dor": "EFO_1"}, 0.4, 150, None)
    assert base != k("abc", {"dor": "EFO_2"}, 0.4, 150, None)
    assert base != k("abd", {"dor": "EFO_1"}, 0.4, 150, None)
    assert base != k("abc", {"dor": "EFO_1"}, 0.5, 150, None)
    assert base != k("abc", {"dor": "EFO_1"}, 0.4, 150, 10)


def test_edicao_pendente_so_com_edicao_de_verdade():
    from mods.evidencias.resumo import edicao_pendente
    assert not edicao_pendente(None)
    assert not edicao_pendente({"edited_rows": {}, "added_rows": [], "deleted_rows": []})
    assert edicao_pendente({"edited_rows": {1: {"alvo": "x"}}, "added_rows": [],
                            "deleted_rows": []})


def test_tabela_de_alvos_vazia_nao_quebra(tmp_path):
    from mods.evidencias.resumo import ler_tabela_de_alvos
    vazio = tmp_path / "v.tsv"
    vazio.write_text("")
    colunas, linhas = ler_tabela_de_alvos(vazio)
    assert "alvo_chembl_id" in colunas and linhas == []
    so_cabecalho = tmp_path / "c.tsv"
    so_cabecalho.write_text("#alvo_chembl_id\talvo\n")
    assert ler_tabela_de_alvos(so_cabecalho) == (["alvo_chembl_id", "alvo"], [])
