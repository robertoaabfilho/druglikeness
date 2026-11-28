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


# ===========================================================================
# Espectro de alvos, decisao 3: o teto do modelo, medido
# ===========================================================================
# "O manifesto da etapa 1 passa a registrar quantos alvos estao fora do espaco
# do modelo e quantos ficaram sem evidencia por isso. E o numero que justifica
# o projeto 2; sem ele, retreinar e opiniao."
#
# A separacao que decide: alvo fora do modelo COM evidencia medida nao custou
# nada ao pipeline (ele ranqueia em 1b, acima de qualquer previsao). Alvo fora
# do modelo SEM evidencia medida e exatamente o que retreinar recuperaria.
from mods.evidencias.resumo import teto_do_modelo


def _alvo(ident, no_modelo="sim"):
    return {"alvo_chembl_id": ident, "alvo": ident, "organismo": "Homo sapiens",
            "no_modelo": no_modelo, "origem": "reactome"}


def _medida(alvo_id):
    return {"alvo_chembl_id": alvo_id, "tipo_evidencia": "medido",
            "pchembl": "8.1"}


def test_teto_conta_alvo_dentro_e_fora_do_modelo():
    t = teto_do_modelo([_alvo("CHEMBL1"), _alvo("CHEMBL9", "nao"),
                        _alvo("CHEMBL8", "nao")], [], n_ligantes=10)
    assert t["alvos_total"] == 3
    assert t["alvos_no_modelo"] == 1
    assert t["alvos_fora_do_modelo"] == 2


def test_teto_separa_fora_do_modelo_com_e_sem_evidencia_medida():
    """O numero que justifica o projeto 2."""
    t = teto_do_modelo([_alvo("CHEMBL9", "nao"), _alvo("CHEMBL8", "nao")],
                       [_medida("CHEMBL9")], n_ligantes=10)
    assert t["alvos_fora_com_evidencia_medida"] == 1
    assert t["alvos_fora_sem_evidencia_medida"] == 1


def test_teto_conta_os_pares_previstos_perdidos():
    """n_ligantes x alvos fora do modelo: os pares que o preditor nao pontua.
    E a grandeza comparavel entre execucoes de tamanhos diferentes."""
    t = teto_do_modelo([_alvo("CHEMBL1"), _alvo("CHEMBL9", "nao")], [],
                       n_ligantes=1000)
    assert t["pares_previstos_possiveis"] == 1000
    assert t["pares_previstos_perdidos"] == 1000


def test_teto_de_lista_toda_no_modelo_nao_perde_par_nenhum():
    t = teto_do_modelo([_alvo("CHEMBL1"), _alvo("CHEMBL2")], [], n_ligantes=5)
    assert t["alvos_fora_do_modelo"] == 0
    assert t["pares_previstos_perdidos"] == 0
    assert t["alvos_fora_sem_evidencia_medida"] == 0


def test_teto_de_lista_vazia_nao_levanta():
    t = teto_do_modelo([], [], n_ligantes=0)
    assert t["alvos_total"] == 0 and t["pares_previstos_perdidos"] == 0


def test_teto_so_conta_evidencia_MEDIDA_como_recuperada():
    """Previsao num alvo fora do modelo nao existe -- mas se existisse na
    tabela por engano, contar como "nao perdeu" esconderia o teto."""
    t = teto_do_modelo([_alvo("CHEMBL9", "nao")],
                       [{"alvo_chembl_id": "CHEMBL9",
                         "tipo_evidencia": "predito", "valor": "0.8"}],
                       n_ligantes=10)
    assert t["alvos_fora_sem_evidencia_medida"] == 1


def test_teto_aceita_a_marca_como_booleano():
    """A lista vem do TSV (`sim`/`nao`) ou de memoria (True/False)."""
    t = teto_do_modelo([{"alvo_chembl_id": "CHEMBL9", "no_modelo": False}], [],
                       n_ligantes=1)
    assert t["alvos_fora_do_modelo"] == 1


def test_teto_conta_os_dobrados_por_familia():
    """Review Focus 2: `colapsar_familia` dobra alvo fora do modelo no
    representante da classe, usando `no_modelo is False` como criterio -- a
    premissa que a decisao 2 derruba. O colapso roda so sobre a lista do Open
    Targets, entao a decisao 2 nao e anulada; mas se os dobrados nao forem
    contados, o numero que justifica o projeto 2 os esconde."""
    alvos = [{"alvo_chembl_id": "CHEMBL1", "no_modelo": "sim",
              "familia_colapsados": "CHEMBL7;CHEMBL8"}]
    t = teto_do_modelo(alvos, [], n_ligantes=1)
    assert t["dobrados_por_familia"] == 2


def test_marca_ausente_nao_e_tomada_como_dentro_do_modelo():
    """ACHADO na execucao de verdade: `input/alvos.txt` e uma lista de IDs e
    NAO tem a coluna `no_modelo`. Tomar a marca ausente como "esta no modelo"
    fazia o teto sair zero sempre -- o numero que justifica o projeto 2
    mentindo por omissao, que e o pior modo de ele estar errado."""
    alvos = [{"alvo_chembl_id": "CHEMBL9"}]  # sem a coluna, e fora dos rotulos
    t = teto_do_modelo(alvos, [], n_ligantes=10,
                       rotulos=[{"chembl_id": "CHEMBL1", "nome": "A",
                                 "organismo": "Homo sapiens", "n_ativos": 10}])
    assert t["alvos_fora_do_modelo"] == 1
    assert t["alvos_no_modelo"] == 0


def test_marca_ausente_resolvida_pelos_rotulos_quando_o_alvo_existe():
    t = teto_do_modelo([{"alvo_chembl_id": "CHEMBL1"}], [], n_ligantes=10,
                       rotulos=[{"chembl_id": "CHEMBL1", "nome": "A",
                                 "organismo": "Homo sapiens", "n_ativos": 10}])
    assert t["alvos_no_modelo"] == 1 and t["alvos_fora_do_modelo"] == 0


def test_sem_rotulos_a_marca_ausente_fica_indeterminada():
    """Sem o modelo em maos (execucao sem checkpoint), nao se afirma nem o
    contrario: o alvo entra em `alvos_marca_ausente`, e o teto nao finge
    saber."""
    t = teto_do_modelo([{"alvo_chembl_id": "CHEMBL9"}], [], n_ligantes=10,
                       rotulos=None)
    assert t["alvos_marca_ausente"] == 1
    assert t["alvos_no_modelo"] == 0 and t["alvos_fora_do_modelo"] == 0


def test_os_rotulos_vencem_a_marca_e_a_divergencia_e_contada():
    """REVERTE a decisao anterior ("a marca explicita vence os rotulos").

    A pergunta que o teto responde e "o que ESTE modelo nao cobre", e os
    rotulos em maos SAO este modelo; a marca registra o que algum modelo
    anterior cobria. Com a marca vencendo, o manifesto afirmava
    `pares_previstos_perdidos` sobre um alvo que o preditor pontua -- e, no
    outro sentido, `pares_previstos_possiveis` que a execucao nao entrega.
    A divergencia continua visivel, como contagem."""
    rotulos = [{"chembl_id": "CHEMBL1", "nome": "A",
                "organismo": "Homo sapiens", "n_ativos": 10}]
    t = teto_do_modelo([{"alvo_chembl_id": "CHEMBL1", "no_modelo": "nao"}], [],
                       n_ligantes=1, rotulos=rotulos)
    assert t["alvos_no_modelo"] == 1
    assert t["alvos_fora_do_modelo"] == 0
    assert t["alvos_marca_divergente"] == 1


def test_os_pares_possiveis_nao_se_chamam_feitos():
    """Sem checkpoint nenhum par e pontuado, e o numero continua positivo: ele
    mede o que o espaco do modelo PERMITE, nao o que a execucao fez."""
    t = teto_do_modelo([{"alvo_chembl_id": "CHEMBL1", "no_modelo": "sim"}], [],
                       n_ligantes=50)
    assert t["pares_previstos_possiveis"] == 50
    assert "pares_previstos_feitos" not in t


# ===========================================================================
# Correcoes da revisao: o contador do teto estava errado por construcao
# ===========================================================================
def test_alvo_declarado_por_nome_e_organismo_conta_como_dentro(rotulos_falsos):
    """C1: `input/alvos.txt` aceita `alvo + organismo` sem ID -- esta no README
    e em `entrada.py`. O teto resolvia a marca ausente comparando STRING de ID,
    e um alvo sem ID caia em "fora do modelo" mesmo com o modelo cobrindo-o:
    o numero que argumenta A FAVOR de retreinar, errado, em silencio, num modo
    de entrada documentado."""
    alvos = [{"alvo_chembl_id": "", "alvo": "ALVO FALSO 0",
              "organismo": "homo sapiens"}]
    t = teto_do_modelo(alvos, [], n_ligantes=10, rotulos=rotulos_falsos)
    assert t["alvos_no_modelo"] == 1
    assert t["alvos_fora_do_modelo"] == 0


def test_nome_ambiguo_nao_e_fora_do_modelo(rotulos_falsos):
    """Nome que casa com mais de um rotulo e INDECIDIVEL sem o ID -- nao e
    prova de que o modelo nao cobre o alvo. `resumo_alvos` ja separa os dois
    casos; o teto os colapsava num so e inflava a perda."""
    alvos = [{"alvo_chembl_id": "", "alvo": "repetido", "organismo": "x"}]
    rotulos = [{"chembl_id": "CHEMBL1", "nome": "repetido", "organismo": "x",
                "n_ativos": 10},
               {"chembl_id": "CHEMBL2", "nome": "repetido", "organismo": "x",
                "n_ativos": 10}]
    t = teto_do_modelo(alvos, [], n_ligantes=10, rotulos=rotulos)
    assert t["alvos_ambiguos"] == 1
    assert t["alvos_fora_do_modelo"] == 0
    assert t["pares_previstos_perdidos"] == 0


def test_medida_sem_pchembl_nao_conta_como_ranqueada(rotulos_falsos):
    """I5: o README diz que esta coluna sao os alvos que "nao custaram nada"
    porque ranquearam por medida. Medida sem `pchembl` finito NAO ranqueia --
    `faixa.classificar` a manda para a faixa 3 com `medida_sem_pchembl`. Contar
    esse alvo como recuperado esconde uma perda real."""
    alvos = [{"alvo_chembl_id": "CHEMBL9", "no_modelo": "nao"}]
    linhas = [{"alvo_chembl_id": "CHEMBL9", "tipo_evidencia": "medido",
               "pchembl": "", "valor": "31622.8"}]
    t = teto_do_modelo(alvos, linhas, n_ligantes=1, rotulos=rotulos_falsos)
    assert t["alvos_fora_com_evidencia_medida"] == 0
    assert t["alvos_fora_sem_evidencia_medida"] == 1


def test_medida_com_pchembl_conta_como_ranqueada(rotulos_falsos):
    alvos = [{"alvo_chembl_id": "CHEMBL9", "no_modelo": "nao"}]
    linhas = [{"alvo_chembl_id": "CHEMBL9", "tipo_evidencia": "medido",
               "pchembl": "7.9"}]
    t = teto_do_modelo(alvos, linhas, n_ligantes=1, rotulos=rotulos_falsos)
    assert t["alvos_fora_com_evidencia_medida"] == 1


def test_os_rotulos_vencem_a_marca_quando_estao_em_maos(rotulos_falsos):
    """I6: a pergunta que o teto responde e "o que ESTE modelo nao cobre", e os
    rotulos em maos SAO este modelo. A marca registra o que algum modelo
    anterior cobria. Com a marca vencendo, o manifesto afirmava
    `pares_previstos_possiveis` que a execucao provadamente nao entrega."""
    alvos = [{"alvo_chembl_id": "CHEMBL100", "no_modelo": "nao"}]
    t = teto_do_modelo(alvos, [], n_ligantes=10, rotulos=rotulos_falsos)
    assert t["alvos_no_modelo"] == 1
    assert t["alvos_marca_divergente"] == 1


def test_divergencia_no_outro_sentido_tambem_e_contada(rotulos_falsos):
    alvos = [{"alvo_chembl_id": "CHEMBL_INEXISTENTE", "no_modelo": "sim"}]
    t = teto_do_modelo(alvos, [], n_ligantes=10, rotulos=rotulos_falsos)
    assert t["alvos_fora_do_modelo"] == 1
    assert t["alvos_marca_divergente"] == 1


def test_sem_rotulos_a_marca_decide(rotulos_falsos):
    """Sem o modelo em maos a marca e tudo o que ha, e ai ela decide."""
    t = teto_do_modelo([{"alvo_chembl_id": "CHEMBL9", "no_modelo": "nao"}], [],
                       n_ligantes=10, rotulos=None)
    assert t["alvos_fora_do_modelo"] == 1
    assert t["alvos_marca_divergente"] == 0


def test_os_quatro_estados_somam_o_total(rotulos_falsos):
    alvos = [{"alvo_chembl_id": "CHEMBL100"},
             {"alvo_chembl_id": "CHEMBL_X"},
             {"alvo_chembl_id": "", "alvo": "?", "organismo": ""}]
    t = teto_do_modelo(alvos, [], n_ligantes=1, rotulos=rotulos_falsos)
    assert (t["alvos_no_modelo"] + t["alvos_fora_do_modelo"]
            + t["alvos_ambiguos"] + t["alvos_marca_ausente"]) == t["alvos_total"]
