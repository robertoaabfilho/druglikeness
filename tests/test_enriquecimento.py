# -*- coding: utf-8 -*-
"""Open Targets e Reactome, so com respostas gravadas (nenhum teste usa a rede).

As respostas do Reactome em tests/dados/respostas/reactome_data_participants_*
foram aparadas aos campos que o codigo le (dbId, displayName, stId,
databaseName, identifier, name, moleculeType, schemaClass, className); o
original traz milhares de identificadores cruzados por proteina."""
import pytest

from mods.evidencias.open_targets import (
    alvos_associados, alvos_do_chembl, buscar_efo, colapsar_familia,
    escrever_alvos_enriquecidos)
from mods.evidencias.reactome import (
    alvos_expandidos_no_modelo, expandir_vias)

SEMENTE = [{"uniprot": "P35372", "alvo_chembl_id": "CHEMBL233"}]


def test_efo_devolve_candidatos_sem_escolher(respostas):
    cands = buscar_efo("neuropathic pain", obter_json=respostas)
    assert cands[0]["id"] == "EFO_0005762"
    assert len(cands) > 1


def test_colapso_de_familia_guarda_os_descartados():
    # B esta fora do modelo: dobra em A, que fica com o registro.
    alvos = [{"alvo_chembl_id": "A", "classe": "VGCC", "score": 0.605, "no_modelo": True},
             {"alvo_chembl_id": "B", "classe": "VGCC", "score": 0.603, "no_modelo": False},
             {"alvo_chembl_id": "C", "classe": "SLC6", "score": 0.583, "no_modelo": True}]
    saida = colapsar_familia(alvos)
    assert [a["alvo_chembl_id"] for a in saida] == ["A", "C"]
    assert saida[0]["familia_colapsados"] == "B"
    assert saida[0]["familia"] == "VGCC"


def test_colapso_nunca_funde_membros_presentes_no_modelo():
    # NET e SERT: mesma classe, ambos no modelo (duloxetina age nos dois).
    alvos = [{"alvo_chembl_id": "CHEMBL228", "classe": "SLC06", "score": 0.58, "no_modelo": True},
             {"alvo_chembl_id": "CHEMBL222", "classe": "SLC06", "score": 0.57, "no_modelo": True}]
    saida = colapsar_familia(alvos)
    assert [a["alvo_chembl_id"] for a in saida] == ["CHEMBL228", "CHEMBL222"]
    assert not any(a.get("familia_colapsados") for a in saida)


def test_colapso_preserva_usuario_e_alvo_sem_classe():
    alvos = [{"alvo_chembl_id": "A", "classe": "X", "score": 0.9, "no_modelo": True},
             {"alvo_chembl_id": "U", "classe": "X", "score": 0.1, "no_modelo": False,
              "origem": "usuario"},
             {"alvo_chembl_id": "S1", "classe": None, "no_modelo": False},
             {"alvo_chembl_id": "S2", "classe": None, "no_modelo": False}]
    assert [a["alvo_chembl_id"] for a in colapsar_familia(alvos)] == ["A", "U", "S1", "S2"]


def test_colapso_sem_nenhum_membro_no_modelo_deixa_um_representante():
    alvos = [{"alvo_chembl_id": "A", "classe": "X", "score": 0.5, "no_modelo": False},
             {"alvo_chembl_id": "B", "classe": "X", "score": 0.6, "no_modelo": False}]
    saida = colapsar_familia(alvos)
    assert [a["alvo_chembl_id"] for a in saida] == ["B"]
    assert saida[0]["familia_colapsados"] == "A"


def test_alvos_associados_traz_escore_tipos_e_uniprot(respostas):
    alvos = alvos_associados("EFO_0005762", corte=0.5, obter_json=respostas)
    assert len(alvos) == 30
    assert all(a["score"] >= 0.5 for a in alvos)
    assert [a["score"] for a in alvos] == sorted((a["score"] for a in alvos), reverse=True)
    por_gene = {a["gene"]: a for a in alvos}
    assert por_gene["SLC6A2"]["uniprot"] == "P23975"
    assert por_gene["SLC6A4"]["uniprot"] == "P31645"
    assert set(por_gene["SLC6A2"]["datatypes"]) <= {"literature", "clinical", "genetic_association"}


def test_alvos_associados_levanta_com_resposta_de_erro():
    from mods.evidencias.chembl import RespostaInesperada
    with pytest.raises(RespostaInesperada):
        alvos_associados("EFO_X", corte=0.1,
                         obter_json=lambda url, corpo=None: {"errors": [{"message": "x"}]})


def test_via_grande_demais_e_descartada(respostas):
    # 'G alpha (i) signalling events' tem 317 proteinas; 'Opioid Signalling', 90.
    rel = {}
    expandidos = expandir_vias(SEMENTE, teto=150, obter_json=respostas, relatorio=rel)
    vias = {a["via_origem"] for a in expandidos}
    assert "R-HSA-111885" in vias
    assert "R-HSA-418594" not in vias
    assert rel["vias_descartadas"]["R-HSA-418594"] == 317


def test_expansao_registra_a_semente(respostas):
    expandidos = expandir_vias(SEMENTE, teto=150, obter_json=respostas)
    assert all(a["alvo_semente"] == "CHEMBL233" for a in expandidos)
    assert all(a["origem"] == "reactome" for a in expandidos)
    assert "P35372" not in {a["uniprot"] for a in expandidos}


def test_alvo_do_usuario_nunca_desaparece_do_enriquecido(tmp_path):
    # Restricao global: o enriquecimento so acrescenta.
    destino = tmp_path / "alvos_enriquecidos.tsv"
    escrever_alvos_enriquecidos(destino, [
        {"alvo_chembl_id": "CHEMBL222", "origem": "usuario"},
        {"alvo_chembl_id": "CHEMBL233", "origem": "open_targets"},
    ])
    texto = destino.read_text(encoding="utf-8")
    assert "CHEMBL222" in texto and "usuario" in texto


def test_escrita_usa_ponto_decimal_e_booleano_legivel(tmp_path):
    destino = tmp_path / "e.tsv"
    escrever_alvos_enriquecidos(destino, [
        {"alvo_chembl_id": "CHEMBL222", "origem": "open_targets", "score": 0.5,
         "no_modelo": True}])
    cabecalho, linha = destino.read_text(encoding="utf-8").splitlines()
    cols = cabecalho.lstrip("#").split("\t")
    celulas = dict(zip(cols, linha.split("\t")))
    assert celulas["escore_associacao"] == "0.5" and celulas["no_modelo"] == "sim"


def _chembl_falso(mapa):
    def obter(url):
        for acc, alvo in mapa.items():
            if f"accession={acc}" in url:
                return {"targets": [alvo], "page_meta": {"total_count": 1}}
        return {"targets": [], "page_meta": {"total_count": 0}}
    return obter


def test_filtro_do_modelo_mede_o_corte():
    alvo = lambda i: {"target_chembl_id": i, "pref_name": i, "organism": "Homo sapiens",
                      "target_type": "SINGLE PROTEIN"}
    rotulos = [{"chembl_id": "CHEMBL1", "n_ativos": 10}, {"chembl_id": "CHEMBL3", "n_ativos": 5}]
    mapa = {"P1": alvo("CHEMBL1"), "P2": alvo("CHEMBL2"), "P3": alvo("CHEMBL3")}
    exp = [{"uniprot": u, "origem": "reactome", "via_origem": "R", "alvo_semente": "S"}
           for u in ("P1", "P2", "P3", "P4")]
    rel = {}
    saida = alvos_expandidos_no_modelo(exp, rotulos, obter_json_chembl=_chembl_falso(mapa),
                                       ja_presentes={"CHEMBL3"}, relatorio=rel)
    assert [a["alvo_chembl_id"] for a in saida] == ["CHEMBL1"]
    assert saida[0]["n_moleculas_treino"] == 10
    assert (rel["proteinas"], rel["no_modelo"], rel["fora_do_modelo"],
            rel["sem_alvo_chembl"], rel["ja_na_lista"]) == (4, 1, 1, 1, 1)


def test_alvos_do_chembl_marca_modelo_e_registra_falha_por_item():
    from mods.evidencias.chembl import ListaTruncada
    alvo = {"target_chembl_id": "CHEMBL222", "pref_name": "NET",
            "organism": "Homo sapiens", "target_type": "SINGLE PROTEIN"}
    def obter(url):
        if "accession=BAD" in url:
            return {"targets": [alvo] * 2, "page_meta": {"total_count": 9}}
        return _chembl_falso({"P23975": alvo})(url)
    cands = [{"gene": "SLC6A2", "uniprot": "P23975", "score": 0.58},
             {"gene": "XXX", "uniprot": "BAD", "score": 0.5},
             {"gene": "YYY", "uniprot": "", "score": 0.4}]
    saida = alvos_do_chembl(cands, [{"chembl_id": "CHEMBL222", "n_ativos": 3493}],
                            obter_json_chembl=obter)
    assert [a["no_modelo"] for a in saida] == [True, False, False]
    assert saida[0]["n_moleculas_treino"] == 3493
    assert "ChEMBL" in saida[1]["observacao"] and "Swiss-Prot" in saida[2]["observacao"]


def test_http_500_e_falha_de_item_mas_timeout_propaga():
    class Http500(OSError):
        response = object()
    def obter_500(url):
        raise Http500("500 Server Error")
    exp = [{"uniprot": "P1", "origem": "reactome", "via_origem": "R", "alvo_semente": "S"}]
    rel = {}
    assert alvos_expandidos_no_modelo(exp, [], obter_json_chembl=obter_500,
                                      relatorio=rel) == []
    assert "500" in rel["problemas"][0]
    def sem_rede(url):
        raise TimeoutError("timeout")
    with pytest.raises(TimeoutError):
        alvos_expandidos_no_modelo(exp, [], obter_json_chembl=sem_rede)
