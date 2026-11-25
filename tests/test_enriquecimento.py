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
    vias = {v for a in expandidos for v in a["via_origem"].split(";")}
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


def test_via_origem_lista_todas_as_vias_que_trouxeram_a_proteina(respostas):
    # O14610 participa de 'Opioid Signalling' e de 'G-protein activation'.
    expandidos = {a["uniprot"]: a for a in expandir_vias(SEMENTE, teto=150, obter_json=respostas)}
    assert expandidos["O14610"]["via_origem"] == "R-HSA-111885;R-HSA-202040"
    # Proteina de uma via so continua com uma so.
    assert all(";" not in a["via_origem"] or "R-HSA-111885" in a["via_origem"]
               for a in expandidos.values())


def test_duas_proteinas_no_mesmo_alvo_viram_uma_linha_com_rastro_das_duas():
    alvo = {"target_chembl_id": "CHEMBL9", "pref_name": "grupo", "organism": "Homo sapiens",
            "target_type": "PROTEIN FAMILY"}
    exp = [{"uniprot": "P1", "origem": "reactome", "via_origem": "R1", "via_nome": "a", "alvo_semente": "S1"},
           {"uniprot": "P2", "origem": "reactome", "via_origem": "R2", "via_nome": "b", "alvo_semente": "S2"}]
    saida = alvos_expandidos_no_modelo(
        exp, [{"chembl_id": "CHEMBL9", "n_ativos": 1}],
        obter_json_chembl=_chembl_falso({"P1": alvo, "P2": alvo}))
    assert len(saida) == 1
    assert saida[0]["via_origem"] == "R1;R2" and saida[0]["alvo_semente"] == "S1;S2"


def test_alvo_do_reactome_tambem_associado_no_open_targets_e_anotado():
    from mods.evidencias.enriquecimento import enriquecer
    alvo = {"target_chembl_id": "CHEMBL4296", "pref_name": "SCN9A", "organism": "Homo sapiens",
            "target_type": "SINGLE PROTEIN"}
    semente = {"target_chembl_id": "CHEMBL1", "pref_name": "S", "organism": "Homo sapiens",
               "target_type": "SINGLE PROTEIN", "target_components": [{"accession": "PS"}]}
    def chembl(url):
        if "accession=Q15858" in url:
            return {"targets": [alvo], "page_meta": {"total_count": 1}}
        if "target/CHEMBL1.json" in url:
            return {**semente, "pref_name": "S", "organism": "Homo sapiens"}
        return {"targets": [], "page_meta": {"total_count": 0}}
    def ot(url, corpo=None):
        return {"data": {"disease": {"associatedTargets": {"count": 1, "rows": [
            {"score": 0.41, "datatypeScores": [],
             "target": {"id": "E", "approvedSymbol": "SCN9A",
                        "proteinIds": [{"id": "Q15858", "source": "uniprot_swissprot"}]}}]}}}}
    def reactome(url):
        if "mapping" in url:
            return [{"stId": "R-X", "displayName": "x"}]
        return [{"databaseName": "UniProt", "identifier": "Q15858"},
                {"databaseName": "UniProt", "identifier": "PS"}]
    rel = {}
    linhas = enriquecer([{"alvo_chembl_id": "CHEMBL1", "alvo": "", "organismo": ""}],
                        {"dor": "EFO_1"}, [{"chembl_id": "CHEMBL1", "n_ativos": 1},
                                           {"chembl_id": "CHEMBL4296", "n_ativos": 9}],
                        corte_associacao=0.5, teto_via=150, obter_json_chembl=chembl,
                        obter_json_open_targets=ot, obter_json_reactome=reactome,
                        relatorio=rel)
    r = [l for l in linhas if l["alvo_chembl_id"] == "CHEMBL4296"]
    assert len(r) == 1 and r[0]["origem"] == "reactome"
    assert r[0]["score"] == 0.41 and "Open Targets" in r[0]["observacao"]


def _cenario_usuario_tambem_no_open_targets():
    """Usuario fornece CHEMBL4296 (SCN9A); Open Targets e Reactome o devolvem
    de novo, e um segundo sintoma tambem o implica."""
    scn9a = {"target_chembl_id": "CHEMBL4296", "pref_name": "SCN9A",
             "organism": "Homo sapiens", "target_type": "SINGLE PROTEIN"}
    def chembl(url):
        if "target/CHEMBL4296.json" in url:
            return {**scn9a, "target_components": [{"accession": "Q15858"}]}
        if "accession=Q15858" in url:
            return {"targets": [scn9a], "page_meta": {"total_count": 1}}
        return {"targets": [], "page_meta": {"total_count": 0}}
    def ot(url, corpo=None):
        efo = corpo["variables"]["efo"]
        return {"data": {"disease": {"associatedTargets": {"count": 1, "rows": [
            {"score": 0.6 if efo == "EFO_2" else 0.5, "datatypeScores": [],
             "target": {"id": "E", "approvedSymbol": "SCN9A",
                        "proteinIds": [{"id": "Q15858", "source": "uniprot_swissprot"}]}}]}}}}
    def reactome(url):
        if "mapping" in url:
            return [{"stId": "R-X", "displayName": "x"}]
        return [{"databaseName": "UniProt", "identifier": "Q15858"}]
    return chembl, ot, reactome


def test_enriquecimento_so_acrescenta_alvo_do_usuario_fica_uma_vez_e_primeiro():
    from mods.evidencias.enriquecimento import enriquecer
    chembl, ot, reactome = _cenario_usuario_tambem_no_open_targets()
    linhas = enriquecer(
        [{"alvo_chembl_id": "CHEMBL4296", "alvo": "", "organismo": ""}],
        {"dor": "EFO_1", "nevralgia": "EFO_2"},
        [{"chembl_id": "CHEMBL4296", "n_ativos": 7}],
        corte_associacao=0.4, teto_via=150, obter_json_chembl=chembl,
        obter_json_open_targets=ot, obter_json_reactome=reactome)
    assert [l["alvo_chembl_id"] for l in linhas] == ["CHEMBL4296"]
    assert linhas[0]["origem"] == "usuario"
    assert linhas[0]["alvo"] == "SCN9A"


def test_sintomas_que_implicam_o_mesmo_alvo_se_acumulam():
    from mods.evidencias.enriquecimento import enriquecer
    chembl, ot, reactome = _cenario_usuario_tambem_no_open_targets()
    # Sem o alvo do usuario: o Open Targets o traz por dois sintomas.
    linhas = enriquecer([], {"dor": "EFO_1", "nevralgia": "EFO_2"},
                        [{"chembl_id": "CHEMBL4296", "n_ativos": 7}],
                        corte_associacao=0.4, teto_via=150, obter_json_chembl=chembl,
                        obter_json_open_targets=ot, obter_json_reactome=reactome)
    assert len(linhas) == 1 and linhas[0]["origem"] == "open_targets"
    assert linhas[0]["sintoma"] == "dor;nevralgia"
    assert linhas[0]["efo_id"] == "EFO_1;EFO_2"
    assert linhas[0]["score"] == 0.6
    # Com o alvo do usuario, ele ganha o rastro sem perder a origem.
    linhas = enriquecer([{"alvo_chembl_id": "CHEMBL4296", "alvo": "", "organismo": ""}],
                        {"dor": "EFO_1", "nevralgia": "EFO_2"},
                        [{"chembl_id": "CHEMBL4296", "n_ativos": 7}],
                        corte_associacao=0.4, teto_via=150, obter_json_chembl=chembl,
                        obter_json_open_targets=ot, obter_json_reactome=reactome)
    assert linhas[0]["origem"] == "usuario" and linhas[0]["sintoma"] == "dor;nevralgia"


def test_a_suite_bloqueia_a_rede():
    import socket
    with pytest.raises(RuntimeError, match="rede"):
        socket.create_connection(("93.184.216.34", 80), timeout=1)
