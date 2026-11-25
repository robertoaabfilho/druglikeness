# -*- coding: utf-8 -*-
import pytest

from mods.evidencias.chembl import (AlvoAmbiguo, ListaTruncada, alvos_por_uniprot,
                                    classe_de_proteina, resolver_alvo_por_nome,
                                    resolver_ligante, uniprot_do_alvo)


def test_uniprot_mapeia_para_o_single_protein(respostas):
    achados = alvos_por_uniprot("P23975", "Homo sapiens", obter_json=respostas)
    assert achados[0]["target_chembl_id"] == "CHEMBL222"
    assert achados[0]["target_type"] == "SINGLE PROTEIN"


def test_preferencia_por_single_protein_sem_filtro_do_servidor():
    # Sem single protein no servidor, cai no consulta sem filtro e ordena:
    # o grupo vem primeiro na API, a proteina deve sair primeiro.
    def falso(url):
        if "target_type" in url:
            return {"targets": []}
        return {"targets": [{"target_chembl_id": "G", "target_type": "SELECTIVITY GROUP"},
                            {"target_chembl_id": "P", "target_type": "SINGLE PROTEIN"}]}
    assert [a["target_chembl_id"] for a in alvos_por_uniprot("X", "", obter_json=falso)] == ["P", "G"]


def test_uniprot_consulta_single_protein_no_servidor_primeiro():
    pedidos = []

    def falso(url):
        pedidos.append(url)
        return {"targets": [{"target_chembl_id": "P", "target_type": "SINGLE PROTEIN"}]}
    alvos_por_uniprot("X", "", obter_json=falso)
    assert len(pedidos) == 1 and "target_type=SINGLE%20PROTEIN" in pedidos[0]


def test_alvo_ambiguo_levanta_com_os_candidatos(respostas):
    # Nunca escolher sozinho entre ortologo, complexo e grupo.
    with pytest.raises(AlvoAmbiguo) as erro:
        resolver_alvo_por_nome("noradrenaline transporter", "", obter_json=respostas)
    assert len(erro.value.candidatos) > 1
    for c in erro.value.candidatos:
        assert c["target_chembl_id"] and c["organism"]


def test_ligante_exato_tem_prioridade(respostas):
    assert resolver_ligante("ZEUITGRIYCTCEM-KRWDZBQOSA-N", "",
                            obter_json=respostas) == "CHEMBL1175"


def test_bloco_com_estereo_igual_resolve(respostas):
    # Exato falha (ultima camada diferente); o bloco acha o mesmo estereo.
    assert resolver_ligante("ZEUITGRIYCTCEM-KRWDZBQOSA-M", "",
                            obter_json=respostas) == "CHEMBL1175"


def test_bloco_com_estereo_diferente_devolve_none(respostas):
    # Existem varios estereoisomeros no bloco e nenhum e este: outro
    # enantiomero seria outra molecula.
    assert resolver_ligante("", "ZEUITGRIYCTCEM-AAAAAAAAAA-N",
                            obter_json=respostas) is None


def _molecula(ident, chave):
    return {"molecule_chembl_id": ident,
            "molecule_structures": {"standard_inchi_key": chave}}


def test_candidato_unico_de_estereo_diferente_nao_resolve():
    def falso(url):
        if "istartswith" in url:
            return {"molecules": [_molecula("CHEMBL9", "AAAAAAAAAAAAAA-BBBBBBBBBB-N")]}
        return {"molecules": []}
    assert resolver_ligante("AAAAAAAAAAAAAA-CCCCCCCCCC-N", "", obter_json=falso) is None


def test_bloco_com_estereo_plano_resolve_se_inequivoco():
    def falso(url):
        if "istartswith" in url:
            return {"molecules": [_molecula("CHEMBL9", "AAAAAAAAAAAAAA-UHFFFAOYSA-N")]}
        return {"molecules": []}
    # Pedido com estereo, candidato plano
    assert resolver_ligante("AAAAAAAAAAAAAA-CCCCCCCCCC-N", "", obter_json=falso) == "CHEMBL9"
    # Pedido plano, candidato com estereo
    def falso2(url):
        if "istartswith" in url:
            return {"molecules": [_molecula("CHEMBL8", "AAAAAAAAAAAAAA-BBBBBBBBBB-N")]}
        return {"molecules": []}
    assert resolver_ligante("AAAAAAAAAAAAAA-UHFFFAOYSA-N", "", obter_json=falso2) == "CHEMBL8"


def test_bloco_plano_com_varios_candidatos_nao_chuta():
    def falso(url):
        if "istartswith" in url:
            return {"molecules": [_molecula("A", "AAAAAAAAAAAAAA-BBBBBBBBBB-N"),
                                  _molecula("B", "AAAAAAAAAAAAAA-CCCCCCCCCC-N")]}
        return {"molecules": []}
    assert resolver_ligante("AAAAAAAAAAAAAA-UHFFFAOYSA-N", "", obter_json=falso) is None


def test_lista_truncada_levanta():
    def falso(url):
        return {"page_meta": {"total_count": 99}, "molecules": [_molecula("A", "X-Y-N")]}
    with pytest.raises(ListaTruncada):
        resolver_ligante("AAAAAAAAAAAAAA-CCCCCCCCCC-N", "", obter_json=falso)

    def falso_alvos(url):
        return {"page_meta": {"total_count": 99}, "targets": [{"target_chembl_id": "A"}]}
    with pytest.raises(ListaTruncada):
        resolver_alvo_por_nome("x", "", obter_json=falso_alvos)


def test_caminho_inverso_alvo_para_uniprot(respostas):
    assert uniprot_do_alvo("CHEMBL222", obter_json=respostas) == "P23975"


def test_classe_de_proteina_para_o_colapso(respostas):
    assert classe_de_proteina("CHEMBL222", obter_json=respostas)
