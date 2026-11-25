# -*- coding: utf-8 -*-
import pytest

from mods.evidencias.chembl import (AlvoAmbiguo, alvos_por_uniprot,
                                    classe_de_proteina, resolver_alvo_por_nome,
                                    resolver_ligante, uniprot_do_alvo)


def test_uniprot_mapeia_para_o_single_protein(respostas):
    achados = alvos_por_uniprot("P23975", "Homo sapiens", obter_json=respostas)
    assert achados[0]["target_chembl_id"] == "CHEMBL222"
    assert achados[0]["target_type"] == "SINGLE PROTEIN"
    # O grupo de seletividade vem depois, nao some.
    assert "CHEMBL2096990" in [a["target_chembl_id"] for a in achados]


def test_alvo_ambiguo_levanta_com_os_candidatos(respostas):
    # Nunca escolher sozinho entre ortologo, complexo e grupo.
    with pytest.raises(AlvoAmbiguo) as erro:
        resolver_alvo_por_nome("noradrenaline transporter", "", obter_json=respostas)
    assert len(erro.value.candidatos) > 1


def test_ligante_sem_chembl_id_resolve_pelo_bloco_de_conectividade(respostas):
    # O ChEMBL indexa por forma salina; a parental so encontra pelo bloco.
    assert resolver_ligante("", "ZEUITGRIYCTCEM-KRWDZBQOSA-N",
                            obter_json=respostas) is not None


def test_ligante_exato_tem_prioridade(respostas):
    assert resolver_ligante("ZEUITGRIYCTCEM-KRWDZBQOSA-N", "",
                            obter_json=respostas) == "CHEMBL1175"


def test_ligante_com_estereo_desconhecido_nao_chuta(respostas):
    # O bloco existe com varios estereoisomeros e nenhum e este: None.
    assert resolver_ligante("", "ZEUITGRIYCTCEM-AAAAAAAAAA-N",
                            obter_json=respostas) is None


def test_bloco_com_um_unico_candidato_resolve():
    def falso(url):
        if "istartswith" in url:
            return {"molecules": [{"molecule_chembl_id": "CHEMBL9",
                                   "molecule_structures": {
                                       "standard_inchi_key": "AAAAAAAAAAAAAA-BBBBBBBBBB-N"}}]}
        return {"molecules": []}
    assert resolver_ligante("AAAAAAAAAAAAAA-CCCCCCCCCC-N", "", obter_json=falso) == "CHEMBL9"


def test_caminho_inverso_alvo_para_uniprot(respostas):
    assert uniprot_do_alvo("CHEMBL222", obter_json=respostas) == "P23975"


def test_classe_de_proteina_para_o_colapso(respostas):
    assert classe_de_proteina("CHEMBL222", obter_json=respostas)
