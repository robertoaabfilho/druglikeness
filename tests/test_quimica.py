# -*- coding: utf-8 -*-
"""A camada quimica: molecula, forma parental, InChIKey e normalizacao de campo.

A forma parental e a chave de identidade do projeto: e ela que reconhece que
sal e base livre sao o mesmo composto, e e por ela que a deduplicacao e o
cruzamento com o prospector acontecem.
"""
import pytest

from mods import quimica

# Cafeina: SMILES curto que o RDKit le em qualquer versao.
CAFEINA = "CN1C=NC2=C1C(=O)N(C)C(=O)N2C"
# Sal de sodio do diclofenaco: tem contra-ion, logo a forma parental tem de
# descartar o fragmento inorganico e devolver so o acido.
DICLOFENACO_SAL = "[Na+].[O-]C(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl"
DICLOFENACO_ACIDO = "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl"


def test_molecula_le_smiles_valido():
    assert quimica.molecula(CAFEINA) is not None


def test_molecula_devolve_none_em_vez_de_levantar():
    """SMILES ilegivel nao derruba a leitura: ele vira None e a linha e
    JULGADA depois, nunca descartada em silencio na leitura."""
    assert quimica.molecula("nao sou smiles") is None
    assert quimica.molecula("") is None


def test_forma_parental_reconhece_sal_e_base_livre_como_o_mesmo_composto():
    """A chave de cruzamento entre as duas bases. Sem isso, o sal de sodio e o
    acido livre contariam como dois compostos diferentes."""
    chave_sal = quimica.inchikey_de(
        quimica.forma_parental(quimica.molecula(DICLOFENACO_SAL)))
    chave_acido = quimica.inchikey_de(
        quimica.forma_parental(quimica.molecula(DICLOFENACO_ACIDO)))
    assert chave_sal
    assert chave_sal == chave_acido


def test_inchikey_da_molecula_como_esta_difere_do_parental_quando_ha_sal():
    """As duas colunas existem por um motivo: `inchikey` e o que a molecula
    declara, `inchikey_parental` e a identidade. No sal elas divergem."""
    mol = quimica.molecula(DICLOFENACO_SAL)
    assert quimica.inchikey_de(mol) != quimica.inchikey_de(
        quimica.forma_parental(mol))


def test_forma_parental_de_none_e_none():
    assert quimica.forma_parental(None) is None


def test_inchikey_de_none_e_vazio_nao_excecao():
    assert quimica.inchikey_de(None) == ""


@pytest.mark.parametrize("entrada", [None, "", "  ", "nan", "NaN", "none",
                                     "-", "null"])
def test_texto_trata_ausencia_como_vazio(entrada):
    """O CSV do prospector escreve ausencia de varias formas; todas sao a
    mesma coisa, e nenhuma pode virar o texto literal 'nan' num rotulo."""
    assert quimica.texto(entrada) == ""


def test_texto_preserva_o_valor_e_tira_o_espaco():
    assert quimica.texto("  1,8-cineol  ") == "1,8-cineol"
