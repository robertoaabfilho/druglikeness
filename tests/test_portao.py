# -*- coding: utf-8 -*-
"""Testes dos modulos do portao (etapa 2), sem CLI.

O que cada teste prende e a MEDICAO que decidiu o desenho, nao a
implementacao: se alguem reescrever `decisao.py` e a amitriptilina voltar a
ser reprovada por PAINS, estes testes falham.
"""

import math

import pytest

from mods.portao import bula, decisao, esquema, expressao, propriedades

# SMILES do conjunto-controle, os mesmos do desenho da etapa 2.
AMITRIPTILINA = "CN(C)CCC=C1c2ccccc2CCc2ccccc21"
ASPIRINA = "CC(=O)Oc1ccccc1C(=O)O"


# --------------------------------------------------------------------------
# propriedades: alerta NOMEADO, e SMILES ilegivel que nao apaga o ligante
# --------------------------------------------------------------------------

def test_alerta_estrutural_e_nomeado_e_nao_apenas_contado():
    """"BRENK = 1" nao se confere; "BRENK: nitro_group" se confere. E a razao
    de este modulo existir ao lado das colunas do ADMET-AI."""
    linha, problema = propriedades.linha_de_propriedades(
        {"nome": "nitrobenzeno", "smiles": "O=[N+]([O-])c1ccccc1"})
    assert problema is None
    assert "nitro_group" in linha["alerta_brenk"]


def test_pains_pega_a_amitriptilina():
    """A medicao que decidiu que alerta estrutural MARCA: reprovar por PAINS
    apagaria o positivo de referencia do pipeline."""
    linha, _ = propriedades.linha_de_propriedades(
        {"nome": "amitriptilina", "smiles": AMITRIPTILINA})
    assert linha["alerta_pains"], "PAINS deveria pegar a amitriptilina"


def test_smiles_ilegivel_mantem_a_identidade_e_relata():
    """O ligante nao desaparece da saida so porque nao pode ser avaliado."""
    linha, problema = propriedades.linha_de_propriedades(
        {"nome": "ruim", "smiles": "C1CC", "chembl_id": "CHEMBL1"})
    assert problema == "SMILES que o RDKit nao le"
    assert linha["ligante"] == "ruim" and linha["ligante_chembl_id"] == "CHEMBL1"
    assert "qed" not in linha


def test_qed_nao_separa_tratar_de_causar():
    """QED e coluna, nunca criterio: no controle os dois QED mais altos sao do
    grupo que CAUSA. Aqui basta que a amitriptilina (trata) nao seja
    separavel da aspirina por QED alto."""
    amitriptilina, _ = propriedades.linha_de_propriedades(
        {"nome": "a", "smiles": AMITRIPTILINA})
    aspirina, _ = propriedades.linha_de_propriedades(
        {"nome": "b", "smiles": ASPIRINA})
    assert amitriptilina["qed"] > aspirina["qed"]


# --------------------------------------------------------------------------
# decisao: a politica, e o que `indeterminado` NUNCA faz
# --------------------------------------------------------------------------

def test_indeterminado_nunca_reprova_com_politica_nenhuma():
    """Ausencia de dado nao e evidencia de problema. E o caso de quase todo
    produto natural na bula e de todo alvo nao humano na expressao."""
    politica = decisao.montar_politica(promovidos=decisao.CRITERIOS)
    aprovado, motivos = decisao.aplicar(
        {c: ("indeterminado", "sem dado") for c in decisao.CRITERIOS}, politica)
    assert aprovado is True and motivos == ""


def test_marca_nao_reprova_na_politica_padrao():
    vereditos = {"alerta_estrutural": ("marca", "styrene_A(13)")}
    aprovado, motivos = decisao.aplicar(vereditos, decisao.POLITICA_PADRAO)
    assert aprovado is True and motivos == ""


def test_promover_um_criterio_passa_a_reprovar_o_mesmo_item():
    vereditos = {"alerta_estrutural": ("marca", "styrene_A(13)")}
    politica = decisao.montar_politica(promovidos=["alerta_estrutural"])
    aprovado, motivos = decisao.aplicar(vereditos, politica)
    assert aprovado is False
    assert "alerta_estrutural: styrene_A(13)" in motivos


def test_rebaixar_alerta_de_bula_para_de_reprovar():
    vereditos = {"alerta_bula": ("reprova", "neuropathy@boxed_warning")}
    assert decisao.aplicar(vereditos, decisao.POLITICA_PADRAO)[0] is False
    politica = decisao.montar_politica(rebaixados=["alerta_bula"])
    assert decisao.aplicar(vereditos, politica)[0] is True


def test_criterio_desconhecido_levanta_em_vez_de_ignorar():
    """Uma flag com erro de digitacao aceita em silencio faria o usuario
    acreditar numa politica que nao esta valendo."""
    with pytest.raises(decisao.CriterioDesconhecido) as e:
        decisao.montar_politica(promovidos=["alerta_estruturall"])
    assert "alerta_estruturall" in str(e.value)


def test_contador_diz_quanto_cada_criterio_reprovaria_se_promovido():
    """O botao e medido antes de ser apertado."""
    contador = decisao.Contador(decisao.CRITERIOS_LIGANTE)
    for veredito in ("marca", "marca", "passa"):
        vereditos = {"alerta_estrutural": (veredito, "x"),
                     "smiles": ("passa", "")}
        aprovado, _ = decisao.aplicar(vereditos, decisao.POLITICA_PADRAO)
        contador.contar(vereditos, aprovado)
    resumo = contador.resumo(decisao.POLITICA_PADRAO)
    assert resumo["itens"] == 3 and resumo["aprovados"] == 3
    assert resumo["por_criterio"]["alerta_estrutural"]["adversos"] == 2
    assert resumo["por_criterio"]["alerta_estrutural"]["politica"] == "marcar"


def test_a_linha_registra_a_politica_em_vigor():
    """Duas execucoes com politicas diferentes nao se confundem."""
    politica = decisao.montar_politica(promovidos=["toxicidade"])
    linha = decisao.linha_de_portao(
        {"ligante": "x"}, {"toxicidade": ("marca", "hERG=95")}, politica)
    assert "toxicidade=reprovar" in linha["politica"]
    assert linha["toxicidade_valor"] == "hERG=95"
    assert linha["aprovado"] is False


# --------------------------------------------------------------------------
# bula: a regra nova, e os quatro falsos positivos que a antiga produzia
# --------------------------------------------------------------------------

def _bula(**secoes):
    return {"campo_do_nome": "generic_name", "total": 1,
            "bulas": [{k: [v] for k, v in secoes.items()}]}


def test_sem_bula_e_indeterminado_e_nunca_sem_alerta():
    veredito, valor = bula.veredito({"bulas": [], "total": 0}, [])
    assert veredito == "indeterminado" and "sem bula" in valor


def test_secao_forte_reprova():
    resultado = _bula(boxed_warning="PERIPHERAL NEUROPATHY has been reported")
    achados = bula.achados(resultado["bulas"], ["neuropathy"])
    veredito, valor = bula.veredito(resultado, achados)
    assert veredito == "reprova" and "boxed_warning" in valor


def test_indicacao_desconta_mesmo_havendo_secao_forte():
    """Pregabalina e duloxetina batem em secao forte, mas sao INDICADAS para a
    condicao. Era 2 dos 4 falsos positivos da regra antiga."""
    resultado = _bula(warnings="neuropathic pain may worsen",
                      indications_and_usage="indicated for neuropathic pain")
    achados = bula.achados(resultado["bulas"], ["neuropathic pain"])
    veredito, valor = bula.veredito(resultado, achados)
    assert veredito == "passa" and "indicado" in valor


def test_secao_fraca_nao_conta_mas_fica_registrada():
    """Parestesia em `adverse_reactions` e quase universal e produziu os outros
    falsos positivos. O achado nao desaparece; ele so nao reprova."""
    resultado = _bula(adverse_reactions="paresthesia was observed in 3%")
    achados = bula.achados(resultado["bulas"], ["paresthesia"])
    assert [a["forca"] for a in achados] == ["fraca"]
    veredito, valor = bula.veredito(resultado, achados)
    assert veredito == "passa" and "fraca" in valor


def test_termo_casa_por_palavra_inteira():
    """Sem isso o falso positivo nao seria nem conferivel."""
    resultado = _bula(warnings="no neuritis-like events")
    assert bula.achados(resultado["bulas"], ["uritis"]) == []


def test_bulas_sao_agregadas_e_qualquer_uma_conta():
    """`limit=1` devolvia uma bula arbitraria entre centenas, e a que traz a
    tarja preta do ciprofloxacino nao e a primeira."""
    resultado = {"campo_do_nome": "generic_name", "total": 2, "bulas": [
        {"warnings": ["nothing relevant here"]},
        {"boxed_warning": ["PERIPHERAL NEUROPATHY"]}]}
    achados = bula.achados(resultado["bulas"], ["neuropathy"])
    assert bula.veredito(resultado, achados)[0] == "reprova"


# --------------------------------------------------------------------------
# expressao: casamento de tecido, e o alvo sem dado
# --------------------------------------------------------------------------

def _linha_de_tecido(tecido, tecido_id, mediana, ancestrais=()):
    return {"tecido_id": tecido_id, "tecido": tecido, "mediana": mediana,
            "maximo": mediana, "unidade": "TPM", "fonte": "gtex",
            "ancestrais": list(ancestrais)}


def test_tecido_nao_casa_por_trecho_do_nome():
    """Medido: um padrao de nome para "nervoso" casou `cortex of kidney`."""
    rim = _linha_de_tecido("cortex of kidney", "UBERON_0001225", 12.0)
    assert expressao.casa_tecido(rim, ["cortex"]) is False
    assert expressao.casa_tecido(rim, ["cortex of kidney"]) is True


def test_tecido_casa_por_ancestralidade():
    cortex = _linha_de_tecido("frontal cortex", "UBERON_0001870", 255.3,
                              ancestrais=["UBERON_0000955"])
    assert expressao.casa_tecido(cortex, ["UBERON_0000955"]) is True


def test_identificador_uberon_e_normalizado():
    assert expressao.normalizar_tecido("uberon:0000955") == "UBERON_0000955"


def test_alvo_sem_expressao_e_indeterminado_e_nao_expressao_zero():
    """Alvo nao humano nao tem GTEx, e isso nao e expressao zero."""
    colunas = expressao.anotar([], ["UBERON_0000955"], corte=1.0)
    assert colunas["expressao_veredito"] == "indeterminado"
    assert colunas["expressao_maxima"] == ""


def test_expressao_baixa_marca_e_explica_o_sert():
    """O portao de ">= 1 TPM em tecido nervoso" reprovaria SERT e NET: os
    alvos da duloxetina e da amitriptilina."""
    sert = _linha_de_tecido("frontal cortex", "UBERON_0001870", 0.45)
    colunas = expressao.anotar([sert], ["UBERON_0001870"], corte=1.0)
    assert colunas["expressao_veredito"] == "marca"
    assert "SERT" in colunas["expressao_observacao"]


def test_sem_tecido_pedido_o_veredito_e_indeterminado():
    cortex = _linha_de_tecido("frontal cortex", "UBERON_0001870", 255.3)
    assert expressao.anotar([cortex], [], corte=1.0)["expressao_veredito"] == (
        "indeterminado")


# --------------------------------------------------------------------------
# esquema: formato reversivel
# --------------------------------------------------------------------------

def test_nan_do_admet_sai_como_celula_vazia():
    assert esquema.formatar(float("nan")) == ""
    assert esquema.formatar(math.inf) == ""


def test_booleano_usa_a_convencao_do_repositorio():
    assert esquema.formatar(True) == "sim" and esquema.formatar(False) == "nao"


def test_escritor_recusa_coluna_fora_do_esquema(tmp_path):
    """Coluna nova sem passar pelo esquema viraria celula perdida."""
    caminho = tmp_path / "t.tsv"
    with esquema.EscritorTSV(caminho, ["a"]) as escritor:
        with pytest.raises(ValueError, match="fora do esquema"):
            escritor.escrever({"a": 1, "b": 2})


def test_ida_e_volta_preserva_o_valor(tmp_path):
    caminho = tmp_path / "t.tsv"
    with esquema.EscritorTSV(caminho, ["nome", "trecho"]) as escritor:
        escritor.escrever({"nome": "x", "trecho": "linha 1\tcom tab\ne quebra"})
    (linha,) = esquema.ler_tsv(caminho)
    assert linha["trecho"] == "linha 1\tcom tab\ne quebra"


def test_colunas_do_portao_vem_da_politica_nao_de_uma_constante():
    colunas = esquema.colunas_portao(["smiles", "alerta_bula"])
    assert colunas[:4] == ["ligante", "ligante_inchikey_parental", "aprovado",
                           "motivos"]
    assert "alerta_bula_veredito" in colunas and "alerta_bula_valor" in colunas
    assert colunas[-1] == "politica"
