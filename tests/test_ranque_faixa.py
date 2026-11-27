# -*- coding: utf-8 -*-
"""Faixa de evidencia e chave de ordem da fase 3.

Desenho: docs/superpowers/specs/2026-11-26-fase3-ranqueamento-design.md
"""

import pytest

from mods.evidencias.esquema import TIPO_MEDIDO, TIPO_PREDITO
from mods.ranque.faixa import FAIXAS, chave, classificar


def _linha(tipo, origem, **extra):
    """Uma linha de `evidencias.tsv`, com os campos que a faixa le."""
    linha = {"#ligante": "x", "ligante_inchikey_parental": "K" * 27,
             "alvo_chembl_id": "CHEMBL228", "alvo": "serotonin transporter",
             "organismo": "Homo sapiens", "origem_alvo": origem,
             "tipo_evidencia": tipo, "valor": "", "unidade": "",
             "pchembl": "", "n_moleculas_treino": "900", "no_modelo": "sim"}
    linha.update(extra)
    return linha


def _med(origem="usuario", pchembl="8.0", **extra):
    return _linha(TIPO_MEDIDO, origem, pchembl=pchembl, unidade="nM", **extra)


def _prev(origem="usuario", valor="0.5", **extra):
    return _linha(TIPO_PREDITO, origem, valor=valor, **extra)


def _uma_de_cada():
    """Listas de evidencia de 1a, 1b, 2a, 2b e 3, nessa ordem."""
    return [[_med(origem="usuario", pchembl="5.0")],
            [_med(origem="open_targets", pchembl="9.9")],
            [_prev(origem="usuario", valor="0.10")],
            [_prev(origem="reactome", valor="0.99")],
            []]


def test_medido_em_alvo_do_usuario_e_faixa_1a():
    r = classificar([_med(origem="usuario", pchembl="8.2")], p_negativo=0.9)
    assert r["faixa"] == "1a" and r["valor"] == 8.2


def test_medido_em_alvo_do_enriquecimento_e_faixa_1b():
    for origem in ("open_targets", "reactome"):
        assert classificar([_med(origem=origem, pchembl="7.0")],
                           p_negativo=None)["faixa"] == "1b"


def test_canal_negativo_nao_entra_no_valor_da_faixa_1():
    # pchembl e medido, p_negativo e previsto: multiplicar misturaria os dois
    a = classificar([_med(origem="usuario", pchembl="8.2")], p_negativo=0.0)
    b = classificar([_med(origem="usuario", pchembl="8.2")], p_negativo=0.99)
    assert a["valor"] == b["valor"] == 8.2


def test_previsto_multiplica_pelo_canal_negativo():
    r = classificar([_prev(origem="usuario", valor="0.90")], p_negativo=0.5)
    assert r["faixa"] == "2a" and r["valor"] == pytest.approx(0.45)


def test_previsto_sem_canal_negativo_ordena_pelo_p_sozinho():
    r = classificar([_prev(origem="usuario", valor="0.90")], p_negativo=None)
    assert r["valor"] == pytest.approx(0.90)


def test_sem_evidencia_em_alvo_afirmativo_e_faixa_3_sem_valor():
    r = classificar([], p_negativo=0.1)
    assert r["faixa"] == "3" and r["valor"] is None


def test_agrega_varios_alvos_por_maximo_nao_por_soma():
    r = classificar([_med(origem="usuario", pchembl="6.0"),
                     _med(origem="usuario", pchembl="7.5")], p_negativo=None)
    assert r["valor"] == 7.5


def test_medida_sem_pchembl_nao_ordena_e_nao_rebaixa_para_faixa_3():
    # 1 das 49 linhas medidas do evidencias_controle.tsv e assim
    r = classificar([_med(origem="usuario", pchembl=""),
                     _prev(origem="usuario", valor="0.80")], p_negativo=None)
    assert r["faixa"] == "2a"
    assert r["medida_sem_pchembl"] == "sim"
    assert r["valor"] == pytest.approx(0.80)


def test_medida_sem_pchembl_sozinha_nao_e_faixa_1():
    r = classificar([_med(origem="usuario", pchembl="")], p_negativo=None)
    assert r["faixa"] == "3" and r["medida_sem_pchembl"] == "sim"


def test_outra_linha_medida_com_pchembl_mantem_na_faixa_1():
    r = classificar([_med(origem="usuario", pchembl=""),
                     _med(origem="usuario", pchembl="9.1")], p_negativo=None)
    assert r["faixa"] == "1a" and r["valor"] == 9.1


def test_alvo_afirmativo_fora_do_modelo_continua_faixa_1():
    # sem previsao nenhuma: medido basta
    r = classificar([_med(origem="usuario", pchembl="8.0", no_modelo="nao")],
                    p_negativo=None)
    assert r["faixa"] == "1a"


def test_chave_ordena_as_faixas_na_ordem_do_desenho():
    # Chave ASCENDENTE: `sorted` sem reverse poe 1a primeiro e a faixa 3 por
    # ultimo. A lista de `_uma_de_cada` ja esta na ordem do desenho.
    ordem = [chave(classificar(l, p_negativo=None), "K") for l in _uma_de_cada()]
    assert ordem == sorted(ordem)


def test_alvo_que_decidiu_o_valor_sai_nomeado():
    # Formato "ID nome": um escore sem o alvo responsavel nao e auditavel, e o
    # ID sozinho nao se le. Quando nao ha ID, sai so o nome.
    r = classificar([_med(origem="usuario", pchembl="6.0",
                          alvo_chembl_id="CHEMBL228", alvo="SERT"),
                     _med(origem="usuario", pchembl="7.5",
                          alvo_chembl_id="CHEMBL222", alvo="NET")],
                    p_negativo=None)
    assert r["alvo"] == "CHEMBL222 NET"


def test_alvo_sem_id_sai_so_com_o_nome():
    r = classificar([_med(alvo_chembl_id="", alvo="cst-II")], p_negativo=None)
    assert r["alvo"] == "cst-II"


def test_faixa_3_nunca_vem_antes_de_uma_faixa_com_evidencia():
    vazia = chave(classificar([], p_negativo=None), "AAA")
    fraca = chave(classificar([_prev(origem="reactome", valor="0.001")],
                              p_negativo=0.99), "ZZZ")
    assert fraca < vazia


def test_empate_de_valor_desempata_por_inchikey_crescente():
    a = chave(classificar([_med(pchembl="8.0")], p_negativo=None), "AAA")
    z = chave(classificar([_med(pchembl="8.0")], p_negativo=None), "ZZZ")
    assert a < z


def test_faixas_sao_as_cinco_do_desenho():
    assert FAIXAS == ("1a", "1b", "2a", "2b", "3")


# --------------------------------------------------------------------------
# O vocabulario tem de ser o da etapa 1
# --------------------------------------------------------------------------



def test_usa_o_vocabulario_da_etapa_1_e_nao_um_sinonimo():
    """A etapa 1 escreve `predito`. Escrever `previsto` aqui matava o caminho
    inteiro: nenhuma linha prevista era reconhecida, as faixas previstas eram
    inalcancaveis, e o composto com previsao saia como "sem evidencia"."""
    linha = _linha(TIPO_PREDITO, "usuario", valor="0.90")
    r = classificar([linha], p_negativo=None)
    assert r["faixa"] == "2a" and r["valor"] == pytest.approx(0.90)


def test_um_sinonimo_nao_e_aceito_em_silencio():
    # "previsto" nao e o token da etapa 1: ele nao pode virar faixa prevista
    r = classificar([_linha("previsto", "usuario", valor="0.90")],
                    p_negativo=None)
    assert r["faixa"] == "3"


def test_os_tokens_vem_do_esquema_da_etapa_1():
    from mods.ranque import faixa as modulo
    assert modulo.MEDIDO == TIPO_MEDIDO
    assert modulo.PREDITO == TIPO_PREDITO


def test_pchembl_nao_finito_nao_ordena():
    """`nan` e `inf` sao texto que converte, e uma chave `nan` rompe a
    transitividade do `sorted` -- o arquivo publicado sai numa ordem que nao e
    a ordem das faixas. O gemeo desta funcao na etapa 1 (`medido.py`) ja tinha
    o guarda de finitude; a copia o havia perdido."""
    for ruim in ("nan", "inf", "-inf", "NaN", "Infinity"):
        r = classificar([_med(pchembl=ruim)], p_negativo=None)
        assert r["faixa"] == "3", ruim
        assert r["medida_sem_pchembl"] == "sim", ruim


def test_p_predito_nao_finito_nao_ordena():
    r = classificar([_prev(valor="nan")], p_negativo=None)
    assert r["faixa"] == "3"


def test_a_ordem_nao_e_rompida_por_um_valor_nao_finito():
    valores = ("9.5", "nan", "inf", "7.0", "5.0")
    chaves = {v: chave(classificar([_med(pchembl=v)], p_negativo=None), "K")
              for v in valores}
    ordem = sorted(valores, key=lambda v: chaves[v])
    # os comparaveis em ordem decrescente de pchembl, os nao finitos no fim
    assert ordem[:3] == ["9.5", "7.0", "5.0"]
    assert set(ordem[3:]) == {"nan", "inf"}
    # e a ordenacao e determinista: sem chave `nan`, `sorted` e transitivo
    assert sorted(chaves.values()) == sorted(sorted(chaves.values()))


# ===========================================================================
# Espectro de alvos: a fase 3 com alvo fora do modelo
# ===========================================================================
# A razao de ser da decisao 2: o alvo entra na lista sem receber previsao, e a
# evidencia MEDIDA nele ranqueia em `1b` -- acima de qualquer previsao, em
# qualquer alvo. Se isto nao valesse, deixar o alvo entrar nao compraria nada.
def test_medida_em_alvo_do_reactome_fora_do_modelo_ranqueia_em_1b():
    r = classificar([_med(origem="reactome", pchembl="7.4", no_modelo="nao")],
                    p_negativo=None)
    assert r["faixa"] == "1b"


def test_essa_medida_vence_previsao_em_alvo_pedido():
    """1b vem antes de 2a: medida em alvo enriquecido ordena acima de previsao
    em alvo declarado. E a hierarquia inteira do ranking num par de linhas."""
    medida_enriquecida = chave(
        classificar([_med(origem="reactome", pchembl="5.0", no_modelo="nao")],
                    p_negativo=None), "K1")
    previsao_pedida = chave(
        classificar([_prev(origem="usuario", valor="0.99")],
                    p_negativo=None), "K2")
    assert medida_enriquecida < previsao_pedida


def test_nada_na_fase_3_le_n_moleculas_treino():
    """O alvo fora do modelo sai com `n_moleculas_treino` VAZIO, e o Review
    Focus do plano temia que algum lugar o lesse como numero.

    Os dois testes que eu tinha escrito aqui nao mediam isso: um repetia o
    caso de evidencia vazia e o outro passava o campo a `classificar`, que
    nunca o le -- passavam por vacuidade. A propriedade de verdade e que a
    coluna e so de ESCRITA em todo o pipeline. Se alguem passar a le-la, este
    teste cai e a tolerancia tem de ser provada no lugar certo."""
    import pathlib
    import re
    leitores = []
    for arquivo in pathlib.Path("mods").rglob("*.py"):
        for numero, linha in enumerate(
                arquivo.read_text(encoding="utf-8").splitlines(), start=1):
            if "n_moleculas_treino" not in linha:
                continue
            # escrita (chave de dict, nome de coluna) nao conta; leitura sim.
            # `d["x"] = v` e atribuicao, nao leitura -- e foi o que o primeiro
            # padrao deste teste confundiu.
            if re.search(r'\[\s*["\']n_moleculas_treino["\']\s*\]\s*=', linha):
                continue
            if re.search(r'(\.get\s*\(|\[)\s*["\']n_moleculas_treino', linha):
                leitores.append(f"{arquivo}:{numero}")
    assert not leitores, (
        "alguem passou a LER n_moleculas_treino: prove que o campo vazio do "
        f"alvo fora do modelo e tolerado ali -- {leitores}")
