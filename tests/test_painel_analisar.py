# -*- coding: utf-8 -*-
"""
T6 -- a aba Analisar: so le, e diz de onde veio.

As funcoes de `vistas` sao puras sobre DataFrame: e por isso que as regras que
importam -- a ordem das faixas, `indeterminado` que nunca vira aprovacao, o
trecho literal da bula que nao pode ser resumido -- podem ser presas aqui, sem
servidor.

Tres estados que a tela tem de distinguir, e nenhum deles e erro: tabela
ausente (nao rodou), tabela vazia (rodou e nada ranqueou) e tabela parcial
(esta rodando agora).
"""
import json

import pandas as pd
import pytest

from mods import caminhos
from mods.painel import vistas
from mods.portao.decisao import VEREDITOS
from mods.ranque.esquema import MOTIVO_FAIXA_3
from mods.ranque.faixa import FAIXAS


def _ranking() -> pd.DataFrame:
    return pd.DataFrame([
        {"posicao": "1", "ligante": "A", "faixa": "1a", "valor": "8.1",
         "alvo": "TRPV1", "pchembl": "8.1", "p_afirmativo": "",
         "portao_aprovado": "sim", "portao_motivos": ""},
        {"posicao": "2", "ligante": "B", "faixa": "2a", "valor": "0.42",
         "alvo": "TRPV1", "pchembl": "", "p_afirmativo": "0.42",
         "portao_aprovado": "sim", "portao_motivos": ""},
        {"posicao": "3", "ligante": "C", "faixa": "1b", "valor": "7.0",
         "alvo": "SCN9A", "pchembl": "7.0", "p_afirmativo": "",
         "portao_aprovado": "nao", "portao_motivos": "bula"},
    ])


# --------------------------------------------------------------------------
# faixas
# --------------------------------------------------------------------------
def test_resumo_das_faixas_na_ordem_da_evidencia():
    """A ordem e 1a, 1b, 2a, 2b, 3 -- medida antes de prevista, alvo do
    usuario antes de alvo enriquecido. Ordenar por contagem poria a faixa
    mais populosa no topo e perderia a hierarquia que ordena o ranking."""
    resumo = vistas.resumo_das_faixas(_ranking())
    assert list(resumo["faixa"]) == ["1a", "1b", "2a"]


def test_resumo_das_faixas_conta_certo():
    resumo = vistas.resumo_das_faixas(_ranking())
    assert resumo.set_index("faixa")["ligantes"]["1a"] == 1


def test_resumo_das_faixas_de_quadro_vazio_nao_levanta():
    resumo = vistas.resumo_das_faixas(pd.DataFrame(columns=["faixa"]))
    assert len(resumo) == 0


def test_filtrar_faixa_vazia_devolve_tudo():
    """Nenhum filtro marcado significa "tudo", nao "nada": uma tela que comeca
    vazia parece resultado vazio."""
    df = _ranking()
    assert len(vistas.filtrar_faixa(df, [])) == len(df)


def test_filtrar_faixa_seleciona():
    assert vistas.filtrar_faixa(_ranking(), ["1a", "1b"])["ligante"].tolist() \
        == ["A", "C"]


def test_legenda_distingue_medido_de_previsto():
    assert "medida" in vistas.legenda_da_faixa("1a").lower()
    assert "prevista" in vistas.legenda_da_faixa("2a").lower()


def test_legenda_distingue_alvo_do_usuario_de_enriquecido():
    assert "declarou" in vistas.legenda_da_faixa("1a").lower()
    assert "enriquec" in vistas.legenda_da_faixa("1b").lower()


def test_toda_faixa_tem_legenda():
    for faixa in FAIXAS:
        assert vistas.legenda_da_faixa(faixa)


def test_legenda_da_faixa_3_diz_que_nao_ranqueia():
    assert "nao ranque" in vistas.legenda_da_faixa("3").lower().replace("ã", "a")


# --------------------------------------------------------------------------
# portao
# --------------------------------------------------------------------------
def test_indeterminado_nao_e_apresentado_como_aprovado():
    """Criterio nao pedido sai `indeterminado`, nunca `passa`. Traduzi-lo para
    "sem alerta" na tela desfaria a distincao que o portao existe para
    manter."""
    df = pd.DataFrame([{"ligante": "A", "aprovado": "indeterminado",
                        "motivos": ""}])
    tela = vistas.portao_para_tela(df)
    assert tela["situacao"].tolist() == ["indeterminado"]


def test_portao_para_tela_preserva_a_coluna_original():
    df = pd.DataFrame([{"ligante": "A", "aprovado": "sim", "motivos": ""}])
    tela = vistas.portao_para_tela(df)
    assert "aprovado" in tela.columns and tela["aprovado"].iloc[0] == "sim"


def test_portao_para_tela_traduz_sim_e_nao():
    df = pd.DataFrame([{"ligante": "A", "aprovado": "sim", "motivos": ""},
                       {"ligante": "B", "aprovado": "nao", "motivos": "bula"}])
    assert vistas.portao_para_tela(df)["situacao"].tolist() == \
        ["aprovado", "reprovado"]


def test_nenhuma_traducao_do_portao_inventa_veredito():
    """Todo rotulo da tela nasce de um veredito que o portao produz."""
    for rotulo in vistas.SITUACOES.values():
        assert rotulo in ("aprovado", "reprovado", *VEREDITOS)


# --------------------------------------------------------------------------
# nao ranqueados e bula
# --------------------------------------------------------------------------
def test_motivos_dos_nao_ranqueados_conta_por_motivo():
    df = pd.DataFrame([{"ligante": "A", "motivo": MOTIVO_FAIXA_3},
                       {"ligante": "B", "motivo": MOTIVO_FAIXA_3},
                       {"ligante": "C", "motivo": "smiles ilegivel"}])
    resumo = vistas.motivos_dos_nao_ranqueados(df)
    assert resumo.set_index("motivo")["ligantes"][MOTIVO_FAIXA_3] == 2


def test_alertas_preservam_o_trecho_literal():
    """O rastro e o produto: sem o trecho e a secao, "tem alerta" nao se
    confere. A negacao em texto livre continua nao resolvida, e e por isso que
    o trecho existe."""
    trecho = "no evidence of peripheral neuropathy was observed"
    df = pd.DataFrame([{"ligante": "A", "termo": "neuropathy",
                        "secao": "adverse_reactions", "forca": "forte",
                        "trecho": trecho, "n_bulas_consultadas": "166",
                        "n_bulas_com_o_termo": "118"}])
    assert vistas.alertas_por_forca(df)["trecho"].iloc[0] == trecho


def test_alertas_mantem_as_duas_contagens():
    """`n_bulas_consultadas` e a janela; `n_bulas_com_o_termo` e o acervo. A
    diferenca e o peso da evidencia: 118 de 166 contra 1 de 211."""
    df = pd.DataFrame([{"ligante": "A", "termo": "t", "secao": "s",
                        "forca": "forte", "trecho": "x",
                        "n_bulas_consultadas": "166",
                        "n_bulas_com_o_termo": "118"}])
    tela = vistas.alertas_por_forca(df)
    assert {"n_bulas_consultadas", "n_bulas_com_o_termo"} <= set(tela.columns)


def test_alertas_ordena_forte_antes_de_fraco():
    df = pd.DataFrame([{"ligante": "A", "forca": "fraca", "termo": "t",
                        "secao": "s", "trecho": "x"},
                       {"ligante": "B", "forca": "forte", "termo": "t",
                        "secao": "s", "trecho": "x"}])
    assert vistas.alertas_por_forca(df)["ligante"].tolist() == ["B", "A"]


# --------------------------------------------------------------------------
# a casca
# --------------------------------------------------------------------------
def _app(tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest
    monkeypatch.setattr(caminhos, "TSV", tmp_path)
    monkeypatch.setattr(caminhos, "EXECUCAO", tmp_path / ".execucao")
    at = AppTest.from_file(str(caminhos.RAIZ / "painel.py"), default_timeout=30)
    at.run()
    return at


def test_analisar_com_ranking_ausente_avisa_e_nao_quebra(tmp_path, monkeypatch):
    at = _app(tmp_path, monkeypatch)
    assert not at.exception
    textos = " ".join(m.value for m in at.info) + " ".join(
        c.value for c in at.caption)
    assert "ainda nao" in textos.lower()


def test_analisar_com_ranking_so_cabecalho_nao_quebra(tmp_path, monkeypatch):
    """Rodou e nada ranqueou: resultado legitimo, nao erro. Com a lista de
    alvos inteira fora do modelo, e o resultado esperado."""
    pasta = tmp_path / "ranking"
    pasta.mkdir()
    (pasta / "ranking.tsv").write_text("#posicao\tligante\tfaixa\n",
                                       encoding="utf-8")
    at = _app(tmp_path, monkeypatch)
    assert not at.exception


def test_analisar_sem_manifesto_avisa_no_lugar_do_selo(tmp_path, monkeypatch):
    """Tabela sem manifesto nao e citavel, e a tela diz isso -- sem esconder a
    tabela, que o pesquisador acabou de gerar."""
    pasta = tmp_path / "ranking"
    pasta.mkdir()
    (pasta / "ranking.tsv").write_text(
        "#posicao\tligante\tfaixa\n1\tA\t1a\n", encoding="utf-8")
    at = _app(tmp_path, monkeypatch)
    assert not at.exception
    assert any("manifesto" in w.value.lower() for w in at.warning)


def test_analisar_mostra_o_manifesto_da_etapa(tmp_path, monkeypatch):
    pasta = tmp_path / "ranking"
    pasta.mkdir()
    (pasta / "ranking.tsv").write_text(
        "#posicao\tligante\tfaixa\n1\tA\t1a\n", encoding="utf-8")
    (pasta / "manifesto.json").write_text(
        json.dumps({"versoes": {"ligantes_md5": "abc123"}}), encoding="utf-8")
    at = _app(tmp_path, monkeypatch)
    assert not at.exception
    assert not any("manifesto" in w.value.lower() for w in at.warning)


# ===========================================================================
# Correcoes da revisao do ramo
# ===========================================================================
def _tmp_de(pasta, nome="ranking.tsv", texto="#posicao\tligante\n1\tA\n"):
    from mods.evidencias import atomico
    tmp = atomico.temporario(pasta / nome)
    tmp.write_text(texto, encoding="utf-8")
    return tmp


def test_analisar_diz_que_a_etapa_esta_rodando(tmp_path, monkeypatch):
    """I2: havia so duas mensagens. A terceira -- "esta rodando agora" --
    recebia o texto de "a etapa pode ter rodado sem a opcao que gera esta
    tabela", que e falso e manda o pesquisador mexer na configuracao."""
    pasta = tmp_path / "ranking"
    pasta.mkdir()
    _tmp_de(pasta)
    at = _app(tmp_path, monkeypatch)
    assert not at.exception
    textos = " ".join(m.value.lower() for m in at.info)
    assert "rodando" in textos or "em curso" in textos
    assert "sem a opcao" not in textos


def test_analisar_avisa_que_a_tabela_mostrada_e_a_anterior(tmp_path,
                                                           monkeypatch):
    """Tabela publicada E execucao reescrevendo a etapa: a tela mostrava os
    numeros velhos sem dizer que estao sendo substituidos."""
    pasta = tmp_path / "ranking"
    pasta.mkdir()
    (pasta / "ranking.tsv").write_text("#posicao\tligante\tfaixa\n1\tA\t1a\n",
                                       encoding="utf-8")
    (pasta / "manifesto.json").write_text(
        json.dumps({"versoes": {"ligantes_md5": "abc"}}), encoding="utf-8")
    _tmp_de(pasta)
    at = _app(tmp_path, monkeypatch)
    assert not at.exception
    avisos = " ".join(w.value.lower() for w in at.warning)
    assert "anterior" in avisos


def test_tabela_sem_cabecalho_nao_derruba_a_tela(tmp_path, monkeypatch):
    """Minor 3: `ler_tsv` levanta ValueError sem o `#` do cabecalho, e o
    traceback subia para a tela. Review Focus 5: nunca stack trace."""
    pasta = tmp_path / "ranking"
    pasta.mkdir()
    (pasta / "ranking.tsv").write_text("posicao\tligante\n1\tA\n",
                                       encoding="utf-8")
    at = _app(tmp_path, monkeypatch)
    assert not at.exception
    assert any("cabe" in e.value.lower() or "formato" in e.value.lower()
               for e in at.error)


def test_tabela_com_linha_de_largura_errada_nao_derruba_a_tela(tmp_path,
                                                              monkeypatch):
    pasta = tmp_path / "ranking"
    pasta.mkdir()
    (pasta / "ranking.tsv").write_text("#a\tb\n1\t2\t3\n", encoding="utf-8")
    at = _app(tmp_path, monkeypatch)
    assert not at.exception
    assert at.error


def test_tabela_de_zero_bytes_nao_e_apresentada_como_vazia(tmp_path,
                                                           monkeypatch):
    """Arquivo de 0 bytes nao e "rodou e nada qualificou": nao e um TSV do
    projeto, e dizer que e resultado legitimo esconde um arquivo truncado."""
    pasta = tmp_path / "ranking"
    pasta.mkdir()
    (pasta / "ranking.tsv").write_bytes(b"")
    at = _app(tmp_path, monkeypatch)
    assert not at.exception
    assert at.error


def test_temporario_da_etapa_sobrevive_a_publicacao(tmp_path, monkeypatch):
    """Minor 2: o `.tmp` pode virar o destino entre o `glob` e o `stat`, e o
    `FileNotFoundError` subia como traceback -- no instante da publicacao, que
    e exatamente quando o pesquisador esta olhando."""
    from mods.painel import leitura
    monkeypatch.setattr(caminhos, "TSV", tmp_path)
    pasta = tmp_path / "etapa1"
    pasta.mkdir()
    tmp = _tmp_de(pasta, "evidencias.tsv")
    stat_real = type(tmp).stat

    def some(self, *a, **k):
        if self.name.endswith(".tmp"):
            raise FileNotFoundError(2, "publicado no meio")
        return stat_real(self, *a, **k)
    monkeypatch.setattr(type(tmp), "stat", some)
    assert leitura.temporario_da_etapa("etapa1") is None
