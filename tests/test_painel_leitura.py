# -*- coding: utf-8 -*-
"""
T1 -- leitura das tabelas publicadas, com cache por (caminho, mtime).

O painel le do sistema de arquivos: e isso que faz execucao de terminal
aparecer no painel e execucao do painel ser lida por um script. Estes testes
prendem o contrato dessa leitura, inclusive o que NAO pode mudar no caminho:
as aspas e as tabulacoes escapadas de dentro de um nome.
"""
import json

import pytest

from mods import caminhos
from mods.painel import leitura


def _tabela(caminho, linhas):
    caminho.parent.mkdir(parents=True, exist_ok=True)
    caminho.write_text("".join(linhas), encoding="utf-8")


def test_ler_tabela_devolve_none_para_arquivo_ausente(tmp_path):
    assert leitura.ler_tabela(tmp_path / "nao_existe.tsv") is None


def test_ler_tabela_preserva_aspas_no_nome(tmp_path):
    """A aspa e parte do nome: `5"-Methoxyhexahydrocurcumin` e um ligante real
    da base, e foi exatamente o que a citacao do CSV corrompeu antes."""
    p = tmp_path / "t.tsv"
    _tabela(p, ['#ligante\tvalor\n', '5"-Methoxyhexahydrocurcumin\t7.1\n'])
    df = leitura.ler_tabela(p)
    assert df["ligante"].iloc[0] == '5"-Methoxyhexahydrocurcumin'


def test_ler_tabela_reverte_o_escape_da_tabulacao(tmp_path):
    """`EscritorTSV` grava `\\t` para uma tabulacao dentro do campo. Um leitor
    que nao reverte devolve a barra invertida literal, e o nome fica errado."""
    p = tmp_path / "t.tsv"
    _tabela(p, ["#ligante\tvalor\n", "Acido\\tascorbico\t7.1\n"])
    assert leitura.ler_tabela(p)["ligante"].iloc[0] == "Acido\tascorbico"


def test_ler_tabela_respeita_o_limite(tmp_path):
    p = tmp_path / "t.tsv"
    _tabela(p, ["#a\tb\n"] + [f"{i}\t{i}\n" for i in range(10)])
    assert len(leitura.ler_tabela(p, limite=3)) == 3


def test_ler_tabela_de_arquivo_so_com_cabecalho_devolve_quadro_vazio(tmp_path):
    """Tabela vazia e resultado legitimo (nada ranqueou), nao ausencia: o
    quadro vem vazio mas COM as colunas, para a tela saber o que mostrar."""
    p = tmp_path / "t.tsv"
    _tabela(p, ["#ligante\tvalor\n"])
    df = leitura.ler_tabela(p)
    assert df is not None and len(df) == 0
    assert list(df.columns) == ["ligante", "valor"]


def test_contar_linhas_nao_conta_o_cabecalho(tmp_path):
    p = tmp_path / "t.tsv"
    _tabela(p, ["#a\tb\n", "1\t2\n", "3\t4\n"])
    assert leitura.contar_linhas(p) == 2


def test_contar_linhas_de_arquivo_ausente_devolve_none(tmp_path):
    assert leitura.contar_linhas(tmp_path / "nao_existe.tsv") is None


def test_ler_manifesto_devolve_o_dicionario(tmp_path):
    p = tmp_path / "manifesto.json"
    p.write_text(json.dumps({"versoes": {"ligantes_md5": "abc"}}),
                 encoding="utf-8")
    assert leitura.ler_manifesto(p)["versoes"]["ligantes_md5"] == "abc"


def test_ler_manifesto_de_json_invalido_devolve_none(tmp_path):
    """Manifesto quebrado nao derruba a aba: a tela avisa que falta, que e o
    mesmo tratamento de manifesto ausente. Tabela sem manifesto nao e
    citavel, e e isso que o painel tem de dizer."""
    p = tmp_path / "manifesto.json"
    p.write_text("{nao e json", encoding="utf-8")
    assert leitura.ler_manifesto(p) is None


def test_estado_da_etapa_e_ausente_sem_arquivo(tmp_path, monkeypatch):
    monkeypatch.setattr(caminhos, "TSV", tmp_path)
    estado = leitura.estado_da_etapa("etapa1")
    assert estado["estado"] == "ausente"
    assert estado["linhas"] is None and estado["manifesto"] is None


def test_estado_da_etapa_pronta_traz_o_manifesto_e_a_contagem(tmp_path,
                                                              monkeypatch):
    monkeypatch.setattr(caminhos, "TSV", tmp_path)
    _tabela(tmp_path / "etapa1" / "evidencias.tsv",
            ["#ligante\tvalor\n", "A\t1\n", "B\t2\n"])
    (tmp_path / "etapa1" / "manifesto.json").write_text(
        json.dumps({"versoes": {"ligantes_md5": "abc"}}), encoding="utf-8")
    estado = leitura.estado_da_etapa("etapa1")
    assert estado["estado"] == "pronta"
    assert estado["linhas"] == 2
    assert estado["manifesto"]["versoes"]["ligantes_md5"] == "abc"
    assert estado["mtime"] is not None


def test_estado_da_etapa_pronta_sem_manifesto_ainda_e_pronta(tmp_path,
                                                             monkeypatch):
    """A tabela existe e o manifesto nao: o estado da tabela e `pronta`, e
    quem avisa que ela nao e citavel e a tela. Confundir os dois esconderia a
    tabela que o pesquisador acabou de gerar."""
    monkeypatch.setattr(caminhos, "TSV", tmp_path)
    _tabela(tmp_path / "portao" / "portao.tsv", ["#ligante\n", "A\n"])
    estado = leitura.estado_da_etapa("portao")
    assert estado["estado"] == "pronta" and estado["manifesto"] is None


def test_estado_da_etapa_cobre_as_tres_etapas(tmp_path, monkeypatch):
    monkeypatch.setattr(caminhos, "TSV", tmp_path)
    for nome in caminhos.ETAPAS:
        assert leitura.estado_da_etapa(nome)["estado"] == "ausente"


def test_estado_da_etapa_recusa_nome_desconhecido():
    """A mesma regra das flags e da configuracao: nome que nao existe e erro,
    nunca uma pasta inventada em silencio."""
    with pytest.raises(ValueError, match="etapa9"):
        leitura.estado_da_etapa("etapa9")


def test_o_modulo_importa_sem_streamlit_rodando():
    """O painel e a casca; a biblioteca tem de ser testavel sem servidor. Se
    `@st.cache_data` exigisse contexto de execucao, este import ja teria
    falhado -- e com ele toda a suite."""
    import importlib
    importlib.reload(leitura)
