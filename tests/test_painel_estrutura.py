# -*- coding: utf-8 -*-
"""
A aba lateral do painel: a estrutura do pipeline.

A regra que a torna confiavel: **ela LE a estrutura do codigo**, nao a repete.
As etapas vem da tabela `ETAPAS` do `main.py`, as faixas de `mods/ranque/
faixa.py`, as pastas de `mods/caminhos.py`, os codigos de saida do
orquestrador. Um texto decorado divergiria do pipeline na primeira mudanca, e
uma documentacao que mente e pior que nenhuma.
"""
from mods.painel import estrutura


def test_as_etapas_vem_da_tabela_do_orquestrador():
    import main
    nomes = [e["id"] for e in estrutura.etapas()]
    assert nomes == [e["id"] for e in main.ETAPAS]


def test_cada_etapa_traz_script_saida_e_descricao():
    for etapa in estrutura.etapas():
        assert etapa["script"] and etapa["saida"] and etapa["descricao"]
        assert etapa["manifesto"]


def test_as_faixas_vem_do_modulo_que_as_atribui():
    from mods.ranque.faixa import FAIXAS
    assert [f["faixa"] for f in estrutura.faixas()] == list(FAIXAS)


def test_cada_faixa_tem_a_legenda_e_diz_se_ranqueia():
    por_faixa = {f["faixa"]: f for f in estrutura.faixas()}
    assert por_faixa["1a"]["ranqueia"] is True
    assert por_faixa["3"]["ranqueia"] is False
    assert all(f["significado"] for f in estrutura.faixas())


def test_as_pastas_vem_de_caminhos():
    from mods import caminhos
    pastas = {p["nome"] for p in estrutura.pastas()}
    assert {"input/", "mods/", "json/", "tsv/", "modelo/"} <= pastas
    assert all(p["papel"] for p in estrutura.pastas())


def test_a_regra_que_desempata_input_de_tsv_esta_dita():
    texto = " ".join(p["papel"] for p in estrutura.pastas()).lower()
    assert "procedencia" in texto or "procedência" in texto


def test_os_codigos_de_saida_vem_do_orquestrador():
    codigos = estrutura.codigos_de_saida()
    assert set(codigos) == {0, 1, 2, 3}
    assert "parcial" in codigos[3].lower()


def test_as_dependencias_vem_do_requirements():
    nomes = estrutura.dependencias()
    assert "rdkit" in nomes and "streamlit" in nomes
    assert not any(n.startswith("#") for n in nomes)


def test_a_estrutura_nao_tem_texto_decorado_das_etapas():
    """Se a descricao das etapas fosse copiada para ca, ela divergiria do
    `main.py` na primeira mudanca. Este teste prende a leitura."""
    import inspect
    fonte = inspect.getsource(estrutura)
    assert "COLETAR" not in fonte.replace("# COLETAR", "")
    assert "import main" in fonte


def test_o_painel_mostra_a_estrutura_na_barra_lateral(tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest

    from mods import caminhos
    monkeypatch.setattr(caminhos, "TSV", tmp_path)
    monkeypatch.setattr(caminhos, "EXECUCAO", tmp_path / ".execucao")
    at = AppTest.from_file(str(caminhos.RAIZ / "painel.py"), default_timeout=30)
    at.run()
    assert not at.exception
    texto = " ".join(str(e.value) for e in at.sidebar.markdown)
    for esperado in ("etapa1", "portao", "ranking", "1a", "2b"):
        assert esperado in texto, f"a barra lateral nao menciona {esperado}"


def test_a_barra_lateral_diz_onde_cada_tabela_e_escrita(tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest

    from mods import caminhos
    monkeypatch.setattr(caminhos, "TSV", tmp_path)
    monkeypatch.setattr(caminhos, "EXECUCAO", tmp_path / ".execucao")
    at = AppTest.from_file(str(caminhos.RAIZ / "painel.py"), default_timeout=30)
    at.run()
    texto = " ".join(str(e.value) for e in at.sidebar.markdown)
    assert "tsv/etapa1/evidencias.tsv" in texto
    assert "manifesto.json" in texto
