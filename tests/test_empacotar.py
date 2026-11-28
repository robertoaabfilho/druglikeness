# -*- coding: utf-8 -*-
"""
O pacote para rodar noutra maquina.

A pergunta que define o conteudo: **o que falta numa maquina limpa para o
pipeline rodar?** O codigo e os testes vem do git; o que o git NAO carrega e
exatamente o que o pacote precisa levar -- e o que ele levar tem de vir com a
procedencia, senao o resultado na outra maquina nao e citavel.
"""
from mods import empacotar


def test_o_pacote_leva_o_codigo_e_os_testes():
    itens = empacotar.conteudo()
    nomes = {i["destino"] for i in itens}
    assert "main.py" in nomes and "painel.py" in nomes
    assert any(n.startswith("mods/") for n in nomes)
    assert any(n.startswith("tests/") for n in nomes)
    assert "requirements.txt" in nomes


def test_o_pacote_leva_o_que_o_git_nao_carrega():
    """O modelo e as entradas ficam fora do repositorio; sem eles a outra
    maquina tem o codigo e nada para rodar."""
    nomes = {i["destino"] for i in empacotar.conteudo()}
    assert any(n.startswith("modelo/") for n in nomes)


def test_o_pacote_leva_a_procedencia_das_entradas():
    nomes = {i["destino"] for i in empacotar.conteudo()}
    assert "input/procedencia.json" in nomes


def test_o_pacote_nao_leva_as_saidas_nem_o_cache():
    """`tsv/` e resultado de uma execucao desta maquina; levar junto faria a
    outra maquina parecer ter rodado o pipeline."""
    nomes = {i["destino"] for i in empacotar.conteudo()}
    assert not any(n.startswith("tsv/") for n in nomes)
    assert not any(n.startswith("json/cache_rede/") for n in nomes)
    assert not any("__pycache__" in n for n in nomes)
    assert not any(n.startswith(".git/") for n in nomes)


def test_o_pacote_declara_o_que_NAO_leva(tmp_path):
    """Um pacote que omite em silencio o que nao cabe faz o pesquisador
    descobrir a falta na outra maquina. O manifesto diz o que ficou fora e
    por que."""
    destino = tmp_path / "p.zip"
    rel = empacotar.montar(destino, incluir_dados_grandes=False)
    # A lista pode ser vazia -- noutra maquina (ou no proprio pacote
    # descompactado) nao ha dado grande no disco para ficar de fora. O que nao
    # pode e um item sem motivo: e o motivo que impede o pesquisador de
    # descobrir a falta do outro lado.
    assert "omitidos" in rel
    assert all(o["motivo"] and o["arquivo"] for o in rel["omitidos"])
    # e nesta maquina, onde os dados grandes existem, eles TEM de ser citados
    from mods import caminhos
    if (caminhos.RAIZ / "input/ligantes.tsv").is_file():
        assert any(o["arquivo"] == "input/ligantes.tsv"
                   for o in rel["omitidos"]), rel["omitidos"]


def test_o_zip_abre_e_traz_o_manifesto(tmp_path):
    import json
    import zipfile
    destino = tmp_path / "p.zip"
    empacotar.montar(destino, incluir_dados_grandes=False)
    with zipfile.ZipFile(destino) as z:
        assert z.testzip() is None
        nomes = z.namelist()
        assert any(n.endswith("PACOTE.json") for n in nomes)
        manifesto = json.loads(z.read(
            [n for n in nomes if n.endswith("PACOTE.json")][0]))
    assert manifesto["commit"] and manifesto["arquivos"] > 0
    assert "omitidos" in manifesto


def test_o_zip_traz_as_instrucoes_de_instalacao(tmp_path):
    import zipfile
    destino = tmp_path / "p.zip"
    empacotar.montar(destino, incluir_dados_grandes=False)
    with zipfile.ZipFile(destino) as z:
        leia = [n for n in z.namelist() if n.endswith("LEIA-ME-PRIMEIRO.md")]
        assert leia
        texto = z.read(leia[0]).decode("utf-8")
    assert "pip install -r requirements.txt" in texto
    assert "pytest" in texto


def test_o_pacote_tem_uma_pasta_raiz_so(tmp_path):
    """Descompactar nao pode espalhar 174 arquivos no diretorio corrente."""
    import zipfile
    destino = tmp_path / "p.zip"
    empacotar.montar(destino, incluir_dados_grandes=False)
    with zipfile.ZipFile(destino) as z:
        raizes = {n.split("/")[0] for n in z.namelist()}
    assert len(raizes) == 1, raizes


def test_o_acervo_do_openfda_nao_entra_no_pacote():
    """ACHADO pelo timeout da suite: o empacotador varria `dados/` inteiro, e
    depois do download o acervo do openFDA sao 3,7 GB ali dentro -- o pacote
    passaria de 31 MB para mais de 3 GB, e tres testes tentavam comprimir isso.

    O acervo e dado GERADO pela maquina, como `tsv/` e o cache de rede: a outra
    maquina baixa o seu com um comando. A mesma regra do `.gitignore`."""
    nomes = {i["destino"] for i in empacotar.conteudo()}
    sobraram = [n for n in nomes if n.startswith("dados/openfda/")
                and not n.endswith("procedencia.json")]
    assert not sobraram, sobraram


def test_mas_a_procedencia_do_acervo_entra():
    """Se o acervo existe nesta maquina, a outra precisa saber de que export
    ele era -- a mesma regra de `input/procedencia.json`."""
    import pytest

    from mods import caminhos
    if not (caminhos.RAIZ / "dados/openfda/procedencia.json").is_file():
        pytest.skip("sem acervo baixado nesta maquina")
    nomes = {i["destino"] for i in empacotar.conteudo()}
    assert "dados/openfda/procedencia.json" in nomes


def test_o_pacote_diz_que_o_acervo_ficou_fora():
    import pytest

    from mods import caminhos
    if not (caminhos.RAIZ / "dados/openfda").is_dir():
        pytest.skip("sem acervo baixado nesta maquina")
    fora = empacotar.omitidos()
    citados = [o for o in fora if "openfda" in o["arquivo"]]
    assert citados, fora
    assert all("openfda_massa" in o["motivo"] for o in citados)
