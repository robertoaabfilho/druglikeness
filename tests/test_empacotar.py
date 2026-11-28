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
    rel = empacotar.montar(destino, completo=False)
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
    empacotar.montar(destino, completo=False)
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
    empacotar.montar(destino, completo=False)
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
    empacotar.montar(destino, completo=False)
    with zipfile.ZipFile(destino) as z:
        raizes = {n.split("/")[0] for n in z.namelist()}
    assert len(raizes) == 1, raizes


def test_o_acervo_do_openfda_nao_entra_no_pacote_ENXUTO():
    """Este teste AFIRMAVA que o acervo nunca entra, e isso deixou de valer: o
    padrao passou a ser o pacote COMPLETO, porque o pacote existe para rodar e
    REPRODUZIR na outra maquina. A exclusao sobrevive apenas no `--enxuto`,
    por pedido explicito."""
    nomes = {i["destino"] for i in empacotar.conteudo(completo=False)}
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
    pasta = caminhos.RAIZ / "dados/openfda"
    # No PACOTE essa pasta existe com a procedencia e NADA mais: nao ha acervo
    # para ficar de fora, e exigir a citacao ali faria a suite do pacote
    # falhar -- e e a suite verde do pacote que prova que ele chegou inteiro.
    tem_acervo = pasta.is_dir() and any(
        a.is_file() and a.name != "procedencia.json" for a in pasta.rglob("*"))
    if not tem_acervo:
        pytest.skip("sem acervo baixado nesta maquina")
    # no pacote ENXUTO: e so la que o acervo fica fora
    fora = empacotar.omitidos(completo=False)
    citados = [o for o in fora if "openfda" in o["arquivo"]]
    assert citados, fora
    # o motivo tem de carregar o COMANDO que resolve, qualquer que seja ele --
    # o teste nao fixa o texto, fixa que ha um comando
    assert all("python" in o["motivo"] for o in citados), citados


# ===========================================================================
# O padrao e TUDO: o pacote tem de rodar e REPRODUZIR na outra maquina
# ===========================================================================
def test_o_padrao_leva_tudo_o_que_existe_nesta_maquina():
    """A correcao de um erro meu: eu escolhia o que valia a pena levar e
    AVISAVA do resto. Aviso nao roda o pipeline na outra maquina.

    Se o dado existe aqui, ele vai. Reprodutibilidade nao se negocia por
    tamanho de arquivo -- e quem quiser um pacote menor pede com `--enxuto`."""
    from mods import caminhos
    nomes = {i["destino"] for i in empacotar.conteudo()}
    for obrigatorio, existe in (
            ("input/ligantes.tsv", (caminhos.RAIZ / "input/ligantes.tsv").is_file()),
            ("input/evidencias_chembl.csv",
             (caminhos.RAIZ / "input/evidencias_chembl.csv").is_file())):
        if existe:
            assert obrigatorio in nomes, obrigatorio


def test_o_padrao_leva_o_indice_do_openfda_se_existe():
    """E o que a etapa 2 LE. Sem ele, a outra maquina nao reproduz o criterio
    que mais reprova -- ela cairia na rede, com outro teto e outra data de
    carga, e o manifesto registraria outra execucao."""
    import pytest

    from mods import caminhos
    indice = caminhos.RAIZ / "dados/openfda/openfda.sqlite"
    if not indice.is_file():
        pytest.skip("sem acervo nesta maquina")
    nomes = {i["destino"] for i in empacotar.conteudo()}
    assert "dados/openfda/openfda.sqlite" in nomes


def test_enxuto_deixa_o_que_e_grande_de_fora_e_diz(tmp_path):
    """O pacote pequeno continua existindo, por pedido explicito -- e ai o que
    ficou fora e dito, com o comando que o traz."""
    nomes = {i["destino"] for i in empacotar.conteudo(completo=False)}
    assert not any(n.startswith("dados/openfda/openfda.sqlite")
                   for n in nomes)
    fora = empacotar.omitidos(completo=False)
    assert all(o["motivo"] for o in fora)


def test_o_cache_de_rede_fica_fora_mesmo_no_completo():
    """A unica exclusao que sobrevive ao `--completo`, e nao e por tamanho: o
    cache e indexado por versao de servico, e servir cache velho noutra
    maquina da resultado de uma versao que nao e a de lá. Correcao, nao
    economia."""
    nomes = {i["destino"] for i in empacotar.conteudo()}
    assert not any(n.startswith("json/cache_rede/") for n in nomes)


def test_o_manifesto_do_pacote_diz_se_e_completo(tmp_path):
    import json
    import zipfile
    destino = tmp_path / "p.zip"
    empacotar.montar(destino, completo=False)
    with zipfile.ZipFile(destino) as z:
        m = json.loads(z.read("druglikeness/PACOTE.json"))
    assert m["completo"] is False


def test_o_leia_me_do_completo_nao_manda_baixar_nada():
    """Num pacote completo, instrucao para baixar o acervo seria ruido -- e,
    pior, faria o pesquisador baixar 1,7 GB que ele ja tem na mao.

    O texto e conferido pela FUNCAO, nao compactando um pacote completo: sao
    ~4 GB, e uma suite que compacta isso tres vezes nao roda. O conteudo do
    completo se verifica pela lista (`conteudo()`), que nao toca em disco."""
    texto = empacotar._leia_me([], [], com_modelo=True, com_acervo=True)
    assert "veio no pacote" in texto
    # a mencao a flag existe, e e para DESACONSELHAR: baixar de novo daria
    # outra carga do openFDA, e o resultado deixaria de ser o mesmo
    assert "não** rode" in texto or "nao** rode" in texto
    assert "acervo do openFDA" in texto


def test_o_leia_me_do_enxuto_manda_baixar_o_acervo():
    texto = empacotar._leia_me(
        [], [{"arquivo": "dados/openfda/", "motivo": "traga com `python "
              "main.py --baixar-acervo`", "mb": 3700.0}],
        com_modelo=True, com_acervo=False)
    assert "--baixar-acervo" in texto


def test_nenhum_teste_compacta_o_pacote_completo():
    """Guarda: `montar(completo=True)` compacta ~4 GB. Se alguem a chamar num
    teste, a suite deixa de rodar -- foi o que aconteceu, duas vezes."""
    import re
    from pathlib import Path
    fonte = Path(__file__).read_text(encoding="utf-8")
    for achado in re.finditer(r"empacotar\.montar\(", fonte):
        trecho = fonte[achado.end():achado.end() + 120]
        assert "completo=False" in trecho, trecho.splitlines()[0]
