# -*- coding: utf-8 -*-
"""
A verificacao de dependencias na partida.

O pipeline comeca conferindo o que precisa e instalando o que falta. Duas
decisoes que atravessam este arquivo:

**Instalar mexe no ambiente de quem roda.** Entao: so o que esta declarado em
`requirements.txt`, nunca uma versao que o arquivo nao pede, e **dizendo** o
que vai instalar antes de instalar. Instalacao silenciosa num ambiente que o
pesquisador compartilha com outro trabalho e a forma mais rapida de quebrar o
outro trabalho.

**O pesado nao entra sozinho.** `torch` sao alguns GB e esta em
`requirements-opcionais.txt` justamente porque o pipeline DEGRADA sem ele: a
etapa 1 roda com evidencia medida e diz o que desligou. Baixar GB sem pedir
seria decidir pelo pesquisador num ponto em que ele tem alternativa.

Nenhum teste aqui instala nada: o instalador e injetado.
"""
import pytest

from mods import dependencias


# --------------------------------------------------------------------------
# o que e exigido, e de onde vem
# --------------------------------------------------------------------------
def test_as_exigencias_vem_do_requirements():
    nomes = {e["pacote"] for e in dependencias.exigidas()}
    assert {"rdkit", "pandas", "numpy", "streamlit", "scipy"} <= nomes


def test_a_exigencia_guarda_a_linha_original_do_arquivo():
    """A linha e o que vai para o `pip install`: `numpy>=2` nao pode virar
    `numpy`, senao a instalacao satisfaz o nome e nao a restricao."""
    por_pacote = {e["pacote"]: e for e in dependencias.exigidas()}
    assert por_pacote["numpy"]["requisito"] == "numpy>=2"


def test_o_modulo_importavel_pode_ter_nome_diferente_do_pacote():
    """`rdkit` instala e importa como `rdkit`, mas a regra geral nao vale para
    todo pacote; o mapa existe para os casos em que diverge."""
    por_pacote = {e["pacote"]: e for e in dependencias.exigidas()}
    assert por_pacote["rdkit"]["modulo"] == "rdkit"


def test_as_opcionais_sao_listadas_separadas():
    nomes = {e["pacote"] for e in dependencias.opcionais()}
    assert "torch" in nomes
    assert "torch" not in {e["pacote"] for e in dependencias.exigidas()}


def test_cada_opcional_diz_o_que_desliga():
    """Quem decide se baixa alguns GB precisa saber o que perde se nao baixar."""
    for item in dependencias.opcionais():
        assert item["desliga"], item


# --------------------------------------------------------------------------
# conferir
# --------------------------------------------------------------------------
def test_conferir_acha_o_que_falta():
    estado = dependencias.conferir(presente=lambda m: m != "scipy")
    assert [f["pacote"] for f in estado["faltando"]] == ["scipy"]
    assert estado["completo"] is False


def test_conferir_com_tudo_presente():
    estado = dependencias.conferir(presente=lambda m: True)
    assert estado["faltando"] == [] and estado["completo"] is True


def test_conferir_relata_a_opcional_ausente_sem_chamar_de_falta():
    """Opcional ausente NAO e falta: o pipeline roda degradado e diz o que
    desligou. Trata-la como falta pararia a execucao por algo que nao para."""
    estado = dependencias.conferir(presente=lambda m: m != "torch")
    assert estado["faltando"] == []
    assert [o["pacote"] for o in estado["opcionais_ausentes"]] == ["torch"]
    assert estado["completo"] is True


# --------------------------------------------------------------------------
# instalar
# --------------------------------------------------------------------------
def test_garantir_instala_so_o_que_falta():
    pedidos = []
    dependencias.garantir(presente=lambda m: m != "scipy",
                          instalador=lambda reqs: pedidos.append(list(reqs)))
    assert pedidos == [["scipy"]]


def test_garantir_passa_o_requisito_com_a_restricao():
    pedidos = []
    dependencias.garantir(presente=lambda m: m != "numpy",
                          instalador=lambda reqs: pedidos.append(list(reqs)))
    assert pedidos == [["numpy>=2"]]


def test_garantir_nao_instala_nada_quando_nao_falta_nada():
    pedidos = []
    dependencias.garantir(presente=lambda m: True,
                          instalador=lambda reqs: pedidos.append(list(reqs)))
    assert pedidos == []


def test_garantir_NUNCA_instala_opcional_sozinho():
    """`torch` sao alguns GB, e o pipeline degrada sem ele. Baixar GB sem
    pedir seria decidir pelo pesquisador onde ele tem alternativa."""
    pedidos = []
    dependencias.garantir(presente=lambda m: m != "torch",
                          instalador=lambda reqs: pedidos.append(list(reqs)))
    assert pedidos == []


def test_garantir_so_confere_quando_pedido():
    """`instalar=False` e o modo de quem nao quer que nada mude no ambiente --
    automacao, ambiente compartilhado, maquina sem rede."""
    pedidos = []
    estado = dependencias.garantir(
        presente=lambda m: m != "scipy", instalar=False,
        instalador=lambda reqs: pedidos.append(list(reqs)))
    assert pedidos == []
    assert estado["completo"] is False
    assert estado["instalou"] == []


def test_garantir_devolve_falha_quando_a_instalacao_nao_resolve():
    """Instalou e o pacote continua ausente: a falha e dita, nao engolida.
    Seguir em frente daria um erro de import no meio da etapa 1."""
    def instalador(reqs):
        return None  # finge instalar e nada muda

    estado = dependencias.garantir(presente=lambda m: m != "scipy",
                                   instalador=instalador)
    assert estado["completo"] is False
    assert any("scipy" in m for m in estado["erros"])


def test_garantir_relata_erro_do_instalador_sem_levantar():
    """Sem rede, o pip falha. O pipeline tem de dizer o que falta e o comando
    que resolve, em vez de subir um traceback de rede."""
    def instalador(reqs):
        raise RuntimeError("Network is unreachable")

    estado = dependencias.garantir(presente=lambda m: m != "scipy",
                                   instalador=instalador)
    assert estado["completo"] is False
    assert any("unreachable" in m for m in estado["erros"])
    assert any("pip install" in m for m in estado["erros"])


def test_o_instalador_de_verdade_usa_o_interpretador_corrente():
    """`pip` do PATH pode ser de outro ambiente -- e instalar no ambiente
    errado parece sucesso e nao resolve nada."""
    import inspect
    fonte = inspect.getsource(dependencias.instalar_com_pip)
    assert "sys.executable" in fonte
    assert "-m" in fonte and "pip" in fonte


def test_o_instalador_de_verdade_nao_atualiza_o_que_nao_foi_pedido():
    """A tarefa e suprir o que falta, nao mexer no que esta instalado e
    funcionando. O teste olha o COMANDO montado, nao a prosa do docstring."""
    comandos = []

    def capturar(comando, **k):
        comandos.append(comando)

        class R:
            returncode = 0
            stdout = stderr = ""
        return R()
    import subprocess
    real = subprocess.run
    subprocess.run = capturar
    try:
        dependencias.instalar_com_pip(["scipy"])
    finally:
        subprocess.run = real
    assert comandos and "--upgrade" not in comandos[0]
    assert comandos[0][:4] == [__import__("sys").executable, "-m", "pip",
                               "install"]


# --------------------------------------------------------------------------
# o relato
# --------------------------------------------------------------------------
def test_o_relato_nomeia_o_que_desliga_sem_a_opcional():
    estado = dependencias.conferir(presente=lambda m: m != "torch")
    texto = "\n".join(dependencias.relato(estado))
    assert "torch" in texto
    assert "degrad" in texto.lower() or "desliga" in texto.lower()


def test_o_relato_de_ambiente_completo_e_uma_linha():
    estado = dependencias.conferir(presente=lambda m: True)
    assert len(dependencias.relato(estado)) == 1


def test_o_relato_diz_o_comando_que_resolve():
    estado = dependencias.conferir(presente=lambda m: m != "rdkit")
    texto = "\n".join(dependencias.relato(estado))
    assert "pip install -r requirements.txt" in texto


# ===========================================================================
# na partida do orquestrador e do painel
# ===========================================================================
def test_o_orquestrador_confere_as_dependencias_antes_de_qualquer_etapa(
        monkeypatch, tmp_path):
    """Antes de qualquer etapa: uma execucao real sao ~23 milhoes de pares, e
    um ModuleNotFoundError no meio dela custa horas por algo que custa
    milissegundos conferir."""
    import main
    ordem = []
    monkeypatch.setattr(main.dependencias, "garantir",
                        lambda **k: ordem.append("dependencias") or
                        {"completo": True, "faltando": [],
                         "opcionais_ausentes": [], "instalou": [],
                         "erros": []})
    monkeypatch.setattr(main, "_precheca",
                        lambda *a, **k: ordem.append("precheca") or
                        ["parar aqui"])
    main.main(["--config", "json/config.json"])
    assert ordem[:2] == ["dependencias", "precheca"]


def test_o_orquestrador_para_com_codigo_2_se_falta_dependencia(monkeypatch):
    """Nao se roda meia etapa com biblioteca faltando: o codigo 2 e `erro`, e
    nada e publicado."""
    import main
    monkeypatch.setattr(main.dependencias, "garantir", lambda **k: {
        "completo": False,
        "faltando": [{"pacote": "scipy", "requisito": "scipy",
                      "modulo": "scipy"}],
        "opcionais_ausentes": [], "instalou": [],
        "erros": ["scipy continua ausente"]})
    assert main.main(["--config", "json/config.json"]) == 2


def test_a_flag_sem_instalar_existe_e_so_confere(monkeypatch):
    """Ambiente compartilhado, automacao, maquina sem rede: quem nao quer que
    nada mude no ambiente tem como dizer."""
    import main
    vistos = {}
    monkeypatch.setattr(main.dependencias, "garantir",
                        lambda **k: vistos.update(k) or
                        {"completo": True, "faltando": [],
                         "opcionais_ausentes": [], "instalou": [],
                         "erros": []})
    monkeypatch.setattr(main, "_precheca", lambda *a, **k: ["parar"])
    main.main(["--sem-instalar", "--config", "json/config.json"])
    assert vistos.get("instalar") is False


def test_o_painel_mostra_o_estado_das_dependencias(tmp_path, monkeypatch):
    """O painel NAO instala: ele e uma interface, e instalar por clique mexe
    no ambiente de quem abriu a aba. Ele mostra o estado e o comando."""
    from streamlit.testing.v1 import AppTest

    from mods import caminhos
    monkeypatch.setattr(caminhos, "TSV", tmp_path)
    monkeypatch.setattr(caminhos, "EXECUCAO", tmp_path / ".execucao")
    at = AppTest.from_file(str(caminhos.RAIZ / "painel.py"), default_timeout=30)
    at.run()
    assert not at.exception
    texto = " ".join(str(e.value) for e in at.sidebar.markdown) + " ".join(
        str(e.value) for e in at.sidebar.caption)
    assert "depend" in texto.lower()


def test_o_painel_nao_instala_nada():
    import pathlib
    fonte = pathlib.Path("painel.py").read_text(encoding="utf-8")
    assert "garantir(" not in fonte, "o painel nao deve instalar dependencias"


# ===========================================================================
# O acervo do openFDA: verificar sempre, baixar so se pedido
# ===========================================================================
def test_acervo_bula_e_configuravel():
    """ACHADO: `--acervo-bula` existia como flag da CLI e NAO estava nas chaves
    do config -- pelo orquestrador o acervo era inalcancavel."""
    from mods import config
    assert "acervo_bula" in config.CHAVES["portao"]


def test_toda_flag_da_etapa_2_que_muda_a_saida_esta_nas_chaves():
    """A guarda inversa da que ja existia. A antiga provava que toda chave tem
    flag; faltava provar que toda flag ALCANCAVEL pelo orquestrador esta como
    chave -- foi por essa fresta que `--acervo-bula` passou."""
    import re
    from pathlib import Path

    from mods import caminhos, config
    fonte = (caminhos.RAIZ / "mods/aplicar_portao.py").read_text(
        encoding="utf-8")
    flags = set(re.findall(r'"--([a-z0-9-]+)"', fonte))
    # Flags que sao do uso direto da CLI, nao da configuracao: a repeticao
    # aponta um manifesto anterior, e --help nao e parametro.
    so_da_cli = {"help", "manifesto"}
    faltando = {f for f in flags - so_da_cli
                if f.replace("-", "_") not in config.CHAVES["portao"]}
    assert not faltando, faltando


def test_o_acervo_ausente_e_conferido_antes_de_rodar(tmp_path, monkeypatch):
    """Se a configuracao aponta um acervo que nao existe, isso tem de aparecer
    ANTES da etapa 1 -- nao depois de horas, quando a etapa 2 comeca."""
    import main
    monkeypatch.setattr(main.caminhos, "RAIZ", tmp_path)
    for nome in ("l.tsv", "e.tsv"):
        (tmp_path / nome).write_text("x", encoding="utf-8")
    cfg = {"portao": {"evidencias": "e.tsv", "ligantes": "l.tsv",
                      "saida": "p", "acervo_bula": "dados/openfda/nao.sqlite"}}
    erros = main._precheca([e for e in main.ETAPAS if e["id"] == "portao"], cfg)
    assert any("acervo_bula" in e for e in erros), erros
    assert any("openfda_massa" in e for e in erros), erros


def test_acervo_presente_nao_gera_erro(tmp_path, monkeypatch):
    import main
    monkeypatch.setattr(main.caminhos, "RAIZ", tmp_path)
    for nome in ("l.tsv", "e.tsv", "a.sqlite"):
        (tmp_path / nome).write_text("x", encoding="utf-8")
    cfg = {"portao": {"evidencias": "e.tsv", "ligantes": "l.tsv",
                      "saida": "p", "acervo_bula": "a.sqlite"}}
    assert main._precheca([e for e in main.ETAPAS if e["id"] == "portao"],
                          cfg) == []


def test_baixar_o_acervo_NAO_e_automatico(tmp_path, monkeypatch):
    """3,7 GB e ~20 min, e ha alternativa: a rede com o teto de taxa. A mesma
    regra do `torch` -- baixar GB sem pedir e decidir pelo pesquisador num
    ponto em que ele tem escolha. E, para uma base de produto natural, o acervo
    compra quase nada: 58 de 60 ligantes nao tem bula nenhuma."""
    import main
    baixou = []
    monkeypatch.setattr(main, "_baixar_acervo",
                        lambda *a, **k: baixou.append(1) or True)
    monkeypatch.setattr(main.caminhos, "RAIZ", tmp_path)
    for nome in ("l.tsv", "e.tsv"):
        (tmp_path / nome).write_text("x", encoding="utf-8")
    cfg = {"portao": {"evidencias": "e.tsv", "ligantes": "l.tsv",
                      "saida": "p", "acervo_bula": "dados/openfda/nao.sqlite"}}
    main._precheca([e for e in main.ETAPAS if e["id"] == "portao"], cfg)
    assert baixou == [], "baixou 3,7 GB sem ninguem pedir"


def test_a_flag_de_baixar_o_acervo_existe():
    """Para que caiba num comando so quando o pesquisador QUER."""
    import main
    args = main._argumentos(["--baixar-acervo"])
    assert args.baixar_acervo is True
    assert main._argumentos([]).baixar_acervo is False


def test_com_a_flag_o_acervo_e_baixado_antes_das_etapas(monkeypatch, tmp_path):
    import main
    ordem = []
    monkeypatch.setattr(main.dependencias, "garantir", lambda **k: {
        "completo": True, "faltando": [], "opcionais_ausentes": [],
        "instalou": [], "erros": []})
    monkeypatch.setattr(main, "_baixar_acervo",
                        lambda *a, **k: ordem.append("acervo") or True)
    monkeypatch.setattr(main, "_precheca",
                        lambda *a, **k: ordem.append("precheca") or ["parar"])
    main.main(["--baixar-acervo", "--config", "json/config.json"])
    assert ordem[:2] == ["acervo", "precheca"]
