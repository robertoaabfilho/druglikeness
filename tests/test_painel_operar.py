# -*- coding: utf-8 -*-
"""
T5 -- a aba Operar: entradas, parametros, o que rodar, acompanhamento.

Duas regras do desenho que estes testes prendem:

* **O upload diz o que substituiu.** Trocar `input/ligantes.tsv` sem dizer que
  havia outro arquivo ali deixaria o pesquisador sem saber qual base gerou o
  resultado que ele tem na tela.
* **A conferencia roda antes da execucao, e e a MESMA da etapa 1**
  (`mods.evidencias.entrada`). Uma validacao propria do painel divergiria da
  do pipeline, e o painel aprovaria entrada que a etapa recusa.
"""
import pytest

from mods import caminhos
from mods.painel import entradas


@pytest.fixture
def input_limpo(tmp_path, monkeypatch):
    monkeypatch.setattr(caminhos, "INPUT", tmp_path / "input")
    (tmp_path / "input").mkdir()
    return tmp_path / "input"


class _Envio:
    """O que `st.file_uploader` devolve: nome e bytes."""

    def __init__(self, nome, dados: bytes):
        self.name = nome
        self._dados = dados

    def getvalue(self) -> bytes:
        return self._dados


LIGANTES = ("#nome\tsmiles\n"
            "Cardamomin\tCOc1cc(O)cc(O)c1C(=O)/C=C/c1ccccc1\n"
            "Capsaicina\tCOc1cc(CNC(=O)CCCCC=CC(C)C)ccc1O\n")


# --------------------------------------------------------------------------
# receber
# --------------------------------------------------------------------------
def test_receber_grava_no_nome_canonico_do_papel(input_limpo):
    """O arquivo vai para `input/ligantes.tsv`, qualquer que seja o nome do
    upload: e esse o caminho que a configuracao aponta."""
    r = entradas.receber(_Envio("minha_base_v7.tsv", LIGANTES.encode()),
                         "ligantes")
    assert r["destino"] == input_limpo / "ligantes.tsv"
    assert r["destino"].read_text(encoding="utf-8") == LIGANTES


def test_receber_diz_que_nao_substituiu_nada(input_limpo):
    r = entradas.receber(_Envio("l.tsv", LIGANTES.encode()), "ligantes")
    assert r["substituiu"] is False and r["md5_anterior"] is None


def test_receber_diz_o_que_substituiu(input_limpo):
    (input_limpo / "ligantes.tsv").write_text("#nome\tsmiles\nA\tCC\n",
                                              encoding="utf-8")
    r = entradas.receber(_Envio("l.tsv", LIGANTES.encode()), "ligantes")
    assert r["substituiu"] is True
    assert r["md5_anterior"] and r["md5_anterior"] != r["md5"]


def test_receber_traz_a_conferencia_do_arquivo_novo(input_limpo):
    r = entradas.receber(_Envio("l.tsv", LIGANTES.encode()), "ligantes")
    assert r["conferencia"]["aceitos"] == 2


def test_receber_recusa_papel_desconhecido(input_limpo):
    with pytest.raises(ValueError, match="qualquer"):
        entradas.receber(_Envio("x.tsv", b"a"), "qualquer")


def test_receber_nao_deixa_o_arquivo_pela_metade(input_limpo, monkeypatch):
    """Escrita atomica: upload interrompido nao destroi a base anterior."""
    antes = "#nome\tsmiles\nA\tCC\n"
    (input_limpo / "ligantes.tsv").write_text(antes, encoding="utf-8")

    def explode(*a, **k):
        raise OSError("disco cheio")
    monkeypatch.setattr(entradas.atomico, "publicar", explode)
    with pytest.raises(OSError):
        entradas.receber(_Envio("l.tsv", LIGANTES.encode()), "ligantes")
    assert (input_limpo / "ligantes.tsv").read_text(encoding="utf-8") == antes


def test_receber_nao_deixa_temporario_para_tras(input_limpo, monkeypatch):
    def explode(*a, **k):
        raise OSError("disco cheio")
    monkeypatch.setattr(entradas.atomico, "publicar", explode)
    with pytest.raises(OSError):
        entradas.receber(_Envio("l.tsv", LIGANTES.encode()), "ligantes")
    assert list(input_limpo.glob("*.tmp")) == []


# --------------------------------------------------------------------------
# conferir
# --------------------------------------------------------------------------
def test_conferir_conta_os_aceitos(input_limpo):
    p = input_limpo / "ligantes.tsv"
    p.write_text(LIGANTES, encoding="utf-8")
    r = entradas.conferir(p, "ligantes")
    assert r["aceitos"] == 2 and r["erros"] == 0


def test_conferir_conta_smiles_ilegivel_sem_levantar(input_limpo):
    """Molecula ruim nao derruba nada: ela e JULGADA, e o painel diz quantas
    sao antes de a execucao comecar."""
    p = input_limpo / "ligantes.tsv"
    p.write_text(LIGANTES + "Quebrado\tnao_e_smiles((((\n", encoding="utf-8")
    r = entradas.conferir(p, "ligantes")
    assert r["aceitos"] == 2
    assert r["problemas"] >= 1


def test_conferir_relata_coluna_obrigatoria_ausente(input_limpo):
    p = input_limpo / "ligantes.tsv"
    p.write_text("#nome\tcor\nA\tazul\n", encoding="utf-8")
    r = entradas.conferir(p, "ligantes")
    assert r["aceitos"] == 0 and r["erros"] >= 1
    assert any("smiles" in m.lower() for m in r["mensagens"])


def test_conferir_alvos_le_a_lista_simples(input_limpo):
    p = input_limpo / "alvos.txt"
    p.write_text("# comentario\nCHEMBL233\nCHEMBL1951\n", encoding="utf-8")
    assert entradas.conferir(p, "alvos")["aceitos"] == 2


def test_conferir_sintomas_ignora_comentario_e_vazio(input_limpo):
    p = input_limpo / "sintomas.txt"
    p.write_text("# c\n\nneuropathic pain\nallodynia\n", encoding="utf-8")
    assert entradas.conferir(p, "sintomas")["aceitos"] == 2


def test_conferir_arquivo_ausente_devolve_erro_nao_excecao(input_limpo):
    r = entradas.conferir(input_limpo / "nao_existe.tsv", "ligantes")
    assert r["aceitos"] == 0 and r["erros"] >= 1


def test_conferir_limita_as_mensagens(input_limpo):
    """Dezessete mil linhas sem nome produzem dezessete mil problemas; a tela
    mostra as primeiras e a contagem."""
    p = input_limpo / "ligantes.tsv"
    p.write_text("#nome\tsmiles\n" + "".join("\tCC\n" for _ in range(50)),
                 encoding="utf-8")
    r = entradas.conferir(p, "ligantes")
    assert len(r["mensagens"]) <= 10 and r["problemas"] >= 50


def test_estado_das_entradas_diz_o_que_falta(input_limpo):
    (input_limpo / "ligantes.tsv").write_text(LIGANTES, encoding="utf-8")
    estado = entradas.estado()
    assert estado["ligantes"]["existe"] is True
    assert estado["alvos"]["existe"] is False


# --------------------------------------------------------------------------
# a casca
# --------------------------------------------------------------------------
def test_painel_abre_sem_tsv_nem_modelo(tmp_path, monkeypatch):
    """A arvore recem-clonada nao tem `modelo/` nem `tsv/`. O painel tem de
    abrir e dizer o que falta -- nao explodir, que foi o que o `main.py` fazia
    antes do pre-checagem."""
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(str(caminhos.RAIZ / "painel.py"), default_timeout=30)
    at.run()
    assert not at.exception


def test_operar_oferece_as_tres_etapas():
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(str(caminhos.RAIZ / "painel.py"), default_timeout=30)
    at.run()
    etapas = at.multiselect(key="etapas_a_rodar")
    assert etapas.options == list(caminhos.ETAPAS)


def test_operar_nao_oferece_rodar_durante_execucao(tmp_path, monkeypatch):
    """Duas abas abertas: a segunda ve a execucao da primeira e nao oferece o
    botao. Recusar depois do clique tambem funcionaria, mas um botao que
    sempre falha e pior que um botao desligado."""
    import os
    import time
    from streamlit.testing.v1 import AppTest
    pasta = tmp_path / ".execucao"
    pasta.mkdir()
    (pasta / "rodando.pid").write_text(f"{os.getpid()} {time.time():.0f} etapa1")
    monkeypatch.setattr(caminhos, "EXECUCAO", pasta)
    at = AppTest.from_file(str(caminhos.RAIZ / "painel.py"), default_timeout=30)
    at.run()
    assert not at.exception
    assert at.button(key="rodar").disabled


def test_parametro_carrega_a_nota_como_ajuda():
    """Cada numero aparece com a medicao que o escolheu ao lado. E o motivo de
    as notas existirem no config: um numero sem o motivo convida a ser mexido."""
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(str(caminhos.RAIZ / "painel.py"), default_timeout=30)
    at.run()
    ajuda = at.get("number_input")
    corte = [w for w in ajuda if w.key == "p_etapa1_corte"]
    assert corte and corte[0].help


def test_estado_nao_recalcula_o_md5_a_cada_chamada(input_limpo):
    """I6: `estado()` lia os 79 MB de `ligantes.tsv` inteiro a cada redesenho
    da aba, contra o que o proprio docstring promete. O md5 pede a mesma chave
    `(caminho, mtime)` que a leitura das tabelas ja usa."""
    p = input_limpo / "ligantes.tsv"
    p.write_text(LIGANTES, encoding="utf-8")
    entradas.estado()  # esquenta o cache
    lidos = {"bytes": 0}
    real = open

    def contando(arquivo, *a, **k):
        f = real(arquivo, *a, **k)
        if str(arquivo).endswith("ligantes.tsv"):
            leitura_original = f.read

            def ler(n=-1):
                dados = leitura_original(n)
                lidos["bytes"] += len(dados)
                return dados
            f.read = ler
        return f
    import builtins
    original = builtins.open
    builtins.open = contando
    try:
        entradas.estado()
    finally:
        builtins.open = original
    assert lidos["bytes"] == 0, f"releu {lidos['bytes']} bytes"


def test_estado_ve_o_arquivo_trocado(input_limpo):
    """O cache nao pode esconder um upload: a chave inclui o mtime."""
    import os
    import time
    p = input_limpo / "ligantes.tsv"
    p.write_text(LIGANTES, encoding="utf-8")
    md5_antes = entradas.estado()["ligantes"]["md5"]
    time.sleep(0.01)
    p.write_text(LIGANTES + "X\tCC\n", encoding="utf-8")
    os.utime(p, (time.time() + 1, time.time() + 1))
    assert entradas.estado()["ligantes"]["md5"] != md5_antes
