# -*- coding: utf-8 -*-
"""Painel rodado pelo AppTest do Streamlit, sem navegador, sem upload real e
sem rede. Cobre so a casca: a logica exibida esta em test_resumo.py."""
from pathlib import Path

import pytest

from mods.evidencias.esquema import COLUNAS_EVIDENCIA

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest

RAIZ = Path(__file__).resolve().parent.parent
DADOS = Path(__file__).parent / "dados"


class _Envio:
    def __init__(self, caminho):
        self.name = Path(caminho).name
        self._caminho = caminho

    def getvalue(self):
        return Path(self._caminho).read_bytes()


@pytest.fixture
def ambiente(tmp_path, monkeypatch):
    rotulos = tmp_path / "rot.tsv"
    rotulos.write_text(
        "target_idx\tchembl_id\tpref_name\torganism\ttarget_type\tactive_molecules\n"
        "0\tCHEMBL222\tNET\tHomo sapiens\tSINGLE PROTEIN\t4000\n"
        "1\tCHEMBL228\tSERT\tHomo sapiens\tSINGLE PROTEIN\t20\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)  # a "pasta de trabalho" do painel e relativa
    for v in ("PP4CKD_PESOS", "EVIDENCIAS_CHEMBL"):
        monkeypatch.delenv(v, raising=False)
    monkeypatch.setenv("PP4CKD_ROTULOS", str(rotulos))
    return tmp_path


def _painel():
    return AppTest.from_file(str(RAIZ / "painel.py"), default_timeout=60)


def test_painel_sem_rotulos_pede_o_arquivo(ambiente, monkeypatch):
    monkeypatch.delenv("PP4CKD_ROTULOS")
    at = _painel().run()
    assert not at.exception
    assert any("target_labels.tsv" in e.value for e in at.error)


def test_painel_sem_pesos_avisa_modo_degradado(ambiente):
    at = _painel().run()
    assert not at.exception
    assert any("Modo degradado" in w.value for w in at.warning)


def test_painel_valida_os_envios_e_mostra_os_avisos(ambiente, monkeypatch):
    envios = {"Ligantes": "ligantes_minimo.tsv", "Alvos": "alvos_minimo.tsv",
              "Sintomas": "sintomas.txt"}

    def falso(self, label, *a, **k):
        return _Envio(DADOS / next(v for p, v in envios.items() if label.startswith(p)))
    monkeypatch.setattr("streamlit.delta_generator.DeltaGenerator.file_uploader", falso)
    at = _painel().run()
    assert not at.exception
    metricas = {m.label: m.value for m in at.metric}
    assert metricas["Ligantes válidos"] == "2"
    assert metricas["Alvos no modelo"] == "1 de 3"
    avisos = " ".join(x.value for x in list(at.warning) + list(at.error))
    assert "sem organismo" in avisos
    assert "não existem no modelo" in avisos
    assert "SMILES ausente ou ilegível" in avisos


def _com_envios(monkeypatch):
    envios = {"Ligantes": "ligantes_minimo.tsv", "Alvos": "alvos_minimo.tsv",
              "Sintomas": "sintomas.txt"}

    def falso(self, label, *a, **k):
        return _Envio(DADOS / next(v for p, v in envios.items() if label.startswith(p)))
    monkeypatch.setattr("streamlit.delta_generator.DeltaGenerator.file_uploader", falso)


def _preparar_lista(ambiente, efos, conteudo):
    """Grava uma lista enriquecida como se uma execucao anterior a tivesse feito,
    na pasta que o painel calcula para estes arquivos + estes EFO."""
    import hashlib
    import json
    from mods.evidencias.resumo import chave_enriquecimento
    chave = hashlib.md5(b"".join(
        (DADOS / n).read_bytes() for n in
        ("ligantes_minimo.tsv", "alvos_minimo.tsv", "sintomas.txt"))).hexdigest()
    enr = ambiente / "painel_saida" / "enriquecimento"
    pasta = enr / chave_enriquecimento(chave, efos, 0.4, 150, None)
    pasta.mkdir(parents=True)
    (pasta / "alvos_enriquecidos.tsv").write_text(conteudo, encoding="utf-8")
    (enr / f"ultimo_{chave}.json").write_text(json.dumps(efos), encoding="utf-8")
    return pasta


LISTA = "#alvo_chembl_id\talvo\torganismo\nCHEMBL222\tNET\tHomo sapiens\n"


def test_painel_mostra_a_lista_enriquecida_dos_mesmos_arquivos_e_efo(ambiente, monkeypatch):
    _com_envios(monkeypatch)
    _preparar_lista(ambiente, {"neuropathic pain": "EFO_1"}, LISTA)
    at = _painel().run()
    assert not at.exception
    assert any("Lista proposta" in s.value for s in at.subheader)
    assert any("EFO_1" in c.value for c in at.caption)


def test_painel_nao_mostra_lista_de_outros_efo_ou_arquivos(ambiente, monkeypatch):
    _com_envios(monkeypatch)
    pasta = _preparar_lista(ambiente, {"neuropathic pain": "EFO_1"}, LISTA)
    # o ponteiro diz EFO_1, mas a lista esta na pasta de outro EFO: nao aparece
    pasta.rename(pasta.parent / "outra_escolha")
    at = _painel().run()
    assert not at.exception
    assert not any("Lista proposta" in s.value for s in at.subheader)


def test_painel_aguenta_lista_enriquecida_vazia(ambiente, monkeypatch):
    _com_envios(monkeypatch)
    _preparar_lista(ambiente, {"neuropathic pain": "EFO_1"}, "")
    at = _painel().run()
    assert not at.exception
    assert any("não acrescentou nenhum alvo" in i.value for i in at.info)


def _chave_de_entrada():
    import hashlib
    return hashlib.md5(b"".join(
        (DADOS / n).read_bytes() for n in
        ("ligantes_minimo.tsv", "alvos_minimo.tsv", "sintomas.txt"))).hexdigest()


def _pasta_da_coleta(alvos=None):
    """A pasta que o painel calcula para a coleta destes arquivos, destes
    parametros e desta lista de alvos (os padroes da barra lateral)."""
    import hashlib
    from mods.evidencias.resumo import assinatura_de_arquivo, chave_coleta
    caminho = alvos or (Path("painel_saida") / "entrada" / "alvos" /
                        "alvos_minimo.tsv")
    assinatura = hashlib.md5(Path(caminho).read_bytes()).hexdigest()
    return Path("painel_saida") / "coleta" / chave_coleta(
        _chave_de_entrada(), assinatura, 0.1, None,
        assinatura_de_arquivo(""), assinatura_de_arquivo(""))


TABELA = "#" + "\t".join(COLUNAS_EVIDENCIA) + "\n"   # tabela valida, so cabecalho


def test_resultado_de_outra_execucao_nao_aparece_como_o_atual(ambiente, monkeypatch):
    """A tabela ficava num caminho fixo (coleta/evidencias.tsv) e era mostrada
    sob "Resultado" independentemente dos arquivos da tela."""
    _com_envios(monkeypatch)
    velha = ambiente / "painel_saida" / "coleta" / "evidencias.tsv"
    velha.parent.mkdir(parents=True)
    velha.write_text(TABELA, encoding="utf-8")
    at = _painel().run()
    assert not at.exception
    assert not any("Resultado" in h.value for h in at.header)


def test_resultado_destas_entradas_aparece_com_a_procedencia(ambiente, monkeypatch):
    import json
    _com_envios(monkeypatch)
    # a coleta roda sobre a copia que o painel salva do envio
    pasta = ambiente / _pasta_da_coleta(
        DADOS / "alvos_minimo.tsv")  # mesmo conteudo da copia salva
    pasta.mkdir(parents=True)
    (pasta / "evidencias.tsv").write_text(TABELA, encoding="utf-8")
    (pasta / "manifesto.json").write_text(json.dumps(
        {"versoes": {"data": "2026-11-25T00:00:00+00:00"},
         "contagens": {"linhas_evidencia": 0, "execucao_parcial": True}}),
        encoding="utf-8")
    at = _painel().run()
    assert not at.exception
    assert any("Resultado" in h.value for h in at.header)
    assert any("2026-11-25" in c.value for c in at.caption)
    assert any("PARCIAL" in w.value for w in at.warning)


def test_tabela_sem_manifesto_e_marcada(ambiente, monkeypatch):
    _com_envios(monkeypatch)
    pasta = ambiente / _pasta_da_coleta(DADOS / "alvos_minimo.tsv")
    pasta.mkdir(parents=True)
    (pasta / "evidencias.tsv").write_text(TABELA, encoding="utf-8")
    at = _painel().run()
    assert not at.exception
    assert any("não é citável" in w.value for w in at.warning)


def test_descartar_edicoes_cria_um_editor_novo_e_sem_pendencia(ambiente, monkeypatch):
    """Nao se verificou que remover o estado do widget do session_state
    descarta as edicoes (o Streamlit o reenvia no rerun). O descarte passa a
    ser uma GERACAO nova do editor: outro widget, estado novo."""
    from mods.evidencias.resumo import edicao_pendente
    _com_envios(monkeypatch)
    _preparar_lista(ambiente, {"neuropathic pain": "EFO_1"}, LISTA)
    at = _painel().run()
    assert not at.exception
    chave = at.session_state["editor_atual"]
    assert chave and chave.endswith("_g0")
    at.session_state[chave] = {"edited_rows": {0: {"alvo": "editado a mao"}},
                               "added_rows": [], "deleted_rows": []}
    at = at.run()
    assert at.session_state["editor_atual"] == chave
    assert edicao_pendente(at.session_state[chave])   # pendencia reconhecida
    at.session_state["geracao_editor"] = 1            # o que o botao faz
    at = at.run()
    nova = at.session_state["editor_atual"]
    assert nova != chave and nova.endswith("_g1")
    assert not edicao_pendente(at.session_state.get(nova))


class _ProcessoTerminado:
    """Trabalho que acabou, para exercitar o portao de redesenho sem subir
    processo nenhum."""
    pid = -1

    def poll(self):
        return 0


def _trabalho(ambiente, pasta):
    import time
    log = ambiente / "coletar.log"
    log.write_text("$ cli\n0 linhas de evidencia\n", encoding="utf-8")
    return {"tipo": "coletar", "proc": _ProcessoTerminado(), "log": log,
            "inicio": time.time(), "saida": pasta / "evidencias.tsv"}


def _edicao(at, chave):
    at.session_state[chave] = {"edited_rows": {0: {"alvo": "editado a mao"}},
                               "added_rows": [], "deleted_rows": []}


def test_aviso_de_edicao_so_aparece_se_um_redesenho_foi_suprimido(ambiente, monkeypatch):
    """O aviso e o botao de descarte (destrutivo) ficavam fora do portao e
    apareciam a cada poll depois do fim do trabalho -- inclusive quando a
    edicao comecou DEPOIS do redesenho, sem nada ter sido suprimido."""
    _com_envios(monkeypatch)
    _preparar_lista(ambiente, {"neuropathic pain": "EFO_1"}, LISTA)
    at = _painel().run()
    chave = at.session_state["editor_atual"]
    # trabalho termina SEM edicao pendente: a pagina e redesenhada, nada suprimido
    at.session_state["trabalho"] = _trabalho(ambiente, ambiente)
    at = at.run()
    assert not at.exception
    assert at.session_state["trabalho"]["concluido"] == 0
    assert "rerun_suprimido" not in at.session_state["trabalho"]
    # o pesquisador comeca a editar DEPOIS: nenhum aviso, nenhum botao de descarte
    _edicao(at, chave)
    at = at.run()
    avisos = " ".join(w.value for w in at.warning)
    assert "edições não aplicadas" not in avisos
    assert not any("descarta" in b.label for b in at.button)


def test_aviso_de_edicao_aparece_quando_a_edicao_estava_em_curso(ambiente, monkeypatch):
    _com_envios(monkeypatch)
    _preparar_lista(ambiente, {"neuropathic pain": "EFO_1"}, LISTA)
    at = _painel().run()
    _edicao(at, at.session_state["editor_atual"])
    at.session_state["trabalho"] = _trabalho(ambiente, ambiente)
    at = at.run()
    assert not at.exception
    assert at.session_state["trabalho"]["rerun_suprimido"] is True
    assert any("edições não aplicadas" in w.value for w in at.warning)
    assert any("descarta" in b.label for b in at.button)
