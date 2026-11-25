# -*- coding: utf-8 -*-
"""Painel rodado pelo AppTest do Streamlit, sem navegador, sem upload real e
sem rede. Cobre so a casca: a logica exibida esta em test_resumo.py."""
from pathlib import Path

import pytest

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
