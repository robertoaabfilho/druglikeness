# -*- coding: utf-8 -*-
"""Escrita atomica das saidas: uma escrita que falha no meio nao pode destruir
o resultado anterior do usuario. Vale para as TRES saidas -- a tabela, a lista
de alvos enriquecidos (que o README manda revisar a mao) e o manifesto (que e
a chave da repeticao)."""
import pathlib

import pytest

from mods.evidencias import atomico
from mods.evidencias.manifesto import escrever_manifesto, ler_manifesto
from mods.evidencias.open_targets import escrever_alvos_enriquecidos


@pytest.fixture
def escrita_que_falha(monkeypatch):
    """Faz a escrita do temporario falhar, sem tocar nos demais arquivos."""
    real = pathlib.Path.write_text

    def falha(self, *a, **k):
        if self.suffix == ".tmp":
            raise OSError("disco cheio")
        return real(self, *a, **k)
    monkeypatch.setattr(pathlib.Path, "write_text", falha)


def test_escrever_texto_preserva_o_anterior(tmp_path, escrita_que_falha):
    destino = tmp_path / "x.tsv"
    destino.write_text("anterior\n", encoding="utf-8")
    with pytest.raises(OSError):
        atomico.escrever_texto(destino, "novo\n")
    assert destino.read_text(encoding="utf-8") == "anterior\n"
    assert list(tmp_path.glob("*.tmp")) == []


def test_manifesto_anterior_sobrevive_a_falha(tmp_path, escrita_que_falha):
    destino = tmp_path / "manifesto.json"
    destino.write_text('{"parametros": {"corte": 0.1}}', encoding="utf-8")
    with pytest.raises(OSError):
        escrever_manifesto(destino, parametros={}, versoes={}, contagens={})
    assert ler_manifesto(destino)["parametros"]["corte"] == 0.1
    assert list(tmp_path.glob("*.tmp")) == []


def test_lista_enriquecida_revisada_sobrevive_a_falha(tmp_path, escrita_que_falha):
    destino = tmp_path / "alvos_enriquecidos.tsv"
    destino.write_text("#alvo_chembl_id\nCHEMBL222\n", encoding="utf-8")
    with pytest.raises(OSError):
        escrever_alvos_enriquecidos(destino, [{"alvo_chembl_id": "CHEMBL233"}])
    assert destino.read_text(encoding="utf-8") == "#alvo_chembl_id\nCHEMBL222\n"
    assert list(tmp_path.glob("*.tmp")) == []


def test_escrita_bem_sucedida_substitui_e_nao_deixa_temporario(tmp_path):
    destino = tmp_path / "alvos_enriquecidos.tsv"
    destino.write_text("velho\n", encoding="utf-8")
    escrever_alvos_enriquecidos(destino, [{"alvo_chembl_id": "CHEMBL233"}])
    texto = destino.read_text(encoding="utf-8")
    assert texto.startswith("#alvo_chembl_id\t") and "CHEMBL233" in texto
    assert list(tmp_path.glob("*.tmp")) == []
