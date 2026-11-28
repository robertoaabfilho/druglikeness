# -*- coding: utf-8 -*-
"""Onde fica cada coisa. Tres garantias: os caminhos saem da raiz do projeto,
importar o modulo nao mexe em disco, e nenhuma constante aponta para fora."""
import subprocess
import sys
from pathlib import Path


def test_caminhos_sao_relativos_a_raiz_do_projeto():
    from mods import caminhos
    assert caminhos.MODS.name == "mods"
    assert caminhos.TSV == caminhos.RAIZ / "tsv"
    assert caminhos.etapa("etapa1") == caminhos.RAIZ / "tsv" / "etapa1"
    assert caminhos.INPUT == caminhos.RAIZ / "input"
    assert caminhos.JSON == caminhos.RAIZ / "json"
    assert caminhos.DADOS == caminhos.RAIZ / "dados"
    assert caminhos.EXEMPLOS == caminhos.RAIZ / "exemplos"
    assert caminhos.MODELO == caminhos.RAIZ / "modelo"
    assert caminhos.CACHE == caminhos.RAIZ / "json" / "cache_rede"


def test_tsv_e_cache_rede_montam_subcaminhos():
    from mods import caminhos
    assert caminhos.tsv("etapa1", "evidencias.tsv") == (
        caminhos.RAIZ / "tsv" / "etapa1" / "evidencias.tsv")
    assert caminhos.cache_rede("ot-24.09") == caminhos.CACHE / "ot-24.09"


def test_importar_caminhos_nao_cria_pasta(tmp_path):
    """Sem efeito colateral: importar o modulo nao mexe em disco. Provado num
    processo limpo cuja RAIZ aponta para uma copia vazia -- se o import criasse
    pasta, ela apareceria aqui."""
    raiz = tmp_path / "projeto"
    (raiz / "mods").mkdir(parents=True)
    origem = Path(__file__).resolve().parent.parent / "mods" / "caminhos.py"
    (raiz / "mods" / "caminhos.py").write_text(
        origem.read_text(encoding="utf-8"), encoding="utf-8")
    r = subprocess.run(
        [sys.executable, "-c", "from mods import caminhos; print(caminhos.RAIZ)"],
        cwd=str(raiz), capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert sorted(p.name for p in raiz.iterdir()) == ["mods"]


def test_nenhum_caminho_aponta_para_fora_do_projeto():
    """Guarda do anti-objetivo: nada de ~/prospector, ~/druglikeness (o antigo)
    nem ~/modelos_baixados em codigo que roda."""
    from mods import caminhos
    fonte = Path(caminhos.__file__).read_text(encoding="utf-8")
    for proibido in ("prospector", "modelos_baixados"):
        assert proibido not in fonte, f"{proibido} citado em caminhos.py"
    assert "druglikeness/" not in fonte.replace("druglikeness2/", "")
