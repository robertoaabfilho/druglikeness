# -*- coding: utf-8 -*-
from mods.coletar_evidencias import main
from mods.evidencias.manifesto import ler_manifesto


def test_execucao_minima_produz_tabela_valida_e_manifesto(tmp_path):
    saida = tmp_path / "evidencias.tsv"
    codigo = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
                   "--alvos", "tests/dados/alvos_minimo.tsv",
                   "--sintomas", "tests/dados/sintomas.txt",
                   "--saida", str(saida), "--sem-rede"])
    assert codigo == 0
    assert saida.read_text(encoding="utf-8").splitlines()[0].startswith("#ligante\t")
    m = ler_manifesto(tmp_path / "manifesto.json")
    assert m["parametros"]["corte"] == 0.1
    assert "pipeline_commit" in m["versoes"]


def test_erro_de_entrada_interrompe_antes_de_escrever(tmp_path):
    saida = tmp_path / "evidencias.tsv"
    codigo = main(["--ligantes", "tests/dados/ligantes_sem_coluna.tsv",
                   "--alvos", "tests/dados/alvos_minimo.tsv",
                   "--sintomas", "tests/dados/sintomas.txt",
                   "--saida", str(saida), "--sem-rede"])
    assert codigo != 0
    assert not saida.exists()


def test_cria_diretorio_de_saida_ausente(tmp_path):
    saida = tmp_path / "novo" / "sub" / "evidencias.tsv"
    codigo = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
                   "--alvos", "tests/dados/alvos_minimo.tsv",
                   "--sintomas", "tests/dados/sintomas.txt",
                   "--saida", str(saida), "--sem-rede"])
    assert codigo == 0
    assert saida.exists() and (saida.parent / "manifesto.json").exists()


def test_execucao_como_script_direto(tmp_path):
    import subprocess
    import sys
    from pathlib import Path
    raiz = Path(__file__).resolve().parent.parent
    saida = tmp_path / "evidencias.tsv"
    r = subprocess.run(
        [sys.executable, "mods/coletar_evidencias.py",
         "--ligantes", "tests/dados/ligantes_minimo.tsv",
         "--alvos", "tests/dados/alvos_minimo.tsv",
         "--sintomas", "tests/dados/sintomas.txt",
         "--saida", str(saida), "--sem-rede"],
        cwd=raiz, capture_output=True, text=True, env={"PATH": "/usr/bin:/bin", "HOME": str(tmp_path)})
    assert r.returncode == 0, r.stderr
    assert saida.exists()
    assert "PARCIAL" in r.stderr
