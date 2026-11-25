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
