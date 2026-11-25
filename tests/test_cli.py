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


ARGS = ["--ligantes", "tests/dados/ligantes_minimo.tsv",
        "--alvos", "tests/dados/alvos_minimo.tsv",
        "--sintomas", "tests/dados/sintomas.txt", "--sem-rede"]


def test_commit_e_hash_de_40_hex_em_execucao_normal(tmp_path):
    import re
    main(ARGS + ["--saida", str(tmp_path / "e.tsv")])
    c = ler_manifesto(tmp_path / "manifesto.json")["versoes"]["pipeline_commit"]
    assert re.fullmatch(r"[0-9a-f]{40}", c)


def test_sem_git_o_commit_e_desconhecido(tmp_path, monkeypatch):
    import subprocess

    def sem_git(*a, **k):
        raise FileNotFoundError("git")
    monkeypatch.setattr(subprocess, "run", sem_git)
    assert main(ARGS + ["--saida", str(tmp_path / "e.tsv")]) == 0
    m = ler_manifesto(tmp_path / "manifesto.json")
    assert m["versoes"]["pipeline_commit"] == "desconhecido"


def test_sintomas_vazio_e_aviso_e_a_execucao_continua(tmp_path, capsys):
    vazio = tmp_path / "sintomas.txt"
    vazio.write_text("# so comentario\n", encoding="utf-8")
    args = [a if a != "tests/dados/sintomas.txt" else str(vazio) for a in ARGS]
    assert main(args + ["--saida", str(tmp_path / "e.tsv")]) == 0
    err = capsys.readouterr().err
    assert "[aviso] sintomas" in err and "NENHUM enriquecimento" in err


def test_limite_invalido_e_rejeitado(tmp_path):
    import pytest
    for ruim in ("0", "-3", "x"):
        with pytest.raises(SystemExit):
            main(ARGS + ["--saida", str(tmp_path / "e.tsv"), "--limite", ruim])
    assert not (tmp_path / "e.tsv").exists()


def test_limite_corta_durante_a_leitura(tmp_path):
    # aspirina e o 1o ligante valido; a linha ruim (sem_estrutura) vem depois
    # do limite e nao pode gerar problema nem entrar na contagem.
    main(ARGS + ["--saida", str(tmp_path / "e.tsv"), "--limite", "1"])
    c = ler_manifesto(tmp_path / "manifesto.json")["contagens"]
    assert c["ligantes"] == 1
    assert c["problemas_ligantes"] == 0


def test_limite_nao_le_linha_malformada_apos_o_ultimo_aceito(tmp_path):
    from mods.evidencias.entrada import ler_ligantes
    arq = tmp_path / "lig.tsv"
    arq.write_text(
        "nome\tsmiles\n"
        "aspirina\tCC(=O)Oc1ccccc1C(=O)O\n"
        "linha_larga\tCCO\textra\tmais\n"      # largura errada, logo apos o 1o
        "etanol\tCCO\n", encoding="utf-8")
    ligantes, problemas = ler_ligantes(arq, limite=1)
    assert [l["nome"] for l in ligantes] == ["aspirina"]
    assert problemas == []
    # sem limite a linha malformada e relatada (prova que o caso e alcancavel)
    _, problemas = ler_ligantes(arq)
    assert any(p.campo == "linha" for p in problemas)


def test_limite_ainda_relata_duplicados_vistos_antes_do_corte(tmp_path):
    from mods.evidencias.entrada import ler_ligantes
    arq = tmp_path / "lig.tsv"
    arq.write_text(
        "nome\tsmiles\n"
        "duloxetina base\tCNCCC(Oc1cccc2ccccc12)c1cccs1\n"
        "cloridrato\tCNCCC(Oc1cccc2ccccc12)c1cccs1.Cl\n"
        "aspirina\tCC(=O)Oc1ccccc1C(=O)O\n", encoding="utf-8")
    ligantes, problemas = ler_ligantes(arq, limite=2)
    assert [l["nome"] for l in ligantes] == ["duloxetina base", "aspirina"]
    assert any("duplicado" in p.mensagem for p in problemas)


def test_leitura_fecha_o_arquivo_quando_falta_coluna(tmp_path, recwarn):
    import gc
    from mods.evidencias.entrada import ler_ligantes
    ler_ligantes("tests/dados/ligantes_sem_coluna.tsv")
    gc.collect()
    assert not [w for w in recwarn if issubclass(w.category, ResourceWarning)]


def test_checkpoint_sem_rotulos_falha_com_mensagem(tmp_path, capsys):
    from mods.coletar_evidencias import main
    rc = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
               "--alvos", "tests/dados/alvos_minimo.tsv",
               "--sintomas", "tests/dados/sintomas.txt",
               "--saida", str(tmp_path / "e.tsv"), "--checkpoint", "x.pt"])
    assert rc == 2
    assert "--checkpoint e --rotulos" in capsys.readouterr().err


def _alvos_do_modelo_falso(tmp_path):
    alvos = tmp_path / "alvos.tsv"
    alvos.write_text("alvo_chembl_id\nCHEMBL100\nCHEMBL101\n", encoding="utf-8")
    return alvos


def _rotulos_tsv(tmp_path, rotulos):
    caminho = tmp_path / "rotulos.tsv"
    linhas = ["target_idx\tchembl_id\tpref_name\torganism\ttarget_type\tactive_molecules"]
    for i, r in enumerate(rotulos):
        linhas.append(f"{i}\t{r['chembl_id']}\t{r['nome']}\t{r['organismo']}"
                      f"\t{r['tipo']}\t{r['n_ativos']}")
    caminho.write_text("\n".join(linhas) + "\n", encoding="utf-8")
    return caminho


def _rodar_com_modelo(tmp_path, modelo, rotulos, corte, nome):
    from mods.evidencias.esquema import ler_evidencias
    saida = tmp_path / nome / "evidencias.tsv"
    rot = _rotulos_tsv(tmp_path, rotulos)
    codigo = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
                   "--alvos", str(_alvos_do_modelo_falso(tmp_path)),
                   "--sintomas", "tests/dados/sintomas.txt", "--saida", str(saida),
                   "--checkpoint", modelo, "--rotulos", str(rot),
                   "--corte", str(corte)])
    assert codigo == 0
    return list(ler_evidencias(saida)), ler_manifesto(saida.parent / "manifesto.json"), rot


def test_checkpoint_grava_linhas_e_manifesto_confere(tmp_path, modelo_falso, rotulos_falsos):
    import hashlib
    todas, m, rot = _rodar_com_modelo(tmp_path, modelo_falso, rotulos_falsos, 0.0, "a")
    assert todas and all(l["tipo_evidencia"] == "predito" for l in todas)
    assert m["contagens"]["linhas_evidencia"] == len(todas)
    assert m["parametros"]["rotulos"] == str(rot)
    assert m["parametros"]["rotulos_md5"] == hashlib.md5(rot.read_bytes()).hexdigest()
    assert m["contagens"]["smiles_ilegiveis_modelo"] == 0
    filtradas, m2, _ = _rodar_com_modelo(tmp_path, modelo_falso, rotulos_falsos, 0.505, "b")
    assert 0 < len(filtradas) < len(todas)
    assert all(float(l["valor"]) >= 0.505 for l in filtradas)
    assert m2["contagens"]["linhas_evidencia"] == len(filtradas)


def test_checkpoint_invalido_nao_deixa_tabela_nem_manifesto(tmp_path, rotulos_falsos, capsys):
    ruim = tmp_path / "ruim.pt"
    ruim.write_bytes(b"isto nao e um checkpoint")
    saida = tmp_path / "o" / "evidencias.tsv"
    codigo = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
                   "--alvos", str(_alvos_do_modelo_falso(tmp_path)),
                   "--sintomas", "tests/dados/sintomas.txt", "--saida", str(saida),
                   "--checkpoint", str(ruim),
                   "--rotulos", str(_rotulos_tsv(tmp_path, rotulos_falsos))])
    assert codigo == 2
    assert not saida.exists()
    assert not (saida.parent / "manifesto.json").exists()
    assert "[erro]" in capsys.readouterr().err


def test_sem_torch_falha_limpo_e_sem_arquivo(tmp_path, rotulos_falsos, monkeypatch, capsys):
    import importlib.util
    monkeypatch.setattr(importlib.util, "find_spec", lambda nome: None)
    saida = tmp_path / "o" / "evidencias.tsv"
    codigo = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
                   "--alvos", str(_alvos_do_modelo_falso(tmp_path)),
                   "--sintomas", "tests/dados/sintomas.txt", "--saida", str(saida),
                   "--checkpoint", "x.pt",
                   "--rotulos", str(_rotulos_tsv(tmp_path, rotulos_falsos))])
    assert codigo == 2 and not saida.exists()
    assert "pip install torch rdkit" in capsys.readouterr().err


def test_falha_preserva_saida_preexistente(tmp_path, rotulos_falsos):
    saida = tmp_path / "evidencias.tsv"
    saida.write_bytes(b"#antigo\nconteudo precioso\n")
    codigo = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
                   "--alvos", str(_alvos_do_modelo_falso(tmp_path)),
                   "--sintomas", "tests/dados/sintomas.txt", "--saida", str(saida),
                   "--checkpoint", str(tmp_path / "nao_existe.pt"),
                   "--rotulos", str(_rotulos_tsv(tmp_path, rotulos_falsos))])
    assert codigo == 2
    assert saida.read_bytes() == b"#antigo\nconteudo precioso\n"
    assert not (tmp_path / "manifesto.json").exists()
    assert [p.name for p in tmp_path.glob("*.tmp")] == []


def test_sucesso_substitui_saida_preexistente(tmp_path, modelo_falso, rotulos_falsos):
    from mods.evidencias.esquema import ler_evidencias
    saida = tmp_path / "evidencias.tsv"
    saida.write_bytes(b"#antigo\nlixo\n")
    codigo = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
                   "--alvos", str(_alvos_do_modelo_falso(tmp_path)),
                   "--sintomas", "tests/dados/sintomas.txt", "--saida", str(saida),
                   "--checkpoint", modelo_falso, "--corte", "0",
                   "--rotulos", str(_rotulos_tsv(tmp_path, rotulos_falsos))])
    assert codigo == 0
    assert b"lixo" not in saida.read_bytes()
    assert list(ler_evidencias(saida))
    assert [p.name for p in tmp_path.glob("*.tmp")] == []


def test_rotulos_inexistente_sem_checkpoint_falha_limpo(tmp_path, capsys):
    saida = tmp_path / "evidencias.tsv"
    codigo = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
                   "--alvos", "tests/dados/alvos_minimo.tsv",
                   "--sintomas", "tests/dados/sintomas.txt", "--saida", str(saida),
                   "--rotulos", str(tmp_path / "nao_existe.tsv")])
    assert codigo == 2 and not saida.exists()
    assert "[erro]" in capsys.readouterr().err


def test_falha_com_checkpoint_corrompido_preserva_saida(tmp_path, rotulos_falsos):
    ruim = tmp_path / "ruim.pt"
    ruim.write_bytes(b"\x00garbage\xff" * 50)
    saida = tmp_path / "evidencias.tsv"
    saida.write_bytes(b"#antigo\nprecioso\n")
    codigo = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
                   "--alvos", str(_alvos_do_modelo_falso(tmp_path)),
                   "--sintomas", "tests/dados/sintomas.txt", "--saida", str(saida),
                   "--checkpoint", str(ruim),
                   "--rotulos", str(_rotulos_tsv(tmp_path, rotulos_falsos))])
    assert codigo == 2
    assert saida.read_bytes() == b"#antigo\nprecioso\n"
    assert list(tmp_path.glob("*.tmp")) == []


def test_modo_da_tabela_igual_ao_do_manifesto(tmp_path, modelo_falso, rotulos_falsos):
    import stat
    saida = tmp_path / "o" / "evidencias.tsv"
    codigo = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
                   "--alvos", str(_alvos_do_modelo_falso(tmp_path)),
                   "--sintomas", "tests/dados/sintomas.txt", "--saida", str(saida),
                   "--checkpoint", modelo_falso,
                   "--rotulos", str(_rotulos_tsv(tmp_path, rotulos_falsos))])
    assert codigo == 0
    modo = lambda p: stat.S_IMODE(p.stat().st_mode)
    assert modo(saida) == modo(saida.parent / "manifesto.json")


def test_diretorio_de_saida_sem_permissao_falha_limpo(tmp_path, capsys):
    import os
    import pytest
    if os.geteuid() == 0:
        pytest.skip("root ignora permissoes de diretorio")
    destino = tmp_path / "ro"
    destino.mkdir()
    destino.chmod(0o500)
    try:
        codigo = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
                       "--alvos", "tests/dados/alvos_minimo.tsv",
                       "--sintomas", "tests/dados/sintomas.txt",
                       "--saida", str(destino / "evidencias.tsv"), "--sem-rede"])
    finally:
        destino.chmod(0o700)
    assert codigo == 2
    err = capsys.readouterr().err
    assert "[erro]" in err and "Traceback" not in err


def test_chmod_que_falha_nao_custa_a_tabela(tmp_path, modelo_falso, rotulos_falsos, monkeypatch):
    import os
    def recusa(*a, **k):
        raise OSError("filesystem sem chmod")
    monkeypatch.setattr(os, "chmod", recusa)
    saida = tmp_path / "evidencias.tsv"
    codigo = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
                   "--alvos", str(_alvos_do_modelo_falso(tmp_path)),
                   "--sintomas", "tests/dados/sintomas.txt", "--saida", str(saida),
                   "--checkpoint", modelo_falso, "--corte", "0",
                   "--rotulos", str(_rotulos_tsv(tmp_path, rotulos_falsos))])
    assert codigo == 0
    assert len(saida.read_text(encoding="utf-8").splitlines()) > 1


def test_enriquecimento_avisa_que_a_evidencia_nao_o_usou(tmp_path, monkeypatch, capsys):
    import json
    import mods.coletar_evidencias as cli
    rotulos = tmp_path / "labels.tsv"
    rotulos.write_text("target_idx\ttid\tchembl_id\tpref_name\ttarget_type\torganism\t"
                       "active_molecules\tis_ckd\n0\t1\tCHEMBL222\tNET\tSINGLE PROTEIN\t"
                       "Homo sapiens\t5\tTrue\n", encoding="utf-8")
    monkeypatch.setattr(cli, "_enriquecer", lambda *a: [
        {"alvo_chembl_id": "CHEMBL222", "origem": "usuario"},
        {"alvo_chembl_id": "CHEMBL233", "origem": "open_targets"}])
    enr = tmp_path / "enr.tsv"
    codigo = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
                   "--alvos", "tests/dados/alvos_minimo.tsv",
                   "--sintomas", "tests/dados/sintomas.txt", "--saida", str(tmp_path / "e.tsv"),
                   "--rotulos", str(rotulos), "--alvos-enriquecidos", str(enr)])
    assert codigo == 0
    erro = capsys.readouterr().err
    assert "2 alvos gravados" in erro and "NAO os usou" in erro and f"--alvos {enr}" in erro
    c = json.loads((tmp_path / "manifesto.json").read_text())["contagens"]
    assert c["alvos_enriquecidos_escritos"] == 2
    assert c["alvos_enriquecidos_usados_na_evidencia"] is False
    # Rodar a coleta com o proprio arquivo enriquecido muda o flag e some o aviso.
    saida2 = tmp_path / "e2.tsv"
    codigo = main(["--ligantes", "tests/dados/ligantes_minimo.tsv", "--alvos", str(enr),
                   "--sintomas", "tests/dados/sintomas.txt", "--saida", str(saida2),
                   "--rotulos", str(rotulos), "--alvos-enriquecidos", str(enr)])
    assert codigo == 0
    c = json.loads((tmp_path / "manifesto.json").read_text())["contagens"]
    assert c["alvos_enriquecidos_usados_na_evidencia"] is True
