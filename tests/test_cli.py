# -*- coding: utf-8 -*-
from mods.coletar_evidencias import SAIDA_PARCIAL, main
from mods.evidencias.manifesto import ler_manifesto

# As fixtures tests/dados/*_minimo.* tem, de proposito, uma linha ruim cada:
# toda execucao com elas e PARCIAL, e e assim que se prova que o codigo de
# saida distingue "completa" de "incompleta" (automacao nao le o stderr).


def test_execucao_minima_produz_tabela_valida_e_manifesto(tmp_path):
    saida = tmp_path / "evidencias.tsv"
    codigo = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
                   "--alvos", "tests/dados/alvos_minimo.tsv",
                   "--sintomas", "tests/dados/sintomas.txt",
                   "--saida", str(saida), "--sem-rede"])
    assert codigo == SAIDA_PARCIAL
    assert saida.read_text(encoding="utf-8").splitlines()[0].startswith("#ligante\t")
    m = ler_manifesto(tmp_path / "manifesto.json")
    assert m["parametros"]["corte"] == 0.1
    assert "pipeline_commit" in m["versoes"]
    assert m["versoes"]["python"] and m["versoes"]["rdkit"] not in ("", "ausente")
    assert m["contagens"]["execucao_parcial"] is True


def test_entrada_limpa_sai_com_zero_e_manifesto_diz_completa(tmp_path):
    lig = tmp_path / "lig.tsv"
    lig.write_text("nome\tsmiles\naspirina\tCC(=O)Oc1ccccc1C(=O)O\n", encoding="utf-8")
    alv = tmp_path / "alvos.tsv"
    alv.write_text("alvo_chembl_id\nCHEMBL222\n", encoding="utf-8")
    codigo = main(["--ligantes", str(lig), "--alvos", str(alv), "--sintomas",
                   "tests/dados/sintomas.txt", "--saida", str(tmp_path / "e.tsv"),
                   "--sem-rede"])
    assert codigo == 0
    assert ler_manifesto(tmp_path / "manifesto.json")[
        "contagens"]["execucao_parcial"] is False


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
    assert codigo == SAIDA_PARCIAL
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
    # Codigo PROPRIO para execucao parcial: antes era 0, e automacao nao tinha
    # como distinguir uma coleta completa de uma que descartou linhas de entrada.
    assert r.returncode == SAIDA_PARCIAL, r.stderr
    assert saida.exists()
    assert "PARCIAL" in r.stderr


ARGS = ["--ligantes", "tests/dados/ligantes_minimo.tsv",
        "--alvos", "tests/dados/alvos_minimo.tsv",
        "--sintomas", "tests/dados/sintomas.txt", "--sem-rede"]


def test_commit_e_hash_de_40_hex_em_execucao_normal(tmp_path):
    """Num checkout de git, o manifesto carrega o commit.

    Fora de um checkout -- e o PACOTE para outra maquina e exatamente isso --
    nao ha commit a registrar, e `desconhecido` e a resposta honesta. Exigir o
    hash ali faria a suite do pacote falhar, e e a suite verde que prova que o
    pacote chegou inteiro."""
    import re
    import subprocess
    main(ARGS + ["--saida", str(tmp_path / "e.tsv")])
    versoes = ler_manifesto(tmp_path / "manifesto.json")["versoes"]
    em_git = subprocess.run(["git", "rev-parse", "--git-dir"],
                            capture_output=True).returncode == 0
    if not em_git:
        assert versoes["pipeline_commit"] == "desconhecido"
        return
    # Arvore suja ganha sufixo: hash limpo de codigo modificado seria uma
    # afirmacao falsa de reprodutibilidade.
    assert re.fullmatch(r"[0-9a-f]{40}(\+sujo)?", versoes["pipeline_commit"])
    assert versoes["pipeline_repositorio"]


def test_arvore_suja_e_marcada_e_o_repositorio_de_fora_fica_visivel(monkeypatch):
    from pathlib import Path
    # A procedencia e compartilhada entre as etapas (mods/evidencias/procedencia.py).
    import mods.evidencias.procedencia as cli
    respostas = {}
    monkeypatch.setattr(cli, "_git", lambda args, cwd: respostas.get(args[0]))
    respostas = {"rev-parse": "a" * 40, "status": ""}
    assert cli.commit_do_pipeline() == "a" * 40
    respostas = {"rev-parse": "a" * 40, "status": " M mods/evidencias/entrada.py"}
    assert cli.commit_do_pipeline() == "a" * 40 + "+sujo"
    respostas = {"rev-parse": "a" * 40, "status": None}
    assert cli.commit_do_pipeline() == "a" * 40 + "+estado_desconhecido"
    # Pacote embutido: o commit e de outro repositorio, e isso tem de aparecer.
    respostas = {"rev-parse": str(Path("/outro/repo"))}
    dados = cli.repositorio_do_pipeline()
    assert dados["pipeline_repositorio"] == "/outro/repo"
    assert "embutido em outro repositorio" in dados["pipeline_repositorio_aviso"]


def test_sem_git_o_commit_e_desconhecido(tmp_path, monkeypatch):
    import subprocess

    def sem_git(*a, **k):
        raise FileNotFoundError("git")
    monkeypatch.setattr(subprocess, "run", sem_git)
    assert main(ARGS + ["--saida", str(tmp_path / "e.tsv")]) == SAIDA_PARCIAL
    m = ler_manifesto(tmp_path / "manifesto.json")
    assert m["versoes"]["pipeline_commit"] == "desconhecido"


def test_sintomas_vazio_e_aviso_e_a_execucao_continua(tmp_path, capsys):
    vazio = tmp_path / "sintomas.txt"
    vazio.write_text("# so comentario\n", encoding="utf-8")
    args = [a if a != "tests/dados/sintomas.txt" else str(vazio) for a in ARGS]
    assert main(args + ["--saida", str(tmp_path / "e.tsv")]) == SAIDA_PARCIAL
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
    assert codigo == SAIDA_PARCIAL
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
    assert codigo == SAIDA_PARCIAL
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
    assert codigo == SAIDA_PARCIAL
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
    assert codigo == SAIDA_PARCIAL
    assert len(saida.read_text(encoding="utf-8").splitlines()) > 1


def test_enriquecimento_avisa_que_a_evidencia_nao_o_usou(tmp_path, monkeypatch, capsys):
    import json
    import mods.coletar_evidencias as cli
    rotulos = tmp_path / "labels.tsv"
    rotulos.write_text("target_idx\ttid\tchembl_id\tpref_name\ttarget_type\torganism\t"
                       "active_molecules\tis_ckd\n0\t1\tCHEMBL222\tNET\tSINGLE PROTEIN\t"
                       "Homo sapiens\t5\tTrue\n", encoding="utf-8")
    monkeypatch.setattr(cli, "_abrir_cache", lambda *a: (None, {}, None, None))
    monkeypatch.setattr(cli, "_enriquecer", lambda *a: [
        {"alvo_chembl_id": "CHEMBL222", "origem": "usuario"},
        {"alvo_chembl_id": "CHEMBL233", "origem": "open_targets"}])
    enr = tmp_path / "enr.tsv"
    codigo = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
                   "--alvos", "tests/dados/alvos_minimo.tsv",
                   "--sintomas", "tests/dados/sintomas.txt", "--saida", str(tmp_path / "e.tsv"),
                   "--rotulos", str(rotulos), "--alvos-enriquecidos", str(enr)])
    assert codigo == SAIDA_PARCIAL
    erro = capsys.readouterr().err
    assert "2 alvos gravados" in erro and "NAO os usou" in erro and f"--alvos {enr}" in erro
    c = json.loads((tmp_path / "manifesto.json").read_text())["contagens"]
    assert c["alvos_enriquecidos_escritos"] == 2
    assert c["alvos_enriquecidos_usados_na_evidencia"] is False
    # Passo 2 do README: coletar COM a lista enriquecida, sem --alvos-enriquecidos.
    # O flag vem da procedencia dos alvos lidos, nao da flag da linha de comando:
    # antes, este passo registrava 'false', o contrario da verdade.
    saida2 = tmp_path / "dois" / "e2.tsv"
    codigo = main(["--ligantes", "tests/dados/ligantes_minimo.tsv", "--alvos", str(enr),
                   "--sintomas", "tests/dados/sintomas.txt", "--saida", str(saida2),
                   "--rotulos", str(rotulos)])
    assert codigo == SAIDA_PARCIAL
    c = json.loads((saida2.parent / "manifesto.json").read_text())["contagens"]
    assert c["alvos_enriquecidos_usados_na_evidencia"] is True
    assert c["origem_dos_alvos_da_evidencia"] == {"usuario": 1, "open_targets": 1}


def test_alvos_e_alvos_enriquecidos_no_mesmo_arquivo_e_recusado(tmp_path, capsys):
    lista = tmp_path / "alvos.tsv"
    lista.write_text("alvo_chembl_id\nCHEMBL222\n", encoding="utf-8")
    antes = lista.read_bytes()
    rc = main(["--ligantes", "tests/dados/ligantes_minimo.tsv", "--alvos", str(lista),
               "--sintomas", "tests/dados/sintomas.txt", "--saida", str(tmp_path / "e.tsv"),
               "--rotulos", "tests/dados/alvos_minimo.tsv",
               "--alvos-enriquecidos", str(lista)])
    assert rc == 2
    assert lista.read_bytes() == antes   # a revisao do usuario fica de pe
    assert "mesmo arquivo" in capsys.readouterr().err


def test_saida_e_alvos_enriquecidos_no_mesmo_arquivo_e_recusado(tmp_path, capsys):
    saida = tmp_path / "e.tsv"
    rc = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
               "--alvos", "tests/dados/alvos_minimo.tsv",
               "--sintomas", "tests/dados/sintomas.txt", "--saida", str(saida),
               "--rotulos", "tests/dados/alvos_minimo.tsv",
               "--alvos-enriquecidos", str(saida)])
    assert rc == 2 and not saida.exists()
    assert "uma apagaria a outra" in capsys.readouterr().err


def _args_minimos(tmp_path, *extra):
    return ["--ligantes", "tests/dados/ligantes_minimo.tsv",
            "--alvos", "tests/dados/alvos_minimo.tsv",
            "--sintomas", "tests/dados/sintomas.txt",
            "--saida", str(tmp_path / "evidencias.tsv"), "--sem-rede", *extra]


def test_falha_ao_preparar_o_manifesto_nao_publica_a_tabela(tmp_path, monkeypatch, capsys):
    """Tabela sem manifesto nao e citavel: era o que sobrava quando a escrita
    do manifesto falhava DEPOIS da troca da tabela."""
    from mods.evidencias import atomico
    real = atomico.temporario

    def falha(destino):
        if str(destino).endswith("manifesto.json"):
            raise OSError("disco cheio")
        return real(destino)
    monkeypatch.setattr(atomico, "temporario", falha)
    saida = tmp_path / "evidencias.tsv"
    saida.write_bytes(b"#antigo\nprecioso\n")
    (tmp_path / "manifesto.json").write_text('{"antigo": true}', encoding="utf-8")
    assert main(_args_minimos(tmp_path)) == 2
    assert saida.read_bytes() == b"#antigo\nprecioso\n"
    assert (tmp_path / "manifesto.json").read_text(encoding="utf-8") == '{"antigo": true}'
    assert list(tmp_path.glob("*.tmp")) == []
    assert "Traceback" not in capsys.readouterr().err


def test_falha_na_lista_enriquecida_nao_publica_a_tabela(tmp_path, monkeypatch, capsys):
    import mods.coletar_evidencias as cli
    from mods.evidencias import atomico
    rotulos = tmp_path / "labels.tsv"
    rotulos.write_text("target_idx\tchembl_id\tpref_name\torganism\ttarget_type\t"
                       "active_molecules\n0\tCHEMBL222\tNET\tHomo sapiens\t"
                       "SINGLE PROTEIN\t5\n", encoding="utf-8")
    monkeypatch.setattr(cli, "_abrir_cache", lambda *a: (None, {}, None, None))
    monkeypatch.setattr(cli, "_enriquecer", lambda *a: [{"alvo_chembl_id": "CHEMBL222"}])
    real = atomico.temporario

    def falha(destino):
        if str(destino).endswith("enr.tsv"):
            raise OSError("disco cheio")
        return real(destino)
    monkeypatch.setattr(atomico, "temporario", falha)
    saida = tmp_path / "evidencias.tsv"
    saida.write_bytes(b"#antigo\nprecioso\n")
    rc = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
               "--alvos", "tests/dados/alvos_minimo.tsv",
               "--sintomas", "tests/dados/sintomas.txt", "--saida", str(saida),
               "--rotulos", str(rotulos),
               "--alvos-enriquecidos", str(tmp_path / "enr.tsv")])
    assert rc == 2
    assert saida.read_bytes() == b"#antigo\nprecioso\n"
    assert not (tmp_path / "manifesto.json").exists()
    assert not (tmp_path / "enr.tsv").exists()
    assert "alvos enriquecidos" in capsys.readouterr().err


def test_largura_do_checkpoint_divergente_falha_limpo(tmp_path, modelo_falso,
                                                      rotulos_falsos, capsys):
    """O README manda trocar o checkpoint: par divergente e o caso esperado.
    Antes isso dava zero linhas sem nenhuma mensagem."""
    saida = tmp_path / "evidencias.tsv"
    saida.write_bytes(b"#antigo\nprecioso\n")
    rc = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
               "--alvos", str(_alvos_do_modelo_falso(tmp_path)),
               "--sintomas", "tests/dados/sintomas.txt", "--saida", str(saida),
               "--checkpoint", modelo_falso,
               "--rotulos", str(_rotulos_tsv(tmp_path, rotulos_falsos[:3]))])
    assert rc == 2
    err = capsys.readouterr().err
    assert "espacos de alvos diferentes" in err and "Traceback" not in err
    assert saida.read_bytes() == b"#antigo\nprecioso\n"
    assert not (tmp_path / "manifesto.json").exists()


def test_o_manifesto_da_etapa_1_registra_o_dispositivo(tmp_path):
    """A fase 3 avisa quando o dispositivo muda; a etapa 1 nao registrava
    nenhum. Era a unica etapa que podia pegar a GPU em silencio E a unica cujo
    manifesto nao podia dizer que pegou -- nem o modo de repeticao detectaria.
    Agora ela e CPU por construcao, e o manifesto prova qual foi."""
    import json
    saida = tmp_path / "e.tsv"
    codigo = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
                   "--alvos", "tests/dados/alvos_minimo.tsv",
                   "--sintomas", "exemplos/sintomas.dor_neuropatica.txt",
                   "--saida", str(saida), "--sem-rede"])
    assert codigo in (0, 3)
    manifesto = json.loads(
        (tmp_path / "manifesto.json").read_text(encoding="utf-8"))
    assert manifesto["versoes"]["dispositivo"] == "cpu"


# ===========================================================================
# Espectro de alvos, decisao 3: o teto do modelo no manifesto
# ===========================================================================
def test_o_manifesto_carrega_o_teto_do_modelo(tmp_path):
    """Sem este numero, retreinar o preditor e opiniao. Ele sai medido a cada
    execucao, nao citado de uma execucao antiga."""
    saida = tmp_path / "evidencias.tsv"
    main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
          "--alvos", "tests/dados/alvos_minimo.tsv",
          "--sintomas", "tests/dados/sintomas.txt",
          "--saida", str(saida), "--sem-rede"])
    teto = ler_manifesto(tmp_path / "manifesto.json")["contagens"]["teto_do_modelo"]
    for chave in ("alvos_total", "alvos_no_modelo", "alvos_fora_do_modelo",
                  "alvos_marca_ausente", "alvos_fora_com_evidencia_medida",
                  "alvos_fora_sem_evidencia_medida", "dobrados_por_familia",
                  "pares_previstos_possiveis", "pares_previstos_perdidos"):
        assert chave in teto, chave


def test_o_teto_soma_com_a_contagem_de_alvos(tmp_path):
    """O teto nao pode contar outra lista de alvos que a que a execucao usou."""
    saida = tmp_path / "evidencias.tsv"
    main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
          "--alvos", "tests/dados/alvos_minimo.tsv",
          "--sintomas", "tests/dados/sintomas.txt",
          "--saida", str(saida), "--sem-rede"])
    c = ler_manifesto(tmp_path / "manifesto.json")["contagens"]
    assert c["teto_do_modelo"]["alvos_total"] == c["alvos"]
    # tres estados, nao dois: dentro, fora e marca ausente (lista de IDs sem a
    # coluna `no_modelo`, rodada sem o modelo em maos)
    assert (c["teto_do_modelo"]["alvos_no_modelo"]
            + c["teto_do_modelo"]["alvos_fora_do_modelo"]
            + c["teto_do_modelo"]["alvos_marca_ausente"]) == c["alvos"]


def test_o_manifesto_permite_notar_que_a_regra_mudou(tmp_path):
    """Os mesmos parametros passam a dar outra lista de alvos depois da
    decisao 2. `pipeline_commit` e `teto_do_modelo` sao o que distingue um
    manifesto de antes da mudanca de um de depois -- sem isso, duas execucoes
    "identicas" com listas diferentes seriam indistinguiveis."""
    saida = tmp_path / "evidencias.tsv"
    main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
          "--alvos", "tests/dados/alvos_minimo.tsv",
          "--sintomas", "tests/dados/sintomas.txt",
          "--saida", str(saida), "--sem-rede"])
    m = ler_manifesto(tmp_path / "manifesto.json")
    assert m["versoes"]["pipeline_commit"]
    assert "teto_do_modelo" in m["contagens"]


def test_a_chave_que_mede_a_decisao_2_nao_e_sempre_zero(tmp_path, monkeypatch):
    """I2: `reactome_fora_do_modelo_entraram` lia um contador aninhado e dava 0
    em toda execucao real -- indistinguivel de "a mudanca nao fez nada"."""
    from mods.evidencias.enriquecimento import reactome_fora_do_modelo
    assert reactome_fora_do_modelo(
        {"reactome_modelo": {"fora_do_modelo": 5}}) == 5
    import inspect

    from mods import coletar_evidencias
    fonte = inspect.getsource(coletar_evidencias)
    assert 'contagens_enriquecimento.get("fora_do_modelo"' not in fonte
    assert "reactome_fora_do_modelo(contagens_enriquecimento)" in fonte


def test_o_teto_da_lista_enriquecida_sai_quando_ela_e_escrita(tmp_path):
    """I3 e I7: na execucao que FAZ o enriquecimento, o teto sobre `--alvos`
    mede a lista PRE-enriquecimento -- zero alvos fora do modelo, zero
    dobrados. Quem seguir o README citaria esse numero e concluiria que o teto
    e zero. O segundo teto mede a lista que a execucao acabou de escrever."""
    import inspect

    from mods import coletar_evidencias
    fonte = inspect.getsource(coletar_evidencias)
    assert "teto_da_lista_enriquecida" in fonte
