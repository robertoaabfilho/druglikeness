# -*- coding: utf-8 -*-
"""
Compatibilidade com Windows.

**O que estes testes provam e o que nao provam.** Eles exercitam o CAMINHO do
Windows (`os.name == "nt"`) nesta maquina Linux, por substituicao. Isso pega
erro de logica e de API -- e foi assim que o problema grave abaixo apareceu.
O que eles NAO substituem e rodar a suite num Windows de verdade: o
comportamento do sistema de arquivos (`os.replace` sobre arquivo aberto), do
console (pagina de codigo) e do instalador do RDKit so se verifica la.

O problema grave que esta varredura encontrou
---------------------------------------------
No Windows, `os.kill(pid, 0)` **nao consulta: ele chama `TerminateProcess`**.
O painel perguntava "esta viva?" a cada 2 segundos pelo fragmento de
acompanhamento -- no Windows isso MATARIA a execucao do pipeline. A consulta
tem de ser `OpenProcess` + `GetExitCodeProcess`, que e o que o painel do
prospector ja fazia e que eu havia perdido na copia.
"""
import os
import sys
from pathlib import Path

import pytest

from mods import caminhos
from mods.painel import execucao


# --------------------------------------------------------------------------
# vivacidade: o erro que mataria a execucao
# --------------------------------------------------------------------------
def test_no_windows_a_consulta_de_vivacidade_nao_usa_os_kill(monkeypatch):
    """`os.kill(pid, 0)` no Windows chama TerminateProcess. Se `vivo()` o
    usar, o painel mata o pipeline ao perguntar se ele esta de pe."""
    monkeypatch.setattr(execucao, "_e_windows", lambda: True)
    chamou = []
    monkeypatch.setattr(os, "kill",
                        lambda *a, **k: chamou.append(a) or None)
    monkeypatch.setattr(execucao, "_vivo_nt", lambda pid: True)
    assert execucao.vivo(os.getpid()) is True
    assert chamou == [], "vivo() chamou os.kill no Windows"


def test_vivo_nt_usa_o_codigo_de_saida_do_processo(monkeypatch):
    """No Windows nao ha `/proc`, e e `GetExitCodeProcess` que distingue
    processo rodando (259, STILL_ACTIVE) de processo que ja saiu."""
    monkeypatch.setattr(execucao, "_e_windows", lambda: True)
    estados = {"codigo": 259}

    def falso(pid):
        return estados["codigo"] == 259
    monkeypatch.setattr(execucao, "_vivo_nt", falso)
    assert execucao.vivo(1234) is True
    estados["codigo"] = 0
    assert execucao.vivo(1234) is False


def test_vivo_nt_existe_e_nao_levanta_fora_do_windows():
    """A funcao tem de poder ser importada e chamada em qualquer sistema: a
    suite roda em Linux, e um `import ctypes.windll` no topo quebraria tudo."""
    assert execucao._vivo_nt(999_999_999) is False


def test_o_zumbi_e_conceito_de_posix_e_nao_atrapalha_no_windows(monkeypatch):
    """No Windows um processo que saiu nao fica zumbi esperando o pai: o
    proprio codigo de saida ja diz que ele terminou. `_zumbi` devolve False
    sem `/proc`, e o caminho do Windows nao depende dele."""
    monkeypatch.setattr(execucao, "_e_windows", lambda: True)
    assert execucao._zumbi(os.getpid()) is False


# --------------------------------------------------------------------------
# lancar o processo de modo que ele sobreviva a aba
# --------------------------------------------------------------------------
def test_no_windows_o_lancamento_nao_usa_start_new_session(monkeypatch,
                                                           tmp_path):
    """`start_new_session` e `setsid`, que nao existe no Windows. La o
    equivalente e `CREATE_NEW_PROCESS_GROUP | DETACHED_PROCESS`."""
    monkeypatch.setattr(caminhos, "EXECUCAO", tmp_path / ".execucao")
    monkeypatch.setattr(execucao, "_e_windows", lambda: True)
    vistos = {}

    class Falso:
        def __init__(self, *a, **k):
            vistos.update(k)
            self.pid = os.getpid()
    monkeypatch.setattr(execucao.subprocess, "Popen", Falso)
    monkeypatch.setattr(execucao, "_CRIACAO_NT", 0x00000200 | 0x00000008)
    execucao.iniciar(["etapa1"])
    assert "start_new_session" not in vistos
    assert vistos.get("creationflags") == 0x00000208


def test_no_posix_o_lancamento_continua_com_start_new_session(monkeypatch,
                                                              tmp_path):
    monkeypatch.setattr(caminhos, "EXECUCAO", tmp_path / ".execucao")
    vistos = {}

    class Falso:
        def __init__(self, *a, **k):
            vistos.update(k)
            self.pid = os.getpid()
    monkeypatch.setattr(execucao.subprocess, "Popen", Falso)
    execucao.iniciar(["etapa1"])
    assert vistos.get("start_new_session") is True
    assert "creationflags" not in vistos


def test_encerrar_no_windows_nao_manda_sigterm(monkeypatch, tmp_path):
    """`SIGTERM` no Windows tambem vira TerminateProcess -- que e o que
    queremos ao ENCERRAR, mas pelo caminho explicito, nao por acidente de
    portabilidade. A publicacao e atomica, entao o encerramento abrupto deixa
    o resultado anterior intacto."""
    import time
    monkeypatch.setattr(caminhos, "EXECUCAO", tmp_path / ".execucao")
    (tmp_path / ".execucao").mkdir(parents=True)
    (tmp_path / ".execucao" / "rodando.pid").write_text(
        f"{os.getpid()} {time.time():.0f} etapa1", encoding="utf-8")
    monkeypatch.setattr(execucao, "_e_windows", lambda: True)
    monkeypatch.setattr(execucao, "_vivo_nt", lambda pid: True)
    matou = []
    monkeypatch.setattr(execucao, "_encerrar_nt", lambda pid: matou.append(pid))
    monkeypatch.setattr(os, "kill", lambda *a: pytest.fail(
        "encerrar() usou os.kill no Windows"))
    execucao.encerrar()
    assert matou == [os.getpid()]


# --------------------------------------------------------------------------
# caminhos e URIs
# --------------------------------------------------------------------------
def test_o_indice_do_openfda_abre_por_uri_valida(tmp_path):
    """`f"file:{caminho}?mode=ro"` com um caminho do Windows (`C:\\x\\y`) nao e
    uma URI valida: as barras invertidas e os dois-pontos da unidade quebram.
    `Path.as_uri()` produz `file:///C:/x/y`, que vale nos dois sistemas."""
    import json
    import zipfile

    from mods import openfda_massa as massa
    bruto = tmp_path / "bruto"
    bruto.mkdir()
    with zipfile.ZipFile(bruto / "p.json.zip", "w") as z:
        z.writestr("p.json", json.dumps({"results": [
            {"openfda": {"generic_name": ["X"]}, "effective_time": "20260101",
             "warnings": ["nada"]}]}))
    indice = tmp_path / "i.sqlite"
    massa.indexar(bruto, indice, export_date="x")
    acervo = massa.Acervo(indice)
    assert acervo.versao == "openfda-massa x"
    fonte = Path("mods/openfda_massa.py").read_text(encoding="utf-8")
    assert 'f"file:{self.indice}?mode=ro"' not in fonte, \
        "a URI do SQLite ainda e montada por concatenacao"


def test_nenhum_caminho_e_montado_com_barra_literal():
    """Caminho montado como texto com `/` funciona no Windows por acidente do
    Python, mas quebra ao ser comparado ou exibido. Todo caminho sai de
    `caminhos.py` ou de `pathlib`."""
    import re
    suspeitos = []
    for arquivo in list(Path("mods").rglob("*.py")) + [Path("main.py"),
                                                       Path("painel.py")]:
        for numero, linha in enumerate(
                arquivo.read_text(encoding="utf-8").splitlines(), start=1):
            if linha.strip().startswith("#") or '"""' in linha:
                continue
            # os.path.join com barra literal, ou concatenacao de caminho
            if re.search(r'\+\s*"/"\s*\+', linha):
                suspeitos.append(f"{arquivo}:{numero}")
    assert not suspeitos, suspeitos


def test_o_orquestrador_usa_o_interpretador_corrente():
    """`python` no PATH do Windows pode ser outro (ou a loja da Microsoft).
    `sys.executable` e o unico que garante o mesmo ambiente das dependencias."""
    fonte = Path("main.py").read_text(encoding="utf-8")
    assert "sys.executable" in fonte
    assert '"python"' not in fonte.replace('COMANDO = "python main.py"', "")


# --------------------------------------------------------------------------
# escrita de arquivo
# --------------------------------------------------------------------------
def test_a_publicacao_atomica_tolera_falha_de_permissao(tmp_path, monkeypatch):
    """No Windows `os.chmod` cobre pouco e pode falhar; o ajuste de modo e
    cortesia e nunca pode custar o resultado."""
    from mods.evidencias import atomico
    destino = tmp_path / "x.tsv"
    tmp = atomico.temporario(destino)
    tmp.write_text("conteudo", encoding="utf-8")
    monkeypatch.setattr(os, "chmod", lambda *a, **k: (_ for _ in ()).throw(
        OSError("sem suporte")))
    atomico.publicar(tmp, destino)
    assert destino.read_text(encoding="utf-8") == "conteudo"


def test_toda_leitura_de_texto_declara_encoding():
    """Sem `encoding=`, o Windows usa a pagina de codigo do sistema (cp1252) e
    um nome de ligante com caractere fora dela vira erro ou lixo. O pipeline
    inteiro e UTF-8."""
    import re
    faltando = []
    for arquivo in list(Path("mods").rglob("*.py")) + [Path("main.py"),
                                                       Path("painel.py")]:
        texto = arquivo.read_text(encoding="utf-8")
        for numero, linha in enumerate(texto.splitlines(), start=1):
            if not re.search(r'\bopen\s*\(', linha):
                continue
            if re.search(r'["\'](r|w|a)b["\']|encoding\s*=|z\.open|'
                         r'urlopen|os\.open|gzip\.open|\.open\(interno', linha):
                continue
            # linha de continuacao: o encoding pode estar na proxima
            proxima = texto.splitlines()[numero] if numero < len(
                texto.splitlines()) else ""
            if "encoding" in proxima or '"rb"' in proxima or '"wb"' in proxima:
                continue
            faltando.append(f"{arquivo}:{numero}: {linha.strip()[:70]}")
    assert not faltando, faltando


def test_o_pid_e_lido_de_volta_com_fim_de_linha_do_windows(tmp_path,
                                                           monkeypatch):
    """O arquivo de PID e escrito em modo texto; no Windows o `\\n` vira
    `\\r\\n`. A leitura usa `split()`, que tolera -- este teste prende isso."""
    import time
    monkeypatch.setattr(caminhos, "EXECUCAO", tmp_path / ".execucao")
    (tmp_path / ".execucao").mkdir(parents=True)
    (tmp_path / ".execucao" / "rodando.pid").write_bytes(
        f"{os.getpid()} {time.time():.0f} etapa1\r\n".encode("utf-8"))
    emc = execucao.em_curso()
    assert emc is not None and emc["etapas"] == ["etapa1"]


def test_requirements_nao_fixa_pacote_so_de_linux():
    """Um requisito que so existe em Linux faria `pip install -r` falhar no
    Windows antes de qualquer coisa."""
    texto = Path("requirements.txt").read_text(encoding="utf-8")
    for linha in texto.splitlines():
        nome = linha.strip().lower()
        assert not nome.startswith(("uvloop", "python-prctl", "pyinotify")), nome
