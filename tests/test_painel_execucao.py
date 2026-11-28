# -*- coding: utf-8 -*-
"""
T3 -- disparar o main.py, saber se esta de pe, medir o ritmo.

Tres propriedades que estes testes prendem:

1. O painel NAO reimplementa a sequencia das etapas: ele monta flags do
   `main.py` (`--somente`, `--a-partir-de`) a partir da tabela `ETAPAS` dele.
2. Duas abas abertas no mesmo painel nao lancam duas execucoes sobre as
   mesmas saidas. A segunda encontra a primeira pelo arquivo de PID e recusa.
3. PID e reciclado pelo sistema operacional. O arquivo guarda PID *e*
   instante de inicio, e expira -- sem o teto, um PID reaproveitado trancaria
   o painel para sempre.
"""
import os
import time

import pytest

from mods import caminhos
from mods.painel import execucao


@pytest.fixture
def execucao_limpa(tmp_path, monkeypatch):
    """A pasta de estado da execucao numa pasta temporaria."""
    monkeypatch.setattr(caminhos, "TSV", tmp_path)
    monkeypatch.setattr(caminhos, "EXECUCAO", tmp_path / ".execucao")
    return tmp_path / ".execucao"


def _pid(pasta, pid, inicio, etapas="etapa1"):
    pasta.mkdir(parents=True, exist_ok=True)
    (pasta / "rodando.pid").write_text(f"{pid} {inicio:.0f} {etapas}",
                                       encoding="utf-8")
    return pasta / "rodando.pid"


# --------------------------------------------------------------------------
# vivacidade
# --------------------------------------------------------------------------
def test_vivo_e_falso_para_pid_inexistente():
    assert execucao.vivo(999_999_999) is False


def test_vivo_e_falso_para_pid_invalido():
    assert execucao.vivo(0) is False
    assert execucao.vivo(-1) is False


def test_vivo_e_verdadeiro_para_o_proprio_processo():
    assert execucao.vivo(os.getpid()) is True


def test_em_curso_e_none_sem_arquivo(execucao_limpa):
    assert execucao.em_curso() is None


def test_em_curso_ve_execucao_de_outra_aba(execucao_limpa):
    """A aba foi fechada e reaberta, ou outra aba disparou: o painel reencontra
    a execucao pelo arquivo, nao pela memoria da sessao."""
    _pid(execucao_limpa, os.getpid(), time.time(), "etapa1,portao")
    emc = execucao.em_curso()
    assert emc["pid"] == os.getpid()
    assert emc["etapas"] == ["etapa1", "portao"]


def test_em_curso_ignora_e_limpa_pid_antigo_demais(execucao_limpa):
    """PID vivo, inicio ha 25 h: e um PID reciclado, nao a nossa execucao."""
    arq = _pid(execucao_limpa, os.getpid(), time.time() - 25 * 3600)
    assert execucao.em_curso() is None
    assert not arq.exists()


def test_em_curso_ignora_e_limpa_pid_morto(execucao_limpa):
    arq = _pid(execucao_limpa, 999_999_999, time.time())
    assert execucao.em_curso() is None
    assert not arq.exists()


def test_em_curso_ignora_arquivo_corrompido(execucao_limpa):
    execucao_limpa.mkdir(parents=True)
    (execucao_limpa / "rodando.pid").write_text("lixo", encoding="utf-8")
    assert execucao.em_curso() is None


# --------------------------------------------------------------------------
# disparo
# --------------------------------------------------------------------------
class _ProcessoFalso:
    def __init__(self, *a, **k):
        self.argv = list(a[0])
        self.pid = os.getpid()
        self.kwargs = k


@pytest.fixture
def popen_falso(monkeypatch):
    """Captura o argv em vez de lancar processo de verdade."""
    capturados = []

    def falso(*a, **k):
        p = _ProcessoFalso(*a, **k)
        capturados.append(p)
        return p
    monkeypatch.setattr(execucao.subprocess, "Popen", falso)
    return capturados


def test_iniciar_recusa_segunda_execucao_simultanea(execucao_limpa,
                                                    popen_falso):
    _pid(execucao_limpa, os.getpid(), time.time())
    with pytest.raises(execucao.JaRodando):
        execucao.iniciar(["etapa1"])
    assert not popen_falso  # nada foi lancado


def test_iniciar_recusa_etapa_desconhecida(execucao_limpa, popen_falso):
    with pytest.raises(ValueError, match="etapa9"):
        execucao.iniciar(["etapa9"])


def test_iniciar_recusa_lista_vazia(execucao_limpa, popen_falso):
    with pytest.raises(ValueError):
        execucao.iniciar([])


def test_iniciar_de_uma_etapa_usa_somente(execucao_limpa, popen_falso):
    execucao.iniciar(["portao"])
    assert popen_falso[0].argv[-2:] == ["--somente", "portao"]


def test_iniciar_de_sufixo_contiguo_usa_a_partir_de(execucao_limpa,
                                                    popen_falso):
    execucao.iniciar(["portao", "ranking"])
    assert "--a-partir-de" in popen_falso[0].argv
    assert popen_falso[0].argv[-1] == "portao"


def test_iniciar_de_todas_nao_passa_flag_de_etapa(execucao_limpa, popen_falso):
    execucao.iniciar(list(caminhos.ETAPAS))
    argv = popen_falso[0].argv
    assert "--somente" not in argv and "--a-partir-de" not in argv


def test_iniciar_recusa_conjunto_nao_contiguo(execucao_limpa, popen_falso):
    """`main.py` sabe dizer "a partir de" e "so uma"; nao sabe dizer "a
    primeira e a terceira". Rodar a do meio sem ter sido pedida misturaria
    duas execucoes."""
    with pytest.raises(ValueError):
        execucao.iniciar(["etapa1", "ranking"])
    assert not popen_falso


def test_iniciar_chama_o_main_do_projeto(execucao_limpa, popen_falso):
    execucao.iniciar(["etapa1"])
    argv = popen_falso[0].argv
    assert argv[0] == os.sys.executable
    assert argv[1].endswith("main.py")
    assert popen_falso[0].kwargs["cwd"] == str(caminhos.RAIZ)


def test_iniciar_repassa_retomar_e_config(execucao_limpa, popen_falso):
    execucao.iniciar(["etapa1"], retomar=True, config="outro.json")
    argv = popen_falso[0].argv
    assert "--retomar" in argv
    assert argv[argv.index("--config") + 1] == "outro.json"


def test_iniciar_grava_pid_e_abre_o_log(execucao_limpa, popen_falso):
    estado = execucao.iniciar(["etapa1"])
    assert (execucao_limpa / "rodando.pid").exists()
    assert estado["pid"] == os.getpid()
    assert estado["log"] == execucao_limpa / "main.log"
    assert estado["log"].exists()


def test_iniciar_sobrevive_ao_fechamento_da_aba(execucao_limpa, popen_falso):
    """`start_new_session` tira o processo do grupo do Streamlit: fechar a aba
    (ou o servidor) nao mata uma execucao de horas."""
    execucao.iniciar(["etapa1"])
    assert popen_falso[0].kwargs["start_new_session"] is True


# --------------------------------------------------------------------------
# resultado e ritmo
# --------------------------------------------------------------------------
def test_ultima_execucao_e_none_sem_registro(execucao_limpa):
    assert execucao.ultima_execucao() is None


def test_registrar_fim_guarda_o_codigo_de_saida(execucao_limpa, popen_falso):
    execucao.iniciar(["etapa1"])
    execucao.registrar_fim(3)
    ultima = execucao.ultima_execucao()
    assert ultima["codigo"] == 3
    assert ultima["etapas"] == ["etapa1"]
    assert not (execucao_limpa / "rodando.pid").exists()


def test_cauda_do_log_devolve_as_ultimas_linhas(execucao_limpa):
    execucao_limpa.mkdir(parents=True)
    (execucao_limpa / "main.log").write_text(
        "".join(f"linha {i}\n" for i in range(100)), encoding="utf-8")
    cauda = execucao.cauda_do_log(linhas=5)
    assert cauda.splitlines() == [f"linha {i}" for i in range(95, 100)]


def test_cauda_do_log_sem_log_devolve_vazio(execucao_limpa):
    assert execucao.cauda_do_log() == ""


def test_ritmo_precisa_de_duas_amostras(execucao_limpa):
    execucao.registrar_ritmo(100)
    assert execucao.ritmo_por_hora() is None


def test_ritmo_por_hora_mede_a_taxa_observada(execucao_limpa):
    """Duas amostras, 1800 s de distancia, 500 registros: 1.000 por hora.
    Taxa medida, nao barra falsa -- o painel antigo avisava que a barra dele
    era estimativa."""
    execucao_limpa.mkdir(parents=True)
    agora = time.time()
    (execucao_limpa / "ritmo.csv").write_text(
        f"{agora - 1800:.0f},0\n{agora:.0f},500\n", encoding="utf-8")
    assert execucao.ritmo_por_hora() == pytest.approx(1000.0, rel=0.02)


def test_ritmo_ignora_amostra_velha(execucao_limpa):
    """Amostra de 20 min atras: a execucao parou, e uma taxa que nao esta mais
    acontecendo e pior que nenhuma."""
    execucao_limpa.mkdir(parents=True)
    agora = time.time()
    (execucao_limpa / "ritmo.csv").write_text(
        f"{agora - 3600:.0f},0\n{agora - 1200:.0f},500\n", encoding="utf-8")
    assert execucao.ritmo_por_hora() is None


def test_registrar_ritmo_nao_duplica_amostra_sem_avanco(execucao_limpa):
    execucao.registrar_ritmo(100)
    execucao.registrar_ritmo(100)
    linhas = (execucao_limpa / "ritmo.csv").read_text().strip().splitlines()
    assert len(linhas) == 1


def test_encerrar_sem_execucao_nao_levanta(execucao_limpa):
    execucao.encerrar()  # nada a encerrar: silencio, nao erro


# ===========================================================================
# Correcoes da revisao do ramo
# ===========================================================================
_FILHOS = []  # mantem o Popen vivo: coletado, o subprocess colhe o zumbi


def _filho_que_sai(codigo: int) -> int:
    """Um filho DE VERDADE deste processo, ja terminado e nao colhido -- isto e,
    zumbi. E o estado em que toda execucao do painel termina, porque o
    `main.py` e filho do processo do Streamlit."""
    import subprocess
    processo = subprocess.Popen(
        [os.sys.executable, "-c", f"raise SystemExit({codigo})"])
    _FILHOS.append(processo)
    while True:  # espera o processo sair SEM o colher
        with open(f"/proc/{processo.pid}/stat", encoding="utf-8") as arq:
            if arq.read().rsplit(")", 1)[1].split()[0] == "Z":
                return processo.pid


def test_vivo_e_falso_para_filho_zumbi():
    """C1: `os.kill(pid, 0)` num zumbi devolve sucesso. Sem olhar o estado do
    processo, o painel declara viva toda execucao que terminou -- e o botao
    "Rodar" fica desligado por 24 h."""
    pid = _filho_que_sai(0)
    assert execucao.vivo(pid) is False


def test_colher_registra_o_codigo_de_saida_real(execucao_limpa):
    """C1+I4: quem observa o fim e o painel, e colher o filho e o que da o
    codigo de saida de verdade -- 3 (parcial) tem de chegar na tela como 3."""
    pid = _filho_que_sai(3)
    _pid(execucao_limpa, pid, time.time(), "etapa1,portao")
    assert execucao.colher()["codigo"] == 3
    ultima = execucao.ultima_execucao()
    assert ultima["codigo"] == 3
    assert ultima["etapas"] == ["etapa1", "portao"]
    assert ultima["inicio"] is not None  # I4: nao perde o inicio


def test_colher_libera_o_painel(execucao_limpa):
    pid = _filho_que_sai(0)
    _pid(execucao_limpa, pid, time.time())
    execucao.colher()
    assert execucao.em_curso() is None
    assert not (execucao_limpa / "rodando.pid").exists()


def test_colher_nao_mexe_em_execucao_viva(execucao_limpa):
    _pid(execucao_limpa, os.getpid(), time.time())
    assert execucao.colher() is None
    assert (execucao_limpa / "rodando.pid").exists()


def test_colher_sem_execucao_nao_levanta(execucao_limpa):
    assert execucao.colher() is None


def test_colher_de_processo_que_nao_e_nosso_filho(execucao_limpa):
    """Painel reiniciado: o `main.py` morreu e foi colhido pelo init. Nao da
    para saber o codigo, e o relato diz isso em vez de inventar 0."""
    _pid(execucao_limpa, 999_999_999, time.time())
    registro = execucao.colher()
    assert registro is not None and registro["codigo"] is None


def test_em_curso_nao_apaga_a_tranca_em_erro_de_leitura(execucao_limpa,
                                                        monkeypatch):
    """I5: apagar a tranca em qualquer OSError (EIO, EINTR) libera o botao
    "Rodar" sobre uma execucao que esta de pe."""
    arq = _pid(execucao_limpa, os.getpid(), time.time())
    original = type(arq).read_text

    def falha(self, *a, **k):
        if self.name == "rodando.pid":
            raise OSError(5, "erro de entrada/saida")
        return original(self, *a, **k)
    monkeypatch.setattr(type(arq), "read_text", falha)
    assert execucao.em_curso() is None
    assert arq.exists()  # a tranca ficou


def test_em_curso_apaga_tranca_corrompida(execucao_limpa):
    """Conteudo que nao se le nunca vai se ler: apagar e o unico jeito de o
    painel voltar a funcionar."""
    execucao_limpa.mkdir(parents=True, exist_ok=True)
    arq = execucao_limpa / "rodando.pid"
    arq.write_text("nao e um pid", encoding="utf-8")
    assert execucao.em_curso() is None
    assert not arq.exists()


def test_iniciar_cria_a_tranca_antes_de_lancar(execucao_limpa, monkeypatch):
    """I5: a tranca nasce com O_EXCL ANTES do Popen. Entre conferir e gravar
    havia o Popen inteiro, e duas abas passavam as duas pelo `em_curso() is
    None` -- dois main.py sobre as mesmas saidas."""
    vista = {}

    def falso(*a, **k):
        vista["tranca"] = (caminhos.EXECUCAO / "rodando.pid").exists()
        return _ProcessoFalso(*a, **k)
    monkeypatch.setattr(execucao.subprocess, "Popen", falso)
    execucao.iniciar(["etapa1"])
    assert vista["tranca"] is True


def test_a_tranca_e_exclusiva(execucao_limpa, popen_falso):
    """Duas chamadas simultaneas: a segunda nao cria a tranca, e recusa."""
    execucao.iniciar(["etapa1"])
    (caminhos.EXECUCAO / "rodando.pid").write_text(
        f"{os.getpid()} {time.time():.0f} etapa1", encoding="utf-8")
    with pytest.raises(execucao.JaRodando):
        execucao.iniciar(["portao"])
    assert len(popen_falso) == 1


def test_tranca_sem_pid_e_respeitada_por_um_minuto(execucao_limpa):
    """A tranca nasce antes de o PID existir. Nessa janela ela vale; passada,
    e resto de um lancamento que falhou."""
    execucao_limpa.mkdir(parents=True, exist_ok=True)
    arq = execucao_limpa / "rodando.pid"
    arq.write_text(f"-1 {time.time():.0f} etapa1", encoding="utf-8")
    assert execucao.em_curso() is not None
    arq.write_text(f"-1 {time.time() - 300:.0f} etapa1", encoding="utf-8")
    assert execucao.em_curso() is None


def test_iniciar_nao_trunca_a_tranca_em_lugar(execucao_limpa, popen_falso):
    """Gravar o PID por troca de nome, nao em modo 'w': em modo 'w' existe uma
    janela de 0 bytes em que outra aba le a tranca vazia e a apaga."""
    execucao.iniciar(["etapa1"])
    assert execucao.em_curso()["pid"] == os.getpid()


def test_ritmo_e_separado_por_etapa(execucao_limpa):
    """I3: as tres etapas escreviam no mesmo ritmo.csv, sem coluna de etapa;
    quando a etapa 2 comeca a contar de zero, a taxa mistura duas etapas."""
    execucao.registrar_ritmo(100, etapa="etapa1")
    execucao.registrar_ritmo(5, etapa="portao")
    assert (execucao_limpa / "ritmo_etapa1.csv").exists()
    assert (execucao_limpa / "ritmo_portao.csv").exists()


def test_ritmo_por_hora_le_a_etapa_pedida(execucao_limpa):
    execucao_limpa.mkdir(parents=True, exist_ok=True)
    agora = time.time()
    (execucao_limpa / "ritmo_etapa1.csv").write_text(
        f"{agora - 1800:.0f},0\n{agora:.0f},500\n", encoding="utf-8")
    assert execucao.ritmo_por_hora(etapa="etapa1") == pytest.approx(1000.0,
                                                                    rel=0.02)
    assert execucao.ritmo_por_hora(etapa="portao") is None


def test_cauda_do_log_nao_le_o_log_inteiro(execucao_limpa, monkeypatch):
    """Minor 8: o log cresce por horas e o fragmento o relia a cada 2 s."""
    execucao_limpa.mkdir(parents=True, exist_ok=True)
    grande = execucao_limpa / "main.log"
    grande.write_text("".join(f"linha {i}\n" for i in range(200_000)),
                      encoding="utf-8")
    lidos = {"bytes": 0}
    real = open

    def contando(arquivo, *a, **k):
        f = real(arquivo, *a, **k)
        if str(arquivo).endswith("main.log"):
            leitura_original = f.read

            def ler(n=-1):
                dados = leitura_original(n)
                lidos["bytes"] += len(dados)
                return dados
            f.read = ler
        return f
    monkeypatch.setattr("builtins.open", contando)
    cauda = execucao.cauda_do_log(linhas=5)
    assert cauda.strip().splitlines()[-1] == "linha 199999"
    assert lidos["bytes"] < 200_000, f"leu {lidos['bytes']} bytes"


def test_o_painel_conta_o_temporario_e_nao_a_copia(execucao_limpa, tmp_path):
    """I3: `registrar_ritmo` recebia o numero de linhas da COPIA, que tem teto
    de 2.000. Numa execucao de ~23 milhoes de pares a contagem chegava a 2.000
    em segundos e nunca mais mudava: `registrar_ritmo` saia por
    `contagem == n0`, a segunda amostra nunca era gravada, e o ritmo observado
    -- o mecanismo 4 da especificacao -- ficava `None` para sempre.

    Este teste prende a fonte da contagem no painel, que e onde o erro estava:
    o temporario inteiro, nao a copia limitada."""
    from mods import caminhos as c
    fonte = (c.RAIZ / "painel.py").read_text(encoding="utf-8")
    assert "registrar_ritmo(n, " not in fonte
    assert "contar_linhas(leitura.temporario_da_etapa(" in fonte or \
        "contar_linhas(tmp" in fonte, \
        "o ritmo tem de vir da contagem do temporario, nao do tamanho da copia"
