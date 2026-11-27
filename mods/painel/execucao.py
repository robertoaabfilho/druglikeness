# -*- coding: utf-8 -*-
"""
execucao.py - disparar o pipeline e acompanhar, sem reimplementa-lo.

O painel lanca `python main.py` como processo separado e **nao** repete a
sequencia das etapas: ele monta as flags do orquestrador a partir da tabela
`ETAPAS` dele. Se a sequencia virasse codigo de interface, sairia do alcance
de qualquer script -- e e exatamente o erro que o desenho nomeia.

Tres coisas que nao sao conveniencia
------------------------------------
**O estado da execucao mora em arquivo, nao na sessao.** `st.session_state`
morre quando a aba fecha, e uma execucao destas dura horas. Em
`tsv/.execucao/` ficam `rodando.pid`, `main.log`, `ritmo.csv` e `ultimo.json`:
e isso que faz o painel reencontrar uma execucao que outra aba disparou, ou
que ele mesmo disparou antes de ser reiniciado.

**O arquivo de PID guarda PID *e* instante de inicio, e expira.** PID e
reciclado pelo sistema operacional: sem o teto, um PID reaproveitado por
outro programa qualquer trancaria o painel para sempre, e o remedio seria
apagar um arquivo que o pesquisador nao sabe que existe. O teto e 24 h, nao
as 6 h do painel de onde o mecanismo vem: a etapa 1 real sao ~23 milhoes de
pares, e tempo nao e restricao deste projeto.

**O progresso e taxa medida, nao barra.** `ritmo.csv` guarda (instante,
contagem) e a taxa sai da diferenca entre duas amostras distantes. Uma barra
de progresso precisa saber o total, e aqui o total depende de quantos pares
sobrevivem a cada filtro -- uma barra assim seria invencao apresentada como
medida.
"""
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

from mods import caminhos
from mods.evidencias import atomico

# Teto de idade do arquivo de PID. Ver o docstring: PID e reciclado.
VALIDADE_PID_H = 24
# Idade maxima da amostra mais recente para a taxa ainda valer.
FRESCOR_RITMO_MIN = 15
# Distancia minima entre as duas amostras que formam a taxa.
JANELA_MINIMA_S = 600
# Intervalo minimo entre amostras gravadas.
INTERVALO_AMOSTRA_S = 120
# Quanto vale uma tranca recem-criada que ainda nao tem PID. Passado isso, ela
# e resto de um lancamento que falhou antes de o processo existir.
JANELA_DE_PARTIDA_S = 60
# Quanto do fim do log a cauda le. O log cresce por horas e nunca e rotado;
# `aplicar_portao.py` imprime uma linha por problema de entrada.
CAUDA_BYTES = 64 * 1024


class JaRodando(Exception):
    """Ha execucao de pe. Quem captura mostra qual, e nao lanca a segunda."""


def _pasta() -> Path:
    pasta = caminhos.EXECUCAO
    pasta.mkdir(parents=True, exist_ok=True)
    return pasta


def arquivo_pid() -> Path:
    return caminhos.EXECUCAO / "rodando.pid"


def arquivo_log() -> Path:
    return caminhos.EXECUCAO / "main.log"


def arquivo_ritmo(etapa: str = "") -> Path:
    """Um arquivo de ritmo POR ETAPA. Num arquivo so, a etapa 2 comecando a
    contar de zero com as amostras da etapa 1 ainda na janela faria a taxa
    misturar duas etapas."""
    return caminhos.EXECUCAO / (f"ritmo_{etapa}.csv" if etapa else "ritmo.csv")


def arquivo_ultimo() -> Path:
    return caminhos.EXECUCAO / "ultimo.json"


# Sinalizadores do Windows para destacar o processo: ele continua de pe
# quando o painel (ou o terminal) fecha. `CREATE_NEW_PROCESS_GROUP |
# DETACHED_PROCESS` e o equivalente de `setsid`, que nao existe la.
_CRIACAO_NT = 0x00000200 | 0x00000008
_AINDA_ATIVO = 259  # STILL_ACTIVE, de GetExitCodeProcess


def _e_windows() -> bool:
    """O sistema, numa funcao so.

    Existe como funcao, e nao como `os.name == "nt"` espalhado, para que a
    suite possa exercitar o caminho do Windows nesta maquina sem trocar
    `os.name` globalmente -- o que quebraria o `pathlib`."""
    return os.name == "nt"


def _vivo_nt(pid: int) -> bool:
    """O processo existe e ainda esta rodando, no Windows.

    **Nao se usa `os.kill` aqui.** No Windows `os.kill(pid, 0)` nao consulta:
    ele chama `TerminateProcess`. O painel pergunta "esta viva?" a cada 2
    segundos pelo fragmento de acompanhamento -- com `os.kill` isso MATARIA a
    execucao do pipeline. A consulta certa e `OpenProcess` +
    `GetExitCodeProcess`, e so ela.

    Fora do Windows devolve False: a funcao existe em todo sistema para que a
    suite possa exercitar este caminho, e o `import ctypes.windll` no topo
    quebraria o import em Linux."""
    if not _e_windows() or pid <= 0:
        return False
    try:
        import ctypes
        kernel = ctypes.windll.kernel32  # type: ignore[attr-defined]
    except (ImportError, AttributeError, OSError):
        return False
    PROCESSO_CONSULTA = 0x1000  # PROCESS_QUERY_LIMITED_INFORMATION
    handle = kernel.OpenProcess(PROCESSO_CONSULTA, False, pid)
    if not handle:
        return False
    try:
        codigo = ctypes.c_ulong()
        ok = kernel.GetExitCodeProcess(handle, ctypes.byref(codigo))
        return bool(ok) and codigo.value == _AINDA_ATIVO
    finally:
        kernel.CloseHandle(handle)


def _encerrar_nt(pid: int) -> None:
    """Encerra o processo no Windows.

    Aqui `TerminateProcess` e o que se QUER -- mas pelo caminho explicito, nao
    por acidente de portabilidade do `os.kill`. A publicacao e atomica, entao
    o encerramento abrupto deixa o resultado anterior intacto."""
    try:
        import ctypes
        kernel = ctypes.windll.kernel32  # type: ignore[attr-defined]
    except (ImportError, AttributeError, OSError):
        return
    PROCESSO_ENCERRA = 0x0001  # PROCESS_TERMINATE
    handle = kernel.OpenProcess(PROCESSO_ENCERRA, False, pid)
    if not handle:
        return
    try:
        kernel.TerminateProcess(handle, 1)
    finally:
        kernel.CloseHandle(handle)


def _zumbi(pid: int) -> bool:
    """O processo saiu e ninguem colheu o codigo dele ainda.

    O campo 3 de `/proc/<pid>/stat` e o estado; `Z` e zumbi. O nome do
    executavel esta entre parenteses e pode conter espaco, por isso o corte e
    pelo ULTIMO `)`. Fora do Linux nao ha `/proc`: ali um zumbi so e detectado
    quando `colher` o recolhe."""
    try:
        with open(f"/proc/{pid}/stat", encoding="utf-8") as arq:
            return arq.read().rsplit(")", 1)[1].split()[0] == "Z"
    except (OSError, IndexError):
        return False


def vivo(pid: int) -> bool:
    """O processo ainda existe E ainda esta rodando.

    **Zumbi nao conta como vivo.** `os.kill(pid, 0)` num zumbi devolve
    sucesso: o processo saiu, mas a entrada sobrevive na tabela do sistema ate
    o pai colher o codigo de saida. E o estado em que TODA execucao disparada
    pelo painel termina, porque o `main.py` e filho do processo do Streamlit.
    Tratar zumbi como vivo deixava o painel declarando de pe cada execucao que
    acabou -- e o botao "Rodar" desligado por 24 h."""
    if pid <= 0:
        return False
    if _e_windows():
        return _vivo_nt(pid)
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True  # existe, e de outro usuario
    return not _zumbi(pid)


def _ler_pid() -> dict | None:
    """O conteudo da tranca, sem julgar se a execucao esta viva.

    Conteudo ilegivel e apagado: ele nunca vai passar a se ler, e deixa-lo
    trancaria o painel para sempre. Erro de LEITURA (EIO, EINTR, permissao)
    nao apaga nada -- apagar a tranca por um erro transitorio liberaria o
    botao "Rodar" sobre uma execucao que esta de pe."""
    arq = arquivo_pid()
    try:
        bruto = arq.read_text(encoding="utf-8")
    except FileNotFoundError:
        return None
    except OSError:
        return None  # transitorio: a tranca fica
    try:
        partes = bruto.split()
        pid, inicio = int(partes[0]), float(partes[1])
    except (IndexError, ValueError):
        arq.unlink(missing_ok=True)  # corrompido: nao se le nunca mais
        return None
    return {"pid": pid, "inicio": inicio,
            "etapas": partes[2].split(",") if len(partes) > 2 else [],
            "log": arquivo_log()}


def em_curso() -> dict | None:
    """A execucao de pe, se ha; `None` se nao.

    Limpa a tranca quando ela descreve execucao morta ou expirada. PID `-1` e
    a tranca recem-criada, antes de o processo existir: ela vale por
    `JANELA_DE_PARTIDA_S`, e depois disso e resto de um lancamento que
    falhou."""
    info = _ler_pid()
    if info is None:
        return None
    idade = time.time() - info["inicio"]
    if info["pid"] <= 0:
        if idade < JANELA_DE_PARTIDA_S:
            return info  # lancando agora
        arquivo_pid().unlink(missing_ok=True)
        return None
    if idade >= VALIDADE_PID_H * 3600 or not vivo(info["pid"]):
        arquivo_pid().unlink(missing_ok=True)
        return None
    return info


def colher() -> dict | None:
    """Fecha a execucao que terminou e devolve o registro dela.

    E aqui que o codigo de saida do `main.py` chega a tela. O `main.py` nao
    sabe que foi disparado pelo painel e nao e tarefa dele saber; quem observa
    o fim e a aba. Colher o filho tambem e o que **destrava o painel**: sem
    isso ele fica zumbi e `vivo()` responderia sim para sempre.

    `codigo` sai `None` quando nao da para saber -- painel reiniciado, e o
    `main.py` ja colhido pelo sistema. Inventar 0 diria "completa" sobre uma
    execucao que pode ter falhado."""
    info = _ler_pid()
    if info is None:
        return None
    if info["pid"] > 0 and vivo(info["pid"]):
        return None
    if info["pid"] <= 0 and time.time() - info["inicio"] < JANELA_DE_PARTIDA_S:
        return None
    codigo = None
    if info["pid"] > 0:
        try:
            morto, estado = os.waitpid(info["pid"], os.WNOHANG)
            if morto == info["pid"]:
                codigo = os.waitstatus_to_exitcode(estado)
        except (ChildProcessError, OSError):
            pass  # nao e nosso filho: o sistema ja o colheu
    return registrar_fim(codigo, info)


def _flags_das_etapas(escolhidas: list[str]) -> list[str]:
    """As flags do `main.py` que expressam esse conjunto de etapas.

    A forma mais estreita que serve: uma etapa -> `--somente`; um sufixo
    contiguo da ordem -> `--a-partir-de`; todas -> nenhuma flag. Conjunto nao
    contiguo e erro: o orquestrador nao sabe dizer "a primeira e a terceira",
    e rodar a do meio sem ter sido pedida entregaria resultado de uma etapa
    que o pesquisador nao mandou rodar."""
    ordem = list(caminhos.ETAPAS)
    if not escolhidas:
        raise ValueError("nenhuma etapa escolhida")
    for nome in escolhidas:
        if nome not in ordem:
            raise ValueError(
                f"etapa desconhecida: {nome!r}; as etapas sao "
                f"{', '.join(ordem)}")
    indices = sorted(ordem.index(n) for n in escolhidas)
    if len(indices) == 1:
        return ["--somente", ordem[indices[0]]]
    contiguo = indices == list(range(indices[0], indices[-1] + 1))
    if not contiguo or indices[-1] != len(ordem) - 1:
        raise ValueError(
            "o orquestrador sabe rodar uma etapa (`--somente`) ou da etapa N "
            "ate o fim (`--a-partir-de`); "
            f"{', '.join(ordem[i] for i in indices)} nao e nem um nem outro. "
            "Rode em duas vezes.")
    if indices[0] == 0:
        return []
    return ["--a-partir-de", ordem[indices[0]]]


def iniciar(etapas: list[str], *, config=None, retomar: bool = False) -> dict:
    """Lanca `main.py` em segundo plano e devolve o estado da execucao.

    Levanta `JaRodando` se ja ha uma de pe: duas abas do painel apertando
    "Rodar" lancariam dois processos sobre as mesmas saidas, e a publicacao
    atomica protege cada arquivo, nao a coerencia entre eles."""
    flags = _flags_das_etapas(etapas)  # valida antes de qualquer efeito
    colher()  # fecha execucao anterior que terminou e nao foi registrada
    atual = em_curso()
    if atual is not None:
        raise JaRodando(
            f"ja ha execucao de pe (PID {atual['pid']}, "
            f"{', '.join(atual['etapas']) or 'todas as etapas'})")
    _pasta()
    inicio = time.time()
    # A tranca nasce ANTES do Popen, e exclusiva. Antes, entre conferir e
    # gravar havia o Popen inteiro: duas abas passavam as duas pelo
    # `em_curso() is None` e lancavam dois `main.py` sobre as mesmas saidas.
    # A publicacao atomica protege cada arquivo, nao a coerencia entre eles.
    try:
        fd = os.open(arquivo_pid(), os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    except FileExistsError:
        raise JaRodando("outra aba comecou uma execucao neste instante") from None
    with os.fdopen(fd, "w", encoding="utf-8") as tranca:
        tranca.write(f"-1 {inicio:.0f} {','.join(etapas)}")
    argv = [sys.executable, str(caminhos.RAIZ / "main.py")]
    if config:
        argv += ["--config", str(config)]
    if retomar:
        argv.append("--retomar")
    argv += flags
    log = arquivo_log()
    try:
        with open(log, "a", encoding="utf-8") as saida:
            saida.write(f"\n===== {time.strftime('%Y-%m-%d %H:%M:%S')} "
                        f"{' '.join(argv[1:])}\n")
            saida.flush()
            # Destacar o processo para ele sobreviver ao fechamento da aba:
            # `setsid` no POSIX, sinalizadores de criacao no Windows.
            destacar = ({"creationflags": _CRIACAO_NT} if _e_windows()
                        else {"start_new_session": True})
            processo = subprocess.Popen(
                argv, cwd=str(caminhos.RAIZ), stdout=saida,
                stderr=subprocess.STDOUT, **destacar)
    except BaseException:
        arquivo_pid().unlink(missing_ok=True)  # a tranca nao sobrevive ao erro
        raise
    # Troca de nome, nao modo "w": em "w" existe uma janela de 0 bytes em que
    # outra aba le a tranca vazia.
    atomico.escrever_texto(arquivo_pid(),
                           f"{processo.pid} {inicio:.0f} {','.join(etapas)}")
    for etapa in etapas:  # a taxa e desta execucao
        arquivo_ritmo(etapa).unlink(missing_ok=True)
    return {"pid": processo.pid, "inicio": inicio, "etapas": list(etapas),
            "log": log}


def registrar_fim(codigo: int | None, info: dict | None = None) -> dict:
    """Guarda o resultado da execucao que acabou e libera o painel.

    `info` e o que a tranca dizia; sem ele, le a tranca. Nao se le a tranca
    por `em_curso()`: `em_curso` APAGA a tranca quando o processo esta morto,
    que e exatamente o caso aqui -- e o registro sairia sem `inicio` e sem
    `etapas`."""
    info = info or _ler_pid() or {}
    dados = {"codigo": codigo, "inicio": info.get("inicio"),
             "fim": time.time(), "etapas": info.get("etapas", []),
             "pid": info.get("pid")}
    _pasta()
    atomico.escrever_texto(arquivo_ultimo(),
                           json.dumps(dados, ensure_ascii=False, indent=2))
    arquivo_pid().unlink(missing_ok=True)
    return dados


def ultima_execucao() -> dict | None:
    try:
        return json.loads(arquivo_ultimo().read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def encerrar() -> None:
    """Pede ao processo que termine. Sem execucao de pe, silencio.

    SIGTERM, nao SIGKILL: a etapa grava em temporario e publica por troca de
    nome, e um encerramento ordenado deixa o resultado anterior intacto."""
    atual = em_curso()
    if atual is None:
        return
    if _e_windows():
        _encerrar_nt(atual["pid"])
    else:
        try:
            os.kill(atual["pid"], signal.SIGTERM)
        except (ProcessLookupError, PermissionError):
            pass
    arquivo_pid().unlink(missing_ok=True)


def cauda_do_log(linhas: int = 40) -> str:
    """As ultimas linhas do log, lendo so o fim do arquivo.

    O fragmento de acompanhamento chama isto a cada 2 s, e o log de uma
    execucao de ~23 milhoes de pares nao cabe numa releitura inteira a cada
    vez."""
    try:
        with open(arquivo_log(), "rb") as arq:
            arq.seek(0, os.SEEK_END)
            tamanho = arq.tell()
            arq.seek(max(0, tamanho - CAUDA_BYTES))
            bruto = arq.read()
    except OSError:
        return ""
    texto = bruto.decode("utf-8", errors="replace")
    if len(bruto) == CAUDA_BYTES:
        texto = texto.split("\n", 1)[-1]  # a primeira linha veio cortada
    return "\n".join(texto.splitlines()[-linhas:])


def registrar_ritmo(contagem: int, etapa: str = "") -> None:
    """Grava (instante, contagem) para medir a taxa, no maximo a cada 2 min.

    Sem avanco na contagem nao grava: amostras repetidas achatariam a taxa
    para zero sem que nada tivesse parado."""
    arq = arquivo_ritmo(etapa)
    try:
        ultima = arq.read_text(encoding="utf-8").strip().splitlines()[-1]
        t0, n0 = float(ultima.split(",")[0]), int(ultima.split(",")[1])
    except (OSError, IndexError, ValueError):
        t0, n0 = 0.0, -1
    if contagem == n0 or time.time() - t0 < INTERVALO_AMOSTRA_S:
        return
    _pasta()
    with open(arq, "a", encoding="utf-8") as saida:
        saida.write(f"{time.time():.0f},{contagem}\n")


def ritmo_por_hora(janela_min: int = 60, etapa: str = "") -> float | None:
    """Registros por hora na ultima janela; `None` se nao da para medir.

    Exige duas amostras com `JANELA_MINIMA_S` de distancia, e a mais recente
    com menos de `FRESCOR_RITMO_MIN`. Uma taxa que nao esta mais acontecendo
    e pior que nenhuma: ela faria a tela prometer um fim que nao vem."""
    try:
        linhas = [l.split(",") for l in
                  arquivo_ritmo(etapa).read_text(
                      encoding="utf-8").strip().splitlines()]
        amostras = [(float(t), int(n)) for t, n in linhas]
    except (OSError, ValueError):
        return None
    if len(amostras) < 2 or time.time() - amostras[-1][0] > FRESCOR_RITMO_MIN * 60:
        return None
    corte = amostras[-1][0] - janela_min * 60
    janela = [a for a in amostras if a[0] >= corte]
    if len(janela) < 2:
        return None
    dt = janela[-1][0] - janela[0][0]
    dn = janela[-1][1] - janela[0][1]
    return dn / dt * 3600 if dt >= JANELA_MINIMA_S and dn > 0 else None
