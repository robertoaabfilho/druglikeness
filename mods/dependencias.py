# -*- coding: utf-8 -*-
"""
dependencias.py - conferir o ambiente antes de rodar, e instalar o que falta.

O pipeline comeca por aqui. Uma execucao real sao ~23 milhoes de pares; um
`ModuleNotFoundError` no meio dela custa horas por algo que custa
milissegundos conferir -- a mesma razao de `main._precheca`.

Tres decisoes que este modulo toma, e por que
---------------------------------------------
**So o que `requirements.txt` declara.** A instalacao usa a LINHA do arquivo,
nao o nome do pacote: instalar `numpy` onde o arquivo pede `numpy>=2`
satisfaria o nome e nao a restricao, e o pacote `pp4ckd` de origem exige
`numpy<2` -- e exatamente o conflito que a restricao existe para evitar.

**Diz antes de instalar.** Instalacao silenciosa num ambiente que o
pesquisador compartilha com outro trabalho e a forma mais rapida de quebrar o
outro trabalho.

**O pesado NAO entra sozinho.** `torch` sao alguns GB, e o pipeline **degrada**
sem ele: a etapa 1 roda so com evidencia medida e diz o que desligou. Baixar
GB sem pedir seria decidir pelo pesquisador num ponto em que ele tem
alternativa. O que este modulo faz com uma opcional ausente e **dizer o que
ela desliga**.
"""
import importlib.util
import subprocess
import sys
from pathlib import Path

from mods import caminhos

# Pacote cujo nome de import difere do nome de instalacao. O mapa existe para
# os casos em que diverge; para os demais o nome vale para os dois.
MODULO_DO_PACOTE = {
    "pytest": "pytest",
    "pillow": "PIL",
    "scikit-learn": "sklearn",
    "python-dateutil": "dateutil",
    "pyyaml": "yaml",
}

# O que cada opcional DESLIGA se faltar. Quem decide se baixa alguns GB precisa
# saber o que perde se nao baixar.
DESLIGA_SEM = {
    "torch": "a evidencia PREVISTA (PP4CKD): a etapa 1 roda em modo degradado, "
             "so com evidencia medida, e a fase 3 nao calcula o canal "
             "negativo -- a coluna sai `indeterminado`, nunca 0",
    "admet-ai": "as 104 colunas de ADMET da etapa 2; os criterios que dependem "
                "delas saem `indeterminado`",
}


def _modulo(pacote: str) -> str:
    return MODULO_DO_PACOTE.get(pacote.lower(), pacote.replace("-", "_"))


def _ler(caminho: Path) -> list[dict]:
    if not caminho.is_file():
        return []
    itens = []
    for linha in caminho.read_text(encoding="utf-8").splitlines():
        limpa = linha.split("#", 1)[0].strip()
        if not limpa:
            continue
        pacote = limpa
        for separador in ("<=", ">=", "==", "~=", "!=", "<", ">", "["):
            pacote = pacote.split(separador, 1)[0]
        pacote = pacote.strip()
        if not pacote:
            continue
        itens.append({"pacote": pacote, "requisito": limpa,
                      "modulo": _modulo(pacote)})
    return itens


def exigidas() -> list[dict]:
    """O que `requirements.txt` declara. Sem isto o pipeline nao roda."""
    return _ler(caminhos.RAIZ / "requirements.txt")


def opcionais() -> list[dict]:
    """O que `requirements-opcionais.txt` declara, com o que cada um desliga.

    Nao sao instaladas automaticamente: ver o docstring do modulo."""
    itens = _ler(caminhos.RAIZ / "requirements-opcionais.txt")
    for item in itens:
        item["desliga"] = DESLIGA_SEM.get(
            item["pacote"].lower(),
            "parte do pipeline que depende dela; o codigo diz, na hora, o que "
            "ficou indeterminado")
    return itens


def _presente_de_verdade(modulo: str) -> bool:
    try:
        return importlib.util.find_spec(modulo) is not None
    except (ImportError, ValueError):
        return False


def conferir(*, presente=None) -> dict:
    """O estado do ambiente, sem mexer em nada.

    `completo` e False apenas quando falta algo EXIGIDO: opcional ausente nao
    e falta, e trata-la como falta pararia a execucao por algo que nao para."""
    presente = presente or _presente_de_verdade
    faltando = [d for d in exigidas() if not presente(d["modulo"])]
    ausentes = [o for o in opcionais() if not presente(o["modulo"])]
    return {"faltando": faltando, "opcionais_ausentes": ausentes,
            "completo": not faltando, "instalou": [], "erros": []}


def instalar_com_pip(requisitos) -> None:
    """Instala com o pip DO INTERPRETADOR CORRENTE.

    `sys.executable -m pip`, nunca o `pip` do PATH: o do PATH pode ser de outro
    ambiente, e instalar no ambiente errado parece sucesso e nao resolve nada.

    Sem `--upgrade`: a tarefa e suprir o que falta, nao mexer no que esta
    instalado e funcionando."""
    comando = [sys.executable, "-m", "pip", "install", *requisitos]
    r = subprocess.run(comando, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(
            (r.stderr or r.stdout or "pip falhou sem mensagem").strip()[-800:])


def garantir(*, instalar: bool = True, presente=None, instalador=None,
             avisar=None) -> dict:
    """Confere o ambiente e, se `instalar`, supre o que falta. Nao levanta.

    Devolve o estado com `instalou` (o que foi suprido) e `erros` (o que
    continua faltando, e por que). Quem chama decide o que fazer -- o pipeline
    para com codigo 2, o painel mostra na tela."""
    presente = presente or _presente_de_verdade
    instalador = instalador or instalar_com_pip
    avisar = avisar or (lambda msg: print(msg, file=sys.stderr))
    estado = conferir(presente=presente)
    if estado["completo"] or not instalar:
        return estado

    requisitos = [d["requisito"] for d in estado["faltando"]]
    avisar(f"[info] dependencias ausentes: "
           f"{', '.join(d['pacote'] for d in estado['faltando'])}. "
           f"Instalando com `{Path(sys.executable).name} -m pip install "
           f"{' '.join(requisitos)}`")
    try:
        instalador(requisitos)
    except Exception as erro:  # pip, rede, permissao: nada disso derruba
        estado["erros"].append(
            f"a instalacao falhou: {erro}. Resolva a mao com "
            f"`pip install -r requirements.txt`")
        return estado

    # Conferir DE NOVO: instalou e o pacote continua ausente e falha, nao
    # sucesso. Seguir em frente daria um erro de import no meio da etapa 1.
    importlib.invalidate_caches()
    restantes = [d for d in estado["faltando"] if not presente(d["modulo"])]
    estado["instalou"] = [d["pacote"] for d in estado["faltando"]
                          if d not in restantes]
    estado["faltando"] = restantes
    estado["completo"] = not restantes
    for item in restantes:
        estado["erros"].append(
            f"{item['pacote']} continua ausente depois da instalacao; "
            f"resolva a mao com `pip install -r requirements.txt`")
    if estado["instalou"]:
        avisar(f"[info] instalado: {', '.join(estado['instalou'])}")
    return estado


def relato(estado: dict) -> list[str]:
    """As linhas que o pipeline (ou o painel) mostra sobre o ambiente."""
    if estado["completo"] and not estado["opcionais_ausentes"]:
        return ["[info] dependencias: todas presentes"]
    linhas = []
    if estado["faltando"]:
        linhas.append(
            "[erro] faltam dependencias EXIGIDAS: "
            + ", ".join(d["requisito"] for d in estado["faltando"])
            + ". Instale com `pip install -r requirements.txt`")
    for erro in estado["erros"]:
        linhas.append(f"[erro] {erro}")
    if estado["instalou"]:
        linhas.append(f"[info] instalado agora: "
                      f"{', '.join(estado['instalou'])}")
    for item in estado["opcionais_ausentes"]:
        linhas.append(
            f"[aviso] `{item['pacote']}` nao esta instalado (opcional, "
            f"pesado). Sem ele, desliga {item['desliga']}. Para instalar: "
            f"`pip install -r requirements-opcionais.txt`")
    if estado["completo"] and not linhas:
        linhas.append("[info] dependencias: todas presentes")
    return linhas
