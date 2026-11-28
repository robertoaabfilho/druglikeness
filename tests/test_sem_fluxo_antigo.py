# -*- coding: utf-8 -*-
"""A guarda do anti-objetivo.

A exigencia nao e "nao importar o codigo antigo": e **nao deixar referencia
nenhuma a ele**. Um import, um nome de modulo herdado ou uma frase do tipo
"substituido por" falham estes testes, porque todos convidam o leitor a
procurar um fluxo que nao existe mais nesta arvore.

E tambem: nada que RODE aqui pode apontar para outro projeto. O que vem de
fora entra por copia, uma vez, e fica registrado em `input/procedencia.json`.

Cada guarda e provada contra um caso sintetico que ela TEM de pegar -- um
detector que nunca morde e indistinguivel de um teste vazio.
"""
from pathlib import Path

import pytest

RAIZ = Path(__file__).resolve().parent.parent

# Os nomes do fluxo antigo: modulos, arquivos de configuracao, flags e as
# frases que apontam para ele.
NOMES_DO_FLUXO_ANTIGO = (
    "prever_alvos",
    "importar_prospector",
    "farmaco_completo",
    "ranking_multibase",
    "validar_enriquecimento",
    "checar_alertas_adversos",
    "comparar_metodos_ranking",
    "ranking_embeddings",
    "calibrar_deeppk",
    "validar_alvos_medidos",
    "config.txt",
    "rodar_alertas_adversos",
    "rodar_comparacao_metodos",
    "rodar_validacao_enriquecimento",
    "Deep-PK",
    "deeppk",
    "substituido pela",
    "substituido por",
    "substituída pela",
    "fluxo antigo",
    "fluxo original",
)

# Caminhos de OUTROS projetos, em forma de CAMINHO -- nao a palavra solta.
# "prospector" como palavra e legitimo e as vezes obrigatorio: a licenca do
# dado derivado tem de dizer de onde ele veio. O que nao pode e um caminho
# para a arvore de outro projeto, ou um caminho absoluto na casa de alguem.
OUTROS_PROJETOS = ("pipeline_prospeccao", "modelos_baixados", "/home/")

# Arquivos onde a procedencia E o assunto, logo citar a origem e correto.
ONDE_PROCEDENCIA_E_O_ASSUNTO = (
    "mods/instantaneo.py",
    "tests/test_instantaneo.py",
    "tests/test_sem_fluxo_antigo.py",
    "tests/test_caminhos.py",
    "README.md",
    "input/procedencia.json",
    # Atribuicao de licenca: dizer de qual arvore o dado derivou e OBRIGATORIO.
    "dados/LICENSE",
    "dados/validacao_neuro_prospector.LEIA-ME.md",
    "tests/dados/LICENSE",
)


def _fontes():
    """Todo arquivo de texto versionado da arvore, menos o que o git ignora."""
    for caminho in sorted(RAIZ.rglob("*")):
        if not caminho.is_file():
            continue
        rel = caminho.relative_to(RAIZ).as_posix()
        if rel.startswith((".git/", "tsv/", "modelo/", "json/cache_rede/",
                           ".superpowers/")):
            continue
        if "__pycache__" in rel or ".pytest_cache" in rel:
            continue
        if caminho.suffix not in (".py", ".md", ".json", ".txt", ".tsv", ".cfg",
                                  ".toml", ".ini", ""):
            continue
        try:
            yield rel, caminho.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue


def _mencoes(texto: str, nomes) -> list[str]:
    return [n for n in nomes if n in texto]


def test_a_guarda_do_fluxo_antigo_morde():
    """Prova que o detector nao passa por vacuidade."""
    assert _mencoes("from mods.prever_alvos import carregar_modelo",
                    NOMES_DO_FLUXO_ANTIGO) == ["prever_alvos"]
    assert _mencoes("isto foi substituido pela fase 3",
                    NOMES_DO_FLUXO_ANTIGO) == ["substituido pela"]
    assert _mencoes("texto inocente", NOMES_DO_FLUXO_ANTIGO) == []


def test_nenhum_arquivo_da_arvore_menciona_o_fluxo_antigo():
    achados = []
    for rel, texto in _fontes():
        if rel == "tests/test_sem_fluxo_antigo.py":
            continue  # a lista de nomes proibidos mora aqui
        for nome in _mencoes(texto, NOMES_DO_FLUXO_ANTIGO):
            achados.append(f"{rel}: {nome!r}")
    assert not achados, "referencias ao fluxo antigo:\n  " + "\n  ".join(achados)


def test_a_guarda_de_outro_projeto_morde():
    assert _mencoes("~/prospector/pipeline_prospeccao_v4/resultados/x.csv",
                    OUTROS_PROJETOS) == ["pipeline_prospeccao"]
    assert _mencoes("/home/alguem/modelos_baixados/pp4ckd",
                    OUTROS_PROJETOS) == ["modelos_baixados", "/home/"]
    assert _mencoes("modelo/ecfp4_fold_4.pt", OUTROS_PROJETOS) == []


def test_nenhum_codigo_que_roda_aponta_para_outro_projeto():
    achados = []
    for rel, texto in _fontes():
        if rel in ONDE_PROCEDENCIA_E_O_ASSUNTO:
            continue
        for nome in _mencoes(texto, OUTROS_PROJETOS):
            achados.append(f"{rel}: {nome!r}")
    assert not achados, (
        "caminho de outro projeto fora de onde a procedencia e o assunto:\n  "
        + "\n  ".join(achados))


@pytest.mark.parametrize("modulo", [
    "mods.caminhos", "mods.quimica", "mods.modelo", "mods.metricas",
])
def test_os_modulos_de_infraestrutura_importam_sozinhos(modulo):
    """Autossuficiencia: nada da arvore nova depende da antiga para carregar."""
    __import__(modulo)


def test_input_versiona_a_procedencia_e_nao_os_dados():
    """`input/` guarda instantaneos de 79 MB e 60 MB. O que o git carrega e o
    REGISTRO (origem, data, MD5, contagens), nao os dados -- a mesma regra de
    `modelo/`. Versionar os dois seria 140 MB num repositorio de codigo, e
    versionar um e nao o outro nao teria criterio."""
    import subprocess
    r = subprocess.run(["git", "ls-files", "input/"], cwd=RAIZ,
                       capture_output=True, text=True)
    if r.returncode != 0:
        # Fora de um checkout (o pacote para outra maquina) nao ha indice do
        # git a consultar, e nao ha propriedade equivalente: esta regra e
        # sobre o REPOSITORIO.
        #
        # Eu tinha inventado uma substituta -- "a pasta input/ nao traz dados
        # grandes" -- e ela ficou FALSA quando o pacote passou a levar tudo
        # para garantir reprodutibilidade. Uma propriedade substituta que nao
        # e a mesma propriedade mente duas vezes: aqui, e sobre o pacote.
        import pytest
        pytest.skip("fora de um checkout de git: esta regra e do repositorio")
    assert r.stdout.split() == ["input/procedencia.json"], r.stdout


def test_a_procedencia_nao_carrega_caminho_absoluto_de_destino():
    """`procedencia.json` e versionado: um caminho absoluto ali vaza a pasta
    pessoal de quem rodou e nao resolve na maquina de ninguem mais. O destino
    e relativo a raiz; a ORIGEM continua absoluta, porque e justamente o que
    se quer registrar."""
    import json
    caminho = RAIZ / "input" / "procedencia.json"
    if not caminho.is_file():
        pytest.skip("sem instantaneo nesta arvore")
    dados = json.loads(caminho.read_text(encoding="utf-8"))
    for entrada in dados["entradas"]:
        assert not entrada.get("destino", "").startswith("/"), entrada["papel"]
