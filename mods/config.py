# -*- coding: utf-8 -*-
"""
config.py - o parametro de cada etapa, num lugar so.

`json/config.json` tem uma secao por etapa, e cada chave e o nome da flag da
CLI daquela etapa com `-` trocado por `_`. O orquestrador le este arquivo e
monta a linha de comando; o pesquisador edita um arquivo, nao tres comandos
com 46 flags somadas.

**Nome que nao existe e ERRO.** Secao desconhecida, chave desconhecida e valor
de tipo errado param a execucao com codigo 2, nomeando o que esta errado. E a
mesma regra das flags (`--reprovar criterio_inexistente` tambem e erro), e
existe pelo mesmo motivo: uma politica que o pesquisador acredita estar valendo
e nao esta e pior que uma mensagem de erro.

Uma chave que comece com `_` e COMENTARIO, nao parametro: e assim que o
arquivo carrega, ao lado de cada numero, a medicao que o escolheu. JSON nao
tem comentario, e um numero sem o motivo convida a ser mexido.

**Chave ausente nao vira valor.** O que nao esta no arquivo simplesmente nao
entra na linha de comando, e a CLI aplica o proprio padrao -- que e onde os
padroes moram e estao documentados. Isto evita a segunda fonte de verdade que
fez o arquivo de configuracao antigo acumular vinte chaves.
"""
import json
from difflib import get_close_matches
from pathlib import Path

from mods import caminhos
from mods.evidencias import atomico


class ErroDeConfig(Exception):
    """Configuracao invalida. Quem captura devolve codigo 2."""


# Tipos aceitos por chave, por etapa. A chave e o nome da flag com `_`.
# `tests/test_config.py` prende cada uma destas contra a flag real do script:
# o esquema nao pode descolar das CLIs sem o teste gritar.
CHAVES: dict[str, dict[str, tuple]] = {
    "etapa1": {
        "ligantes": (str,), "alvos": (str,), "sintomas": (str,),
        "saida": (str,), "checkpoint": (str,), "rotulos": (str,),
        "evidencias_chembl": (str,), "alvos_enriquecidos": (str,),
        "cache": (str,), "manifesto": (str,),
        "corte": (int, float), "corte_associacao": (int, float),
        "teto_via": (int,), "limite": (int,),
        "efo": (list,), "sem_rede": (bool,),
    },
    "portao": {
        "evidencias": (str,), "ligantes": (str,), "alvos": (str,),
        "saida": (str,), "cache": (str,), "manifesto": (str,),
        "termos_alerta": (str, list), "tecidos": (str, list),
        "reprovar": (str, list), "marcar": (str, list),
        "corte_expressao": (int, float), "percentil_toxicidade": (int, float),
        "limite": (int,), "limite_bula": (int,), "workers": (int,),
        "gpu": (bool,), "sem_rede": (bool,),
    },
    "ranking": {
        "evidencias": (str,), "ligantes": (str,), "alvos": (str,),
        "portao": (str,), "saida": (str,), "checkpoint": (str,),
        "rotulos": (str,), "manifesto": (str,),
        "manifesto_etapa1": (str,),
        "limite": (int,), "lote": (int,), "minimo_treino": (int,),
        "gpu": (bool,), "sem_canal_negativo": (bool,),
    },
}


def _conferir_secao(etapa: str, valores: dict) -> None:
    conhecidas = CHAVES[etapa]
    for chave, valor in valores.items():
        if chave.startswith("_"):
            continue  # comentario: "_nota_x": "..." documenta a chave ao lado
        if chave not in conhecidas:
            parecida = get_close_matches(chave, conhecidas, n=1, cutoff=0.6)
            dica = f"; quis dizer {parecida[0]!r}?" if parecida else ""
            raise ErroDeConfig(
                f"{etapa}: chave desconhecida {chave!r}{dica}. "
                f"As chaves desta secao sao: "
                f"{', '.join(sorted(conhecidas))}")
        if valor is None:
            continue
        tipos = conhecidas[chave]
        # bool e subclasse de int em Python: um booleano num campo numerico
        # passaria em silencio, e `limite: true` nao e um limite.
        if isinstance(valor, bool) and bool not in tipos:
            raise ErroDeConfig(
                f"{etapa}.{chave}: esperado "
                f"{' ou '.join(t.__name__ for t in tipos)}, "
                f"veio booleano {valor!r}")
        if not isinstance(valor, tipos):
            raise ErroDeConfig(
                f"{etapa}.{chave}: esperado "
                f"{' ou '.join(t.__name__ for t in tipos)}, "
                f"veio {valor!r} ({type(valor).__name__})")


def ler(caminho=None) -> dict:
    """Le e valida a configuracao. Devolve {etapa: {chave: valor}}."""
    caminho = Path(caminho) if caminho else caminhos.CONFIG
    if not caminho.is_file():
        raise ErroDeConfig(f"configuracao nao encontrada: {caminho}")
    try:
        dados = json.loads(caminho.read_text(encoding="utf-8"))
    except json.JSONDecodeError as e:
        raise ErroDeConfig(f"{caminho}: JSON invalido -- {e}") from e
    if not isinstance(dados, dict):
        raise ErroDeConfig(f"{caminho}: a raiz tem de ser um objeto")
    for etapa, valores in dados.items():
        if etapa.startswith("_"):
            continue  # comentario: "_nota": "..."
        if etapa not in CHAVES:
            raise ErroDeConfig(
                f"secao desconhecida {etapa!r}. As secoes sao: "
                f"{', '.join(CHAVES)}")
        if not isinstance(valores, dict):
            raise ErroDeConfig(f"{etapa}: esperado objeto, veio {valores!r}")
        _conferir_secao(etapa, valores)
    return {e: v for e, v in dados.items() if not e.startswith("_")}


def gravar(valores: dict, caminho=None) -> Path:
    """Escreve `valores` sobre a configuracao existente e devolve o caminho.

    E por aqui que o painel grava o formulario. Tres garantias, nessa ordem:

    **As notas sobrevivem.** Cada numero do arquivo tem um `_nota_*` ao lado
    carregando a medicao que o escolheu (`corte` 0.1 e o regime de producao,
    `corte_associacao` 0.4 contra 0.5 sao 44 contra 30 alvos). Um `json.dump`
    do dicionario de valores apagaria todas, e o proximo a mexer no numero
    nao saberia por que ele e aquele. O que se grava e a MESCLA sobre o
    arquivo atual, nunca uma substituicao.

    **A ordem sobrevive.** A ordem do arquivo e a ordem em que o pesquisador
    le, e cada nota esta ao lado da sua chave; reordenar por alfabeto
    separaria as duas.

    **Valida antes de escrever.** Chave, secao ou tipo errado nao chegam ao
    disco: a configuracao e de onde as tres etapas leem, e um arquivo pela
    metade pararia todas.

    Chave com valor `None` e REMOVIDA, nao gravada como `null`: apagar o campo
    no formulario devolve o parametro ao padrao da CLI, que e onde os padroes
    moram e estao documentados."""
    caminho = Path(caminho) if caminho else caminhos.CONFIG
    ler(caminho)  # valida o que ja esta no disco, e levanta se nao existe
    bruto = json.loads(caminho.read_text(encoding="utf-8"))  # com os `_`
    for etapa, mudancas in valores.items():
        if etapa not in CHAVES:
            raise ErroDeConfig(
                f"secao desconhecida {etapa!r}. As secoes sao: "
                f"{', '.join(CHAVES)}")
        _conferir_secao(etapa, mudancas)
    for etapa, mudancas in valores.items():
        secao = bruto.setdefault(etapa, {})
        for chave, valor in mudancas.items():
            if valor is None:
                secao.pop(chave, None)
            else:
                secao[chave] = valor
    atomico.escrever_texto(
        caminho, json.dumps(bruto, ensure_ascii=False, indent=2) + "\n")
    return caminho


def argumentos_de(etapa: str, cfg: dict) -> list[str]:
    """As flags da CLI daquela etapa, na ordem das chaves.

    Chave ausente ou nula NAO entra: a CLI aplica o proprio padrao. Booleano
    falso tambem nao entra -- `--gpu: false` nao e `--gpu`.
    """
    if etapa not in CHAVES:
        raise ErroDeConfig(f"secao desconhecida {etapa!r}")
    args: list[str] = []
    for chave, valor in (cfg.get(etapa) or {}).items():
        if chave.startswith("_"):
            continue  # comentario, nao parametro
        flag = "--" + chave.replace("_", "-")
        if valor is None:
            continue
        if isinstance(valor, bool):
            if valor:
                args.append(flag)
            continue
        if isinstance(valor, list):
            for item in valor:
                args += [flag, str(item)]
            continue
        args += [flag, str(valor)]
    return args
