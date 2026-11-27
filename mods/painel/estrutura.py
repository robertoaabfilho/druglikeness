# -*- coding: utf-8 -*-
"""
estrutura.py - a arquitetura do pipeline, lida do proprio codigo.

A barra lateral do painel mostra o que o pipeline E: as etapas, as faixas de
evidencia, as pastas, os codigos de saida. **Nada disso e texto decorado
aqui.** As etapas vem da tabela `ETAPAS` do `main.py`, as faixas de
`mods/ranque/faixa.py`, as pastas de `mods/caminhos.py`, as dependencias do
`requirements.txt`.

Por que isso importa: uma descricao copiada divergiria do pipeline na primeira
mudanca, e documentacao que mente e pior que nenhuma. Aqui, renomear uma saida
ou acrescentar uma faixa muda a tela sozinho.

O unico texto proprio deste modulo e o PAPEL de cada pasta e o significado de
cada codigo de saida -- os dois sao decisoes de desenho que o codigo nao
carrega em forma legivel, e os dois tem teste que exige que nenhum fique
vazio.
"""
from pathlib import Path

from mods import caminhos
from mods.ranque.esquema import MOTIVO_FAIXA_3
from mods.ranque.faixa import FAIXAS

# Papel de cada pasta. A regra que desempata `input/` de `tsv/`, porque
# ligantes e alvos tambem sao TSV: PROCEDENCIA, nao extensao.
PAPEL_DAS_PASTAS = {
    "input/": "o que VOCE fornece: ligantes, alvos, sintomas, evidencia "
              "medida. Fora do git (so a procedencia e versionada)",
    "mods/": "o codigo: as tres CLIs, a infraestrutura e a biblioteca do "
             "painel",
    "modelo/": "os pesos do preditor e os rotulos dos alvos. Fora do git",
    "json/": "config.json (o parametro de cada etapa) e o cache de rede, uma "
             "subpasta por versao dos servicos",
    "tsv/": "o que a MAQUINA gerou, uma pasta por etapa, cada tabela com o "
            "seu manifesto.json ao lado. A regra que desempata com input/ e "
            "procedencia, nao extensao",
    "dados/": "dado de referencia versionado, com licenca propria",
    "exemplos/": "listas de exemplo para comecar: copie para input/",
    "tests/": "a suite; nenhum teste usa rede",
}

# O que cada codigo de saida significa. O orquestrador devolve o PIOR entre as
# etapas que rodou, e o 3 existe porque automacao nao le o stderr.
CODIGOS_DE_SAIDA = {
    0: "execucao completa",
    3: "parcial: a saida existe, mas linhas de entrada foram descartadas",
    1: "nenhuma entrada utilizavel (nada e escrito)",
    2: "erro (nada e publicado; as saidas anteriores ficam intactas)",
}

# Significado de cada faixa, na linguagem de quem vai citar o resultado.
SIGNIFICADO_DAS_FAIXAS = {
    "1a": "evidencia MEDIDA (ChEMBL) num alvo que VOCE declarou. Ordena por "
          "pchembl",
    "1b": "evidencia MEDIDA num alvo que o enriquecimento acrescentou. Ordena "
          "por pchembl",
    "2a": "evidencia PREVISTA pelo modelo num alvo que VOCE declarou. Ordena "
          "por probabilidade",
    "2b": "evidencia PREVISTA num alvo que o enriquecimento acrescentou",
    "3": f"{MOTIVO_FAIXA_3}: NAO RANQUEIA, e aparece em nao_ranqueados.tsv "
         "com o motivo",
}


def etapas() -> list[dict]:
    """As etapas, na ordem em que rodam, como o orquestrador as declara."""
    import main
    return [dict(e) for e in main.ETAPAS]


def faixas() -> list[dict]:
    """As faixas de evidencia, na ordem da hierarquia.

    `ranqueia` e False so para a faixa 3: entrar na lista nunca foi entrar no
    ranking."""
    return [{"faixa": f, "significado": SIGNIFICADO_DAS_FAIXAS.get(f, ""),
             "ranqueia": f != "3"} for f in FAIXAS]


def pastas() -> list[dict]:
    """As pastas da arvore, com o papel de cada uma."""
    return [{"nome": nome, "papel": papel}
            for nome, papel in PAPEL_DAS_PASTAS.items()]


def codigos_de_saida() -> dict:
    return dict(CODIGOS_DE_SAIDA)


def dependencias() -> list[str]:
    """Os pacotes do `requirements.txt`, sem comentario nem versao."""
    caminho = caminhos.RAIZ / "requirements.txt"
    if not caminho.is_file():
        return []
    nomes = []
    for linha in caminho.read_text(encoding="utf-8").splitlines():
        limpa = linha.strip()
        if not limpa or limpa.startswith("#"):
            continue
        nomes.append(limpa.split("<")[0].split(">")[0].split("=")[0].strip())
    return nomes


def arvore() -> str:
    """A arvore de pastas como texto, para mostrar numa linha de codigo."""
    linhas = ["main.py        orquestrador das tres etapas",
              "painel.py      esta interface"]
    for pasta in pastas():
        linhas.append(f"{pasta['nome']:<15}{pasta['papel'].split(':')[0]}")
    return "\n".join(linhas)


def onde_cada_tabela_vai() -> list[dict]:
    """Tabela e manifesto de cada etapa, pelos caminhos reais."""
    saida = []
    for etapa in etapas():
        saida.append({"etapa": etapa["id"],
                      "tabela": etapa["saida"],
                      "manifesto": etapa["manifesto"]})
    return saida


def versao_do_pipeline() -> str:
    """O commit do repositorio, quando ha um. `desconhecido` fora de um
    checkout -- ausencia nao vira veredito."""
    import subprocess
    try:
        r = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                           cwd=str(caminhos.RAIZ), capture_output=True,
                           text=True, timeout=5)
    except (OSError, subprocess.SubprocessError):
        return "desconhecido"
    return r.stdout.strip() if r.returncode == 0 else "desconhecido"
