#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
empacotar.py - o pacote para rodar o pipeline noutra maquina.

A pergunta que define o conteudo: **o que falta numa maquina limpa?** O codigo
vem do git; o que o git NAO carrega -- os pesos do modelo, os rotulos dos
alvos, as entradas -- e exatamente o que o pacote precisa levar.

Tres regras
-----------
**O que entra vem com a procedencia.** `input/procedencia.json` vai sempre:
sem ele o resultado na outra maquina nao e citavel, porque nao se sabe de que
instantaneo ele veio.

**O que NAO cabe e DITO, nao omitido em silencio.** A base de ligantes tem
75 MB e a evidencia medida do ChEMBL outros 60; por padrao eles ficam fora, e
o manifesto do pacote registra o que ficou fora, o tamanho e o comando que o
traz. Um pacote que omite calado faz o pesquisador descobrir a falta do outro
lado.

**As saidas nao vao.** `tsv/` e resultado de uma execucao DESTA maquina;
levar junto faria a outra parecer ter rodado o pipeline. O cache de rede
tambem fica: ele e indexado por versao de servico, e servir cache velho noutra
maquina e pior que consultar de novo.

Uso
---
    python mods/empacotar.py --saida ../druglikeness-pacote.zip
    python mods/empacotar.py --saida ... --com-dados   # inclui os 135 MB
"""
import argparse
import json
import subprocess
import sys
import time
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from mods import caminhos  # noqa: E402

SAIDA_OK, SAIDA_ERRO = 0, 2

# Pastas e arquivos que o pacote leva. O codigo tambem vai, mesmo vindo do
# git: quem recebe o pacote pode nao ter como clonar.
RAIZES = ("mods", "tests", "exemplos", "dados", "json")
ARQUIVOS = ("main.py", "painel.py", "requirements.txt",
            "requirements-opcionais.txt", "README.md", "LICENSE")

# O que nunca entra, e por que.
EXCLUIR = {
    "__pycache__": "bytecode desta maquina e da versao de Python dela",
    ".pytest_cache": "cache de execucao dos testes",
    ".git": "o repositorio: quem recebe o pacote nao precisa do historico",
    ".superpowers": "area de trabalho de desenvolvimento",
    "cache_rede": "cache indexado por versao de servico; servir cache velho "
                  "noutra maquina e pior que consultar de novo",
}

# Pastas de dado GERADO que ficam fora do pacote, com o comando que as traz. O
# acervo do openFDA sao 3,7 GB (1,7 de particoes e 2,0 de indice): levar junto
# faria o pacote passar de 31 MB para mais de 3 GB, e a outra maquina baixa o
# seu com um comando. A mesma regra de `tsv/` e do cache.
GERADO_FORA = {
    "dados/openfda": "o acervo do openFDA (3,7 GB com o indice); baixe o seu "
                     "com `python mods/openfda_massa.py --baixar --indexar`",
}

# Dados grandes: ficam fora por padrao, e o manifesto diz como trazer.
DADOS_GRANDES = {
    "input/ligantes.tsv":
        "a base de ligantes; traga com `python mods/instantaneo.py "
        "--ligantes-origem <arquivo>`",
    "input/evidencias_chembl.csv":
        "a bioatividade medida do ChEMBL; traga com `--evidencias-origem`",
}


def _excluido(caminho: Path) -> str | None:
    como_texto = caminho.as_posix()
    for prefixo, motivo in GERADO_FORA.items():
        if como_texto.startswith(prefixo + "/"):
            # a procedencia fica: a outra maquina precisa saber de que export
            # era o acervo desta, como em `input/procedencia.json`
            if caminho.name != "procedencia.json":
                return motivo
    for parte in caminho.parts:
        if parte in EXCLUIR:
            return EXCLUIR[parte]
    if caminho.suffix in (".pyc", ".pyo"):
        return EXCLUIR["__pycache__"]
    return None


def conteudo(*, incluir_dados_grandes: bool = False) -> list[dict]:
    """Os arquivos do pacote: `{origem, destino, bytes}`."""
    itens: list[dict] = []
    vistos: set[str] = set()

    def juntar(origem: Path, destino: str):
        if destino in vistos or not origem.is_file():
            return
        vistos.add(destino)
        itens.append({"origem": origem, "destino": destino,
                      "bytes": origem.stat().st_size})

    for nome in ARQUIVOS:
        juntar(caminhos.RAIZ / nome, nome)
    for raiz in RAIZES:
        pasta = caminhos.RAIZ / raiz
        if not pasta.is_dir():
            continue
        for arquivo in sorted(pasta.rglob("*")):
            if not arquivo.is_file():
                continue
            relativo = arquivo.relative_to(caminhos.RAIZ)
            if _excluido(relativo):
                continue
            juntar(arquivo, relativo.as_posix())

    # O modelo: e o que o git nao carrega e sem o que a etapa 1 degrada.
    if caminhos.MODELO.is_dir():
        for arquivo in sorted(caminhos.MODELO.iterdir()):
            if arquivo.is_file():
                juntar(arquivo, f"modelo/{arquivo.name}")

    # A procedencia vai SEMPRE: sem ela o resultado nao e citavel.
    juntar(caminhos.PROCEDENCIA, "input/procedencia.json")
    for alvo in ("alvos.txt", "sintomas.txt"):
        juntar(caminhos.INPUT / alvo, f"input/{alvo}")
    if incluir_dados_grandes:
        for destino in DADOS_GRANDES:
            juntar(caminhos.RAIZ / destino, destino)
    return itens


def omitidos(*, incluir_dados_grandes: bool = False) -> list[dict]:
    """O que existe nesta maquina e NAO foi para o pacote, com o motivo."""
    fora = []
    if not incluir_dados_grandes:
        for destino, motivo in DADOS_GRANDES.items():
            caminho = caminhos.RAIZ / destino
            if caminho.is_file():
                fora.append({"arquivo": destino, "motivo": motivo,
                             "mb": round(caminho.stat().st_size / (1 << 20), 1)})
    if caminhos.TSV.is_dir() and any(caminhos.TSV.iterdir()):
        fora.append({
            "arquivo": "tsv/",
            "motivo": "as saidas sao resultado de uma execucao DESTA maquina; "
                      "levar junto faria a outra parecer ter rodado o pipeline",
            "mb": 0.0})
    if caminhos.CACHE.is_dir() and any(caminhos.CACHE.iterdir()):
        fora.append({"arquivo": "json/cache_rede/",
                     "motivo": EXCLUIR["cache_rede"], "mb": 0.0})
    for prefixo, motivo in GERADO_FORA.items():
        pasta = caminhos.RAIZ / prefixo
        if not pasta.is_dir():
            continue
        tamanho = sum(a.stat().st_size for a in pasta.rglob("*")
                      if a.is_file() and a.name != "procedencia.json")
        if tamanho:
            fora.append({"arquivo": f"{prefixo}/", "motivo": motivo,
                         "mb": round(tamanho / (1 << 20), 1)})
    return fora


def _commit() -> str:
    try:
        r = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                           cwd=str(caminhos.RAIZ), capture_output=True,
                           text=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        return "desconhecido"
    return r.stdout.strip() if r.returncode == 0 else "desconhecido"


def _leia_me(itens: list[dict], fora: list[dict], com_modelo: bool) -> str:
    linhas = [
        "# Pipeline de drogabilidade — pacote para rodar noutra máquina",
        "",
        f"Gerado em {time.strftime('%Y-%m-%d %H:%M')}, do commit "
        f"`{_commit()}`. {len(itens)} arquivos.",
        "",
        "## Instalar",
        "",
        "```bash",
        "python -m venv .venv",
        "source .venv/bin/activate      # Windows: .venv\\Scripts\\activate",
        "pip install -r requirements.txt",
        "python -m pytest -q            # nenhum teste usa rede",
        "```",
        "",
        "A suíte verde é a prova de que o pacote chegou inteiro. Se ela falhar,",
        "não rode o pipeline antes de entender por quê.",
        "",
        "## Rodar",
        "",
        "```bash",
        "python main.py                 # as três etapas",
        "streamlit run painel.py        # a interface",
        "```",
        "",
    ]
    if com_modelo:
        linhas += ["Os pesos do modelo vieram no pacote, em `modelo/`.", ""]
    else:
        linhas += [
            "**O modelo NÃO veio no pacote.** Sem ele a etapa 1 degrada: roda só",
            "com a evidência medida e diz o que falta. Traga com",
            "`python mods/instantaneo.py --modelo-origem <pasta>`.", ""]
    if fora:
        linhas += ["## O que ficou fora, e por quê", ""]
        for item in fora:
            tamanho = f" ({item['mb']:.0f} MB)" if item["mb"] else ""
            linhas.append(f"- **`{item['arquivo']}`**{tamanho} — "
                          f"{item['motivo']}")
        linhas.append("")
    linhas += [
        "## O que é seu e o pipeline não deduz",
        "",
        "`input/alvos.txt` e `input/sintomas.txt` são suas listas. Se não",
        "vieram, copie de `exemplos/` e edite — o pipeline não as adivinha.",
        "",
        "Leia `README.md` para a arquitetura, as etapas e as limitações",
        "medidas. As limitações importam: elas dizem o que o resultado **não**",
        "sustenta.",
        "",
    ]
    return "\n".join(linhas)


def montar(destino, *, incluir_dados_grandes: bool = False,
           raiz_no_zip: str = "druglikeness") -> dict:
    """Escreve o zip e devolve o relatorio.

    Tudo vai dentro de UMA pasta raiz: descompactar nao pode espalhar
    arquivos no diretorio corrente."""
    destino = Path(destino)
    itens = conteudo(incluir_dados_grandes=incluir_dados_grandes)
    fora = omitidos(incluir_dados_grandes=incluir_dados_grandes)
    com_modelo = any(i["destino"].startswith("modelo/") for i in itens)
    manifesto = {
        "gerado_em": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "commit": _commit(),
        "arquivos": len(itens),
        "bytes": sum(i["bytes"] for i in itens),
        "com_modelo": com_modelo,
        "com_dados_grandes": bool(incluir_dados_grandes),
        "omitidos": fora,
    }
    destino.parent.mkdir(parents=True, exist_ok=True)
    temporario = destino.with_name(destino.name + ".parcial")
    try:
        with zipfile.ZipFile(temporario, "w", zipfile.ZIP_DEFLATED,
                             compresslevel=6) as z:
            for item in itens:
                z.write(item["origem"], f"{raiz_no_zip}/{item['destino']}")
            z.writestr(f"{raiz_no_zip}/PACOTE.json",
                       json.dumps(manifesto, ensure_ascii=False, indent=2,
                                  default=str) + "\n")
            z.writestr(f"{raiz_no_zip}/LEIA-ME-PRIMEIRO.md",
                       _leia_me(itens, fora, com_modelo))
        temporario.replace(destino)
    finally:
        temporario.unlink(missing_ok=True)
    manifesto["zip_bytes"] = destino.stat().st_size
    manifesto["destino"] = str(destino)
    return manifesto


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--saida", required=True, metavar="ARQUIVO.zip")
    ap.add_argument("--com-dados", action="store_true",
                    help="inclui a base de ligantes e a evidencia medida "
                         "(~135 MB); por padrao ficam fora e o manifesto diz "
                         "como traze-los")
    args = ap.parse_args(argv)
    try:
        rel = montar(args.saida, incluir_dados_grandes=args.com_dados)
    except OSError as erro:
        print(f"[erro] {erro}", file=sys.stderr)
        return SAIDA_ERRO
    print(f"[info] {rel['destino']}: {rel['arquivos']} arquivos, "
          f"{rel['zip_bytes'] / (1 << 20):.1f} MB compactados "
          f"(commit {rel['commit']})", file=sys.stderr)
    for item in rel["omitidos"]:
        print(f"[info] fora do pacote: {item['arquivo']} -- {item['motivo']}",
              file=sys.stderr)
    return SAIDA_OK


if __name__ == "__main__":
    sys.exit(main())
