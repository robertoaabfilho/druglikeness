# -*- coding: utf-8 -*-
"""
entrada.py

Leitura e validacao das tres entradas da etapa 1: ligantes, alvos e sintomas.

Regras
------
* TSV e o formato. `.csv` e aceito com aviso: 42,8% dos nomes de produto
  natural tem virgula, entao uma virgula nao citada quebra a linha.
* ligantes: colunas minimas `nome` e `smiles`; o resto e recalculado.
* alvos: `alvo_chembl_id`, OU `alvo` + `organismo`. Nome sem organismo e
  erro (nao aviso): "noradrenaline transporter" existe em humano, camundongo
  e rato, e aceitar o nome puro traria o ortologo errado sem ninguem notar.
* Ligantes sao deduplicados por `inchikey_parental`: sal e base livre sao o
  mesmo composto e duplicariam toda linha a jusante.
* Linha com erro e descartada e relatada; nada falha em silencio.
* Se o alvo existe ou nao no espaco do modelo PP4CKD NAO e decidido aqui
  (vem do arquivo de rotulos do modelo, onde for preciso).
"""

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from mods.importar_prospector import forma_parental, inchikey_de, molecula, texto


@dataclass
class Problema:
    linha: int
    campo: str
    mensagem: str
    gravidade: Literal["erro", "aviso"]


def _ler_tabela(caminho) -> tuple[list[str], list[tuple[int, dict]], list[Problema]]:
    """Le TSV (ou CSV, com aviso). Devolve (cabecalho, [(linha, dict)], problemas).
    Linhas com numero de celulas diferente do cabecalho sao descartadas."""
    caminho = Path(caminho)
    problemas: list[Problema] = []
    csv_ = caminho.suffix.lower() == ".csv"
    if csv_:
        problemas.append(Problema(
            1, "arquivo", "arquivo .csv aceito, mas o formato e TSV: nomes com "
            "virgula quebram a leitura; prefira .tsv", "aviso"))
    cabecalho: list[str] = []
    linhas: list[tuple[int, dict]] = []
    with open(caminho, newline="", encoding="utf-8-sig") as f:
        leitor = csv.reader(f, delimiter="," if csv_ else "\t")
        for celulas in leitor:
            if not any(c.strip() for c in celulas):
                continue
            if not cabecalho:
                cabecalho = [c.strip() for c in celulas]
                continue
            if len(celulas) != len(cabecalho):
                dica = " (virgula nao citada em .csv?)" if csv_ else ""
                problemas.append(Problema(
                    leitor.line_num, "linha", f"{len(celulas)} celulas, cabecalho "
                    f"tem {len(cabecalho)}{dica}; linha descartada", "erro"))
                continue
            linhas.append((leitor.line_num, dict(zip(cabecalho, celulas))))
    if not cabecalho:
        problemas.append(Problema(1, "arquivo", "arquivo vazio", "erro"))
    return cabecalho, linhas, problemas


def _exigir(cabecalho: list[str], colunas: list[str], problemas: list[Problema]) -> bool:
    faltam = [c for c in colunas if c not in cabecalho]
    for c in faltam:
        problemas.append(Problema(1, c, f"coluna obrigatoria '{c}' ausente", "erro"))
    return not faltam


def ler_ligantes(caminho) -> tuple[list[dict], list[Problema]]:
    cabecalho, linhas, problemas = _ler_tabela(caminho)
    if cabecalho and not _exigir(cabecalho, ["nome", "smiles"], problemas):
        return [], problemas
    ligantes: list[dict] = []
    vistos: dict[str, str] = {}  # inchikey_parental -> primeiro nome
    descartados: dict[str, list[str]] = {}
    for n, linha in linhas:
        nome, smiles = texto(linha.get("nome")), texto(linha.get("smiles"))
        if not nome:
            problemas.append(Problema(n, "nome", "nome ausente", "erro"))
            continue
        if not smiles:
            problemas.append(Problema(n, "smiles", f"{nome}: smiles ausente", "erro"))
            continue
        mol = molecula(smiles)
        parental = forma_parental(mol)
        chave_parental = inchikey_de(parental)
        if not chave_parental:
            problemas.append(Problema(
                n, "smiles", f"{nome}: smiles que o RDKit nao le: {smiles}", "erro"))
            continue
        registro = dict(linha)
        registro.update(nome=nome, smiles=smiles, inchikey=inchikey_de(mol),
                        inchikey_parental=chave_parental)
        informada = texto(linha.get("inchikey_parental"))
        if informada and informada != chave_parental:
            problemas.append(Problema(
                n, "inchikey_parental", f"{nome}: valor informado {informada} "
                f"difere do recalculado {chave_parental}; usado o recalculado",
                "aviso"))
        if chave_parental in vistos:
            descartados.setdefault(chave_parental, []).append(nome)
            continue
        vistos[chave_parental] = nome
        ligantes.append(registro)
    for chave, nomes in descartados.items():
        problemas.append(Problema(
            0, "nome", f"ligante duplicado pela forma parental ({chave}): mantido "
            f"'{vistos[chave]}', descartado(s): {', '.join(nomes)}", "aviso"))
    return ligantes, problemas


def ler_alvos(caminho) -> tuple[list[dict], list[Problema]]:
    cabecalho, linhas, problemas = _ler_tabela(caminho)
    if cabecalho and "alvo_chembl_id" not in cabecalho and not _exigir(
            cabecalho, ["alvo", "organismo"], problemas):
        return [], problemas
    alvos: list[dict] = []
    for n, linha in linhas:
        chembl_id = texto(linha.get("alvo_chembl_id"))
        nome, organismo = texto(linha.get("alvo")), texto(linha.get("organismo"))
        if not chembl_id:
            if not nome:
                problemas.append(Problema(
                    n, "alvo", "linha sem alvo_chembl_id nem alvo", "erro"))
                continue
            if not organismo:
                problemas.append(Problema(
                    n, "organismo", f"{nome}: organismo obrigatorio quando o alvo "
                    "e dado por nome (o mesmo nome existe em varias especies)",
                    "erro"))
                continue
        registro = dict(linha)
        registro.update(alvo_chembl_id=chembl_id, alvo=nome, organismo=organismo)
        alvos.append(registro)
    return alvos, problemas


def ler_sintomas(caminho) -> list[str]:
    """Um sintoma por linha; linhas em branco e '#' (comentario) sao ignoradas."""
    sintomas = []
    for linha in Path(caminho).read_text(encoding="utf-8-sig").splitlines():
        linha = linha.strip()
        if linha and not linha.startswith("#"):
            sintomas.append(linha)
    return sintomas
