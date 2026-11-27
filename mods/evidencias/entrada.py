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
* Alvos sao deduplicados pela identidade usada a jusante (ID do ChEMBL, ou
  nome+organismo quando nao ha ID): a tabela tem UMA linha por evidencia, e um
  alvo repetido na entrada dobraria cada ensaio medido dele.
* `alvo_chembl_id` e normalizado para maiusculas AQUI, no unico ponto de
  leitura: `chembl222` e `CHEMBL222` sao o mesmo alvo, e a jusante ninguem
  precisa repetir a normalizacao.
* Todo alvo sai com `origem_alvo` (`usuario` por padrao; a coluna `origem` do
  `alvos_enriquecidos.tsv` e aceita com esse nome, que e o da saida). Sem isso
  a procedencia do alvo se perderia ao realimentar a lista enriquecida.
* Linha com erro e descartada e relatada; nada falha em silencio.
* Se o alvo existe ou nao no espaco do modelo PP4CKD NAO e decidido aqui
  (vem do arquivo de rotulos do modelo, onde for preciso).
"""

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from mods.quimica import forma_parental, inchikey_de, molecula, texto


@dataclass
class Problema:
    linha: int
    campo: str
    mensagem: str
    gravidade: Literal["erro", "aviso"]


def _ler_tabela_lazy(caminho):
    """Le TSV (ou CSV, com aviso). Devolve (cabecalho, gerador, problemas): o
    gerador produz (linha, dict) sob demanda, entao quem para cedo nao paga o
    resto do arquivo. `problemas` cresce conforme o gerador avanca. Linhas com
    numero de celulas diferente do cabecalho sao descartadas."""
    caminho = Path(caminho)
    problemas: list[Problema] = []
    csv_ = caminho.suffix.lower() == ".csv"
    if csv_:
        problemas.append(Problema(
            1, "arquivo", "arquivo .csv aceito, mas o formato e TSV: nomes com "
            "virgula quebram a leitura; prefira .tsv", "aviso"))
    def gerar():
        # O arquivo pertence ao gerador: close() do gerador fecha o arquivo.
        # O primeiro item produzido e o cabecalho (lido logo abaixo).
        with open(caminho, newline="", encoding="utf-8-sig") as f:
            # TSV nao tem citacao: aspas fazem parte do nome. Com o dialeto
            # padrao, `"cis" isomer of X` seria lido como `cis isomer of X` --
            # o rotulo de identidade reescrito em silencio. No .csv a citacao
            # continua valendo: e dela que a virgula do nome depende.
            leitor = csv.reader(f, delimiter="," if csv_ else "\t",
                                **({} if csv_ else {"quoting": csv.QUOTE_NONE}))
            cabecalho: list[str] = []
            for celulas in leitor:
                if any(c.strip() for c in celulas):
                    cabecalho = [c.strip() for c in celulas]
                    # Convencao do repositorio: cabecalho '#nome\t...' (um '#' so).
                    if cabecalho[0].startswith("#"):
                        cabecalho[0] = cabecalho[0][1:].strip()
                    break
            yield cabecalho
            if not cabecalho:
                return
            for celulas in leitor:
                if not any(c.strip() for c in celulas):
                    continue
                if len(celulas) != len(cabecalho):
                    dica = " (virgula nao citada em .csv?)" if csv_ else ""
                    problemas.append(Problema(
                        leitor.line_num, "linha", f"{len(celulas)} celulas, "
                        f"cabecalho tem {len(cabecalho)}{dica}; linha descartada",
                        "erro"))
                    continue
                yield leitor.line_num, dict(zip(cabecalho, celulas))

    gerador = gerar()
    cabecalho = next(gerador)
    if not cabecalho:
        problemas.append(Problema(1, "arquivo", "arquivo vazio", "erro"))
    return cabecalho, gerador, problemas


def _ler_tabela(caminho) -> tuple[list[str], list[tuple[int, dict]], list[Problema]]:
    """Como `_ler_tabela_lazy`, mas le tudo e devolve lista."""
    cabecalho, gerador, problemas = _ler_tabela_lazy(caminho)
    return cabecalho, list(gerador), problemas


def _exigir(cabecalho: list[str], colunas: list[str], problemas: list[Problema]) -> bool:
    faltam = [c for c in colunas if c not in cabecalho]
    for c in faltam:
        problemas.append(Problema(1, c, f"coluna obrigatoria '{c}' ausente", "erro"))
    return not faltam


def ler_ligantes(caminho, limite: int | None = None, *,
                 manter_ilegiveis: bool = False
                 ) -> tuple[list[dict], list[Problema]]:
    """`limite`: para a leitura apos N ligantes aceitos; linhas alem disso nao
    sao lidas nem geram problemas.

    `manter_ilegiveis`: mantem na lista o ligante cujo SMILES o RDKit nao le,
    com `inchikey_parental` vazio, e rebaixa o problema dele a `aviso`.

    Padrao FALSO, que e o que a etapa 1 precisa: sem molecula nao ha
    fingerprint nem InChIKey parental, logo nao ha evidencia prevista nem como
    casar evidencia medida -- o ligante seria uma linha que nada pode
    preencher, e descarta-lo com erro de linha e o certo.

    A etapa 2 passa VERDADEIRO, porque ali a pergunta e outra: o portao julga
    o ligante, e um ligante que desaparece da tabela de vereditos nao foi
    reprovado, foi apagado -- que e anti-objetivo declarado do desenho da
    etapa 2. Ele sai com as colunas de propriedade vazias e o criterio
    `smiles` reprovando, com o motivo escrito.

    Vale so para SMILES PRESENTE e ilegivel. Linha sem a coluna `smiles`
    preenchida continua sendo erro em qualquer etapa: sem nome de molecula nem
    estrutura nao ha o que julgar nem o que relatar."""
    cabecalho, linhas, problemas = _ler_tabela_lazy(caminho)
    if cabecalho and not _exigir(cabecalho, ["nome", "smiles"], problemas):
        linhas.close()
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
            # Mantido, o problema e `aviso`: o ligante esta na saida e sera
            # julgado, entao a execucao nao e parcial por causa dele. Descartado,
            # e `erro`, e quem conta erros torna a execucao parcial -- a
            # diferenca entre "nao avaliei este" e "perdi este".
            problemas.append(Problema(
                n, "smiles", f"{nome}: smiles que o RDKit nao le: {smiles}"
                + (" (mantido para julgamento, sem propriedades)"
                   if manter_ilegiveis else ""),
                "aviso" if manter_ilegiveis else "erro"))
            if not manter_ilegiveis:
                continue
            # Sem chave parental nao ha como deduplicar: cada linha ilegivel
            # fica, porque duas delas podem ser moleculas diferentes.
            registro = dict(linha)
            registro.update(nome=nome, smiles=smiles, inchikey="",
                            inchikey_parental="")
            ligantes.append(registro)
            if limite is not None and len(ligantes) >= limite:
                break
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
        if limite is not None and len(ligantes) >= limite:
            break  # nao puxa a proxima linha: ela nao pode gerar problema
    linhas.close()
    for chave, nomes in descartados.items():
        problemas.append(Problema(
            0, "nome", f"ligante duplicado pela forma parental ({chave}): mantido "
            f"'{vistos[chave]}', descartado(s): {', '.join(nomes)}", "aviso"))
    return ligantes, problemas


_CABECALHOS_DE_ALVO = {"alvo_chembl_id", "alvo", "organismo"}
_RE_CHEMBL = re.compile(r"CHEMBL\d+", re.IGNORECASE)

# Procedencias que caracterizam um alvo VINDO do enriquecimento (as que
# `escrever_alvos_enriquecidos` grava em `origem`); `usuario` e o padrao.
ORIGENS_DE_ENRIQUECIMENTO = ("open_targets", "reactome")
ORIGEM_PADRAO = "usuario"


def _origem_do_alvo(linha: dict) -> str:
    """Procedencia do alvo sob um nome so, `origem_alvo`. Aceita tambem a
    coluna `origem`, que e como a etapa grava o `alvos_enriquecidos.tsv`: sem
    isso realimentar aquele arquivo perderia a procedencia de cada alvo."""
    bruta = texto(linha.get("origem_alvo")) or texto(linha.get("origem"))
    return bruta.strip().lower() or ORIGEM_PADRAO


def _identidade(alvo: dict):
    """Como o alvo e reconhecido a jusante: pelo ID, ou pelo par
    nome+organismo quando nao ha ID."""
    if alvo["alvo_chembl_id"]:
        return ("id", alvo["alvo_chembl_id"])
    return ("nome", alvo["alvo"].casefold(), alvo["organismo"].casefold())


def _deduplicar_alvos(alvos: list[dict], problemas: list[Problema]) -> list[dict]:
    """Mantem a primeira ocorrencia de cada alvo e relata o que foi dobrado.
    `evidencias.tsv` tem uma linha por evidencia: um alvo repetido dobraria
    cada ensaio medido dele, e a duplicata seria indistinguivel de dois
    ensaios de verdade."""
    unicos: list[dict] = []
    vistos: dict = {}
    dobrados: dict = {}
    for alvo in alvos:
        chave = _identidade(alvo)
        if chave in vistos:
            dobrados.setdefault(chave, []).append(
                alvo["alvo_chembl_id"] or f"{alvo['alvo']} ({alvo['organismo']})")
            continue
        vistos[chave] = alvo
        unicos.append(alvo)
    for chave, repetidos in dobrados.items():
        rotulo = (vistos[chave]["alvo_chembl_id"]
                  or f"{vistos[chave]['alvo']} ({vistos[chave]['organismo']})")
        problemas.append(Problema(
            0, "alvo", f"alvo repetido na entrada ({rotulo}): mantida a primeira "
            f"ocorrencia, dobrada(s) {len(repetidos)}; sem isso cada evidencia "
            "medida dele sairia repetida", "aviso"))
    return unicos


def _e_lista_de_alvos(caminho) -> bool:
    """Formato 'um alvo por linha' (o do preditor de alvos), em vez de tabela.
    Tabela: arquivo .csv, ou primeira linha de conteudo com tabulacao, ou cuja
    primeira linha (dividida em tabulacao E virgula) tem um nome de coluna."""
    if Path(caminho).suffix.lower() == ".csv":
        return False
    with open(caminho, encoding="utf-8-sig") as f:
        for linha in f:
            if not linha.strip():
                continue
            if "\t" in linha:
                return False
            celulas = re.split(r"[,\t]", linha.strip().lstrip("#"))
            return not any(c.strip().lower() in _CABECALHOS_DE_ALVO for c in celulas)
    return False


def _ler_lista_de_alvos(caminho) -> tuple[list[dict], list[Problema]]:
    """'ID  # comentario' ou 'nome | organismo[, organismo]'; '#' comenta."""
    alvos: list[dict] = []
    problemas: list[Problema] = []
    with open(caminho, encoding="utf-8-sig") as f:
        for n, linha in enumerate(f, start=1):
            linha = linha.split("#", 1)[0].strip()
            if not linha:
                continue
            termo, _, organismos = linha.partition("|")
            termo = termo.strip()
            if _RE_CHEMBL.fullmatch(termo):
                alvos.append({"alvo_chembl_id": termo.upper(), "alvo": "",
                              "organismo": "", "origem_alvo": ORIGEM_PADRAO})
                continue
            lista = [o.strip() for o in organismos.split(",") if o.strip()]
            if not lista:
                problemas.append(Problema(
                    n, "organismo", f"{termo}: organismo obrigatorio quando o alvo "
                    "e dado por nome ('nome | organismo'); o mesmo nome existe em "
                    "varias especies", "erro"))
                continue
            alvos.extend({"alvo_chembl_id": "", "alvo": termo, "organismo": o,
                          "origem_alvo": ORIGEM_PADRAO} for o in lista)
    if not alvos and not problemas:
        problemas.append(Problema(1, "arquivo", "arquivo vazio", "erro"))
    return _deduplicar_alvos(alvos, problemas), problemas


def ler_alvos(caminho) -> tuple[list[dict], list[Problema]]:
    if _e_lista_de_alvos(caminho):
        return _ler_lista_de_alvos(caminho)
    cabecalho, linhas, problemas = _ler_tabela(caminho)
    if cabecalho and "alvo_chembl_id" not in cabecalho and not _exigir(
            cabecalho, ["alvo", "organismo"], problemas):
        return [], problemas
    alvos: list[dict] = []
    for n, linha in linhas:
        chembl_id = texto(linha.get("alvo_chembl_id")).upper()
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
        registro.update(alvo_chembl_id=chembl_id, alvo=nome, organismo=organismo,
                        origem_alvo=_origem_do_alvo(linha))
        alvos.append(registro)
    return _deduplicar_alvos(alvos, problemas), problemas


def ler_sintomas(caminho) -> list[str]:
    """Um sintoma por linha; linhas em branco e '#' (comentario) sao ignoradas."""
    sintomas = []
    for linha in Path(caminho).read_text(encoding="utf-8-sig").splitlines():
        linha = linha.strip()
        if linha and not linha.startswith("#"):
            sintomas.append(linha)
    return sintomas
