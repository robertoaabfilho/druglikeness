# -*- coding: utf-8 -*-
"""
medido.py

Evidencia MEDIDA: ensaios do ChEMBL lidos de `evidencias_chembl.csv` (uma linha
por ensaio, nao por par ligante/alvo).

Regras
------
* O alvo casa pelo PAR (nome, organismo): o arquivo nao traz ID de alvo e o
  mesmo nome existe em varias especies. Nunca so pelo nome.
* UM alvo por par (nome, organismo). Dois pedidos para o mesmo alvo (um por
  ID, outro por nome) sao dobrados e relatados: a tabela tem uma linha por
  evidencia, e a linha repetida pareceria um segundo ensaio.
* O ligante casa pela InChIKey completa (estereo incluso).
* Todo registro e mantido e rotulado; nada e descartado aqui. A fase 2 decide.
  Rotulos: `relacao` (um `>` e limite, nao medida), `triagem_hts`, `unidade`,
  `tipo_ensaio` e, em `procedencia`, `estereoisomero` quando a molecula
  ensaiada nao e a do ligante.
* `correspondencia_ligante` traz a `correspondencia` da fonte: "exata" ou
  "estereoisomero". Um estereoisomero NAO e uma correspondencia de identidade:
  e outro composto, e a atividade medida nele nao e a do ligante. As linhas
  ficam (a fase 2 decide) mas filtraveis por esta coluna, sem ler prosa.
  Linhas previstas deixam a coluna vazia.
* `valor` e `unidade` NAO sao comparaveis entre linhas (120 nM e 78% nao se
  ordenam). `pchembl` e a unica coluna comparavel e fica vazia quando a fonte
  nao a traz; nunca e calculada nem preenchida com outra coisa.
* `confianca_ensaio` fica vazia: o arquivo de origem nao a traz.
* Erro de resolucao de UM ligante ou alvo (`ListaTruncada`, `RespostaInesperada`,
  rede) vira problema em `contagens["problemas"]` e o resto segue.
* Leitura em fluxo: o arquivo tem ~60 MB e a saida pode ter milhoes de linhas.
"""

import csv
import math
import unicodedata
from pathlib import Path
from typing import Callable, Iterator

from mods.evidencias.chembl import ListaTruncada, RespostaInesperada
from mods.evidencias.esquema import formatar_atividade

# Falhas de uma consulta isolada: registradas e puladas, nunca abortam a carga.
ERROS_DE_ITEM = (ListaTruncada, RespostaInesperada, LookupError, OSError)

csv.field_size_limit(1 << 24)


def e_triagem(descricao: str) -> bool:
    """Varredura HTS de painel, em vez de ensaio dirigido a um alvo."""
    return "PUBCHEM_BIOASSAY" in descricao or "qHTS" in descricao


COLUNAS_OBRIGATORIAS = ("InChIKey", "alvo", "organismo_alvo", "medida", "valor",
                        "pchembl", "correspondencia", "descricao_ensaio")


def validar_cabecalho(caminho_csv: str) -> None:
    """ValueError nomeando as colunas que faltam: com o arquivo errado toda
    linha cairia como 'sem alvo' e o resultado pareceria 'sem evidencia'."""
    with open(Path(caminho_csv).expanduser(), encoding="utf-8-sig", newline="") as f:
        cabecalho = next(csv.reader(f), [])
    faltam = [c for c in COLUNAS_OBRIGATORIAS if c not in cabecalho]
    if faltam:
        raise ValueError(f"{caminho_csv}: nao parece evidencias_chembl.csv; "
                         f"colunas ausentes: {', '.join(faltam)}")


def _inteiro_se_inteiro(n: float):
    return int(n) if n.is_integer() else n


def _texto_numerico(texto: str) -> str:
    """'2021.0' -> '2021'; texto que nao e numero fica como veio."""
    n = _numero(texto)
    return str(_inteiro_se_inteiro(n)) if n is not None else (texto or "").strip()


def _correspondencia(texto: str) -> str:
    """'exata', 'estereoisomero' ou 'desconhecida' (nunca chuta estereoisomero).
    Aceita a forma acentuada em NFC ou NFD."""
    base = unicodedata.normalize("NFKD", texto or "")
    base = "".join(ch for ch in base if not unicodedata.combining(ch))
    base = base.strip().casefold()
    return base if base in ("exata", "estereoisomero") else "desconhecida"


def _numero(texto: str):
    """float finito, ou None. Aceita so ponto decimal, como a fonte."""
    try:
        n = float((texto or "").strip())
    except ValueError:
        return None
    return n if math.isfinite(n) else None


def _chave_alvo(nome: str, organismo: str) -> tuple[str, str]:
    return (nome or "").strip().casefold(), (organismo or "").strip().casefold()


def _indexar_alvos(alvos, rotulos, resolver_alvo, problemas, avisos=None):
    """{(nome, organismo): alvo resolvido}; alvo sem nome e resolvido pelo ID.

    UM alvo por par (nome, organismo). Dois pedidos que resolvem para o mesmo
    alvo (um por ID, outro por nome) dobrariam cada ensaio dele, e a linha
    repetida seria indistinguivel de dois ensaios de verdade numa tabela de uma
    linha por evidencia. O ID chega normalizado da leitura (`ler_alvos`)."""
    por_id = {r["chembl_id"]: r for r in rotulos or [] if r.get("chembl_id")}
    indice: dict[tuple[str, str], dict] = {}
    for alvo in alvos:
        chembl_id = (alvo.get("alvo_chembl_id") or "").strip()
        nome, organismo, tipo = alvo.get("alvo", ""), alvo.get("organismo", ""), ""
        if chembl_id:
            r = por_id.get(chembl_id)
            if r is not None:
                nome, organismo, tipo = r["nome"], r["organismo"], r.get("tipo", "")
            elif not (nome and organismo) and resolver_alvo is not None:
                try:
                    achado = resolver_alvo(chembl_id)
                except ERROS_DE_ITEM as e:
                    problemas.append(f"alvo {chembl_id}: {type(e).__name__}: {e}")
                    continue
                nome, organismo = achado["nome"], achado["organismo"]
                tipo = achado.get("tipo", "")
        if not (nome and organismo):
            problemas.append(f"alvo {chembl_id or nome or '?'}: sem nome e organismo "
                             "para casar com o arquivo; nenhuma linha medida dele")
            continue
        chave = _chave_alvo(nome, organismo)
        if chave in indice:
            ja = indice[chave]
            # Deduplicacao deliberada, nao falha: vai para `avisos`, que nao
            # torna a execucao parcial.
            (problemas if avisos is None else avisos).append(
                f"alvo {chembl_id or nome}: e o mesmo alvo que "
                f"{ja['alvo_chembl_id'] or ja['alvo']} ({nome}, {organismo}); "
                "usado uma vez so, para nao repetir cada ensaio medido dele")
            continue
        indice[chave] = {
            "alvo_chembl_id": chembl_id, "alvo": nome, "organismo": organismo,
            "tipo_alvo": tipo, "origem_alvo": alvo.get("origem_alvo", "")}
    return indice


def linhas_medidas(caminho_csv: str, ligantes: list[dict], alvos: list[dict], *,
                   rotulos: list[dict] | None = None,
                   resolver_alvo: Callable[[str], dict] | None = None,
                   resolver_ligante: Callable[[str, str], str | None] | None = None,
                   contagens: dict | None = None) -> Iterator[dict]:
    """Uma linha por ensaio medido de um ligante sobre um alvo pedidos.

    `rotulos`: rotulos do modelo, para dar nome/organismo a alvo dado so por ID.
    `resolver_alvo(chembl_id) -> {nome, organismo, tipo}` e
    `resolver_ligante(inchikey, parental) -> chembl_id | None`: opcionais (rede);
    suas falhas por item entram em `contagens["problemas"]`.
    `contagens` recebe: linhas_lidas, linhas_emitidas, sem_alvo_utilizavel,
    ligante_ou_alvo_nao_pedido, problemas (lista de textos), avisos (avisos
    que nao tornam a execucao parcial, como alvo pedido duas vezes)."""
    c = contagens if contagens is not None else {}
    for k in ("linhas_lidas", "linhas_emitidas", "sem_alvo_utilizavel",
              "ligante_ou_alvo_nao_pedido"):
        c.setdefault(k, 0)
    problemas = c.setdefault("problemas", [])
    avisos = c.setdefault("avisos", [])

    indice_alvos = _indexar_alvos(alvos, rotulos, resolver_alvo, problemas, avisos)
    por_chave: dict[str, dict] = {}
    for campo in ("inchikey", "inchikey_parental"):  # a chave exata tem prioridade
        for lig in ligantes:
            chave = (lig.get(campo) or "").strip().upper()
            if chave:
                por_chave.setdefault(chave, lig)
    ids_ligante: dict[str, str] = {}

    def id_do_ligante(lig: dict) -> str:
        if lig.get("chembl_id"):
            return lig["chembl_id"]
        ik, parental = lig.get("inchikey", ""), lig.get("inchikey_parental", "")
        if resolver_ligante is None:
            return ""
        if (ik, parental) not in ids_ligante:
            try:
                ids_ligante[(ik, parental)] = resolver_ligante(ik, parental) or ""
            except ERROS_DE_ITEM as e:
                problemas.append(f"ligante {lig.get('nome', ik)}: "
                                 f"{type(e).__name__}: {e}")
                ids_ligante[(ik, parental)] = ""
        return ids_ligante[(ik, parental)]

    validar_cabecalho(caminho_csv)
    if not indice_alvos or not por_chave:
        return
    desconhecidas: set[str] = set()
    with open(Path(caminho_csv).expanduser(), encoding="utf-8-sig", newline="") as f:
        for linha in csv.DictReader(f):
            c["linhas_lidas"] += 1
            nome = (linha.get("alvo") or "").strip()
            if not nome or nome == "Unchecked":
                c["sem_alvo_utilizavel"] += 1
                continue
            lig = por_chave.get((linha.get("InChIKey") or "").strip().upper())
            alvo = indice_alvos.get(_chave_alvo(nome, linha.get("organismo_alvo", "")))
            if lig is None or alvo is None:
                c["ligante_ou_alvo_nao_pedido"] += 1
                continue

            medida = (linha.get("medida") or "").strip()
            relacao = (linha.get("relacao") or "").strip()
            unidade = (linha.get("unidade") or "").strip()
            valor = _numero(linha.get("valor"))
            pchembl = _numero(linha.get("pchembl"))
            bruta = (linha.get("correspondencia") or "").strip()
            corresp = _correspondencia(bruta)
            exata = corresp == "exata"
            if corresp == "desconhecida" and bruta not in desconhecidas:
                desconhecidas.add(bruta)
                problemas.append(f"correspondencia nao reconhecida {bruta!r}: linhas "
                                 "marcadas 'desconhecida', sem ID de ligante")
            ensaiada = (linha.get("molecula_ensaiada") or "").strip()
            descricao = linha.get("descricao_ensaio") or ""
            procedencia = "medido (ChEMBL, evidencias_chembl.csv"
            if exata:
                procedencia += ")"
                id_lig = lig.get("chembl_id") or ensaiada or id_do_ligante(lig)
            elif corresp == "desconhecida":
                procedencia += "; correspondencia desconhecida)"
                id_lig = lig.get("chembl_id") or ""
            else:
                # Outra molecula (estereoisomero) foi ensaiada: o ID dela nao e o do ligante.
                procedencia += f"; estereoisomero da molecula ensaiada {ensaiada})"
                id_lig = id_do_ligante(lig)

            c["linhas_emitidas"] += 1
            yield {
                "ligante": lig.get("nome", ""),
                "ligante_inchikey_parental": lig.get("inchikey_parental", ""),
                "ligante_chembl_id": id_lig,
                "correspondencia_ligante": corresp,
                "alvo_chembl_id": alvo["alvo_chembl_id"],
                "alvo": alvo["alvo"],
                "organismo": alvo["organismo"],
                "tipo_alvo": alvo["tipo_alvo"],
                "origem_alvo": alvo["origem_alvo"],
                "tipo_evidencia": "medido",
                "fonte": "ChEMBL",
                "atividade_relatada": formatar_atividade(
                    medida, relacao, valor, unidade),
                "medida": medida,
                "relacao": relacao,
                "valor": "" if valor is None else _inteiro_se_inteiro(valor),
                "unidade": unidade,
                "pchembl": "" if pchembl is None else _inteiro_se_inteiro(pchembl),
                "tipo_ensaio": (linha.get("tipo_ensaio") or "").strip(),
                "ensaio_descricao": descricao,
                "triagem_hts": "sim" if e_triagem(descricao) else "nao",
                "confianca_ensaio": "",
                "n_moleculas_treino": "",
                "documento": (linha.get("documento_chembl") or "").strip(),
                "pmid": _texto_numerico(linha.get("PMID")),
                "ano": _texto_numerico(linha.get("ano")),
                "procedencia": procedencia,
            }
