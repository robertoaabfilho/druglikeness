# -*- coding: utf-8 -*-
"""
faixa.py

Em que faixa de evidencia um composto cai, e com que valor ele ordena
DENTRO dela.

Desenho: `docs/superpowers/specs/2026-11-26-fase3-ranqueamento-design.md`.

As cinco faixas
---------------
    1a  evidencia MEDIDA em alvo de `origem_alvo = usuario`   -> max pchembl
    1b  evidencia MEDIDA em alvo do enriquecimento            -> max pchembl
    2a  so PREVISTA, em alvo do usuario       -> max p x (1 - p_negativo)
    2b  so PREVISTA, em alvo do enriquecimento -> max p x (1 - p_negativo)
    3   nenhuma evidencia em alvo afirmativo  -> NAO ranqueia (valor None)

Tres decisoes que vivem aqui, com o motivo de cada uma
------------------------------------------------------
* **O canal negativo nao multiplica a faixa 1.** `pchembl` e medido e
  `p_negativo` e previsto; multiplica-los poria medido e previsto no mesmo
  numero, o que a etapa 1 proibe. E o ganho medido (AUC 0,806) foi obtido
  inteiramente sobre valores previstos: e um numero da faixa 2, e aplica-lo a
  `pchembl` estenderia a outra grandeza uma regra medida noutra.
* **A procedencia e sub-faixa, nao peso.** O requisito 2 diz que pedir e
  declarar, logo o alvo do usuario vale mais que o do enriquecimento. Dizer
  isso com um peso exigiria uma constante que n = 12 nao calibra; dizer com
  sub-faixa nao exige nenhuma.
* **`max` e nao soma** para agregar varios alvos: e o que foi medido, e soma
  premiaria o composto promiscuo -- que a medicao do canal humano mostrou ser
  enganoso (AUC 0,389 sozinho e invertido).

Linha medida sem `pchembl`
--------------------------
Caso presente em 1 das 49 linhas medidas do `evidencias_controle.tsv`. Ela nao
e comparavel, logo nao ordena -- mas rebaixar o composto para a faixa 3 diria
"sem evidencia" sobre um composto que FOI medido. Entao ela nao entra no valor,
o composto cai para a faixa prevista se nao houver outra medida comparavel, e
`medida_sem_pchembl` registra o caso.
"""

from mods.evidencias.entrada import ORIGENS_DE_ENRIQUECIMENTO

FAIXAS = ("1a", "1b", "2a", "2b", "3")
ORIGENS_DO_USUARIO = ("usuario",)

MEDIDO = "medido"
PREVISTO = "previsto"


def _numero(texto):
    """float de uma celula de TSV, ou None quando ela nao converte.

    Celula vazia e celula com texto dao o mesmo resultado de proposito: as
    duas significam "nao da para comparar", e e isso que decide a faixa."""
    try:
        return float((texto or "").strip())
    except (TypeError, ValueError):
        return None


def _melhor(linhas, campo):
    """(valor, linha) do maior `campo` comparavel, ou (None, None).

    Devolve tambem a linha para que o alvo que decidiu o valor possa ser
    nomeado na saida -- um escore sem o alvo responsavel nao e auditavel."""
    melhor_valor = None
    melhor_linha = None
    for linha in linhas:
        valor = _numero(linha.get(campo))
        if valor is None:
            continue
        if melhor_valor is None or valor > melhor_valor:
            melhor_valor, melhor_linha = valor, linha
    return melhor_valor, melhor_linha


def _nome_do_alvo(linha) -> str:
    if linha is None:
        return ""
    id_alvo = (linha.get("alvo_chembl_id") or "").strip()
    nome = (linha.get("alvo") or "").strip()
    return f"{id_alvo} {nome}".strip() if id_alvo else nome


def classificar(linhas: list[dict], *, p_negativo: float | None) -> dict:
    """Faixa e valor de UM composto, a partir das linhas de evidencia dele.

    `p_negativo` e o agregado do canal negativo, ou None quando a execucao nao
    tem canal (sem pesos do PP4CKD, SMILES ilegivel, ou `--sem-canal-negativo`).
    None NAO e zero: zero afirmaria "nenhuma atividade indesejada", e ausencia
    de dado nao e afirmacao. Sem canal, a faixa prevista ordena pelo `p` so.
    """
    por_faixa = {}
    for faixa, tipo, origens in (("1a", MEDIDO, ORIGENS_DO_USUARIO),
                                 ("1b", MEDIDO, ORIGENS_DE_ENRIQUECIMENTO),
                                 ("2a", PREVISTO, ORIGENS_DO_USUARIO),
                                 ("2b", PREVISTO, ORIGENS_DE_ENRIQUECIMENTO)):
        por_faixa[faixa] = [
            l for l in linhas
            if (l.get("tipo_evidencia") or "").strip() == tipo
            and (l.get("origem_alvo") or "").strip() in origens]

    medidas = por_faixa["1a"] + por_faixa["1b"]
    sem_pchembl = any(_numero(l.get("pchembl")) is None for l in medidas)

    resultado = {"faixa": "3", "valor": None, "alvo": "", "p_afirmativo": "",
                 "pchembl": "", "medida_sem_pchembl": "sim" if sem_pchembl else "",
                 "observacao": ""}

    for faixa in ("1a", "1b"):
        pchembl, linha = _melhor(por_faixa[faixa], "pchembl")
        if pchembl is not None:
            resultado.update(faixa=faixa, valor=pchembl, pchembl=pchembl,
                             alvo=_nome_do_alvo(linha))
            return resultado

    for faixa in ("2a", "2b"):
        p, linha = _melhor(por_faixa[faixa], "valor")
        if p is None:
            continue
        valor = p if p_negativo is None else p * (1.0 - p_negativo)
        resultado.update(faixa=faixa, valor=valor, p_afirmativo=p,
                         alvo=_nome_do_alvo(linha))
        if p_negativo is None:
            resultado["observacao"] = "sem canal negativo"
        return resultado

    if sem_pchembl:
        resultado["observacao"] = "so medida nao comparavel"
    return resultado


def chave(resultado: dict, inchikey: str) -> tuple:
    """Chave ASCENDENTE de ordenacao: `sorted(..., key=...)` sem `reverse`.

    `(indice da faixa, -valor, inchikey)`. O `inchikey` entra para que o
    empate seja estavel: sem ele, duas execucoes sobre a mesma entrada podem
    dar arquivos diferentes, e a repeticao byte a byte deixa de valer.

    Ascendente, e nao descendente com o `inchikey` invertido, porque inverter
    a string caractere a caractere para manter o desempate crescente seria
    codigo sem leitor.
    """
    valor = resultado.get("valor")
    return (FAIXAS.index(resultado["faixa"]),
            float("inf") if valor is None else -float(valor),
            inchikey or "")
