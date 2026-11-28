# -*- coding: utf-8 -*-
"""
parametros.py - a conversao entre o formulario do painel e `json/config.json`.

Esta camada existe por causa de um achado da revisao: era o unico pedaco do
painel que ficou dentro do Streamlit em vez de virar biblioteca, e foi ali que
nasceu o unico bug grave de dados do ramo.

**Chave `(str, list)` nao se itera.** `termos_alerta`, `tecidos`, `reprovar` e
`marcar` aceitam string separada por virgula -- e e assim que as CLIs as
documentam (`T[,T]`). Juntar com `", ".join(valor)` sobre uma string a separa
CARACTERE A CARACTERE: `"alerta_bula_forte"` virava
`["a","l","e","r","t","a",...]`, e alerta de bula e um dos dois unicos
criterios que REPROVAM. A politica que o pesquisador acredita estar valendo
deixava de valer.

**Grava-se o que mudou, e nada mais.** O formulario ve as 46 chaves; o arquivo
tem as que o pesquisador escolheu. Gravar o formulario inteiro acrescentava
`sem_rede: false` e removia `efo: []` a cada clique, deixando `_nota_efo`
documentando uma chave que saiu -- a segunda fonte de verdade que o
`config.py` existe para nao ter.
"""

AUSENTE = object()  # distingue "nao esta no arquivo" de "esta e vale None"


def para_texto(valor) -> str:
    """O valor do arquivo como o texto que o campo mostra.

    String sai inteira. Lista sai juntada por `, `. `None` sai vazio."""
    if valor is None:
        return ""
    if isinstance(valor, str):
        return valor
    if isinstance(valor, (list, tuple)):
        return ", ".join(str(v) for v in valor)
    return str(valor)


def de_texto(texto: str, *, anterior=None):
    """O texto do campo de volta ao valor, **no tipo que o arquivo tinha**.

    Quem estava gravado como string continua string: reescrever como lista
    mudaria o arquivo de quem nunca pediu. Sem valor anterior, lista -- e a
    forma mais explicita, e as CLIs aceitam as duas."""
    texto = (texto or "").strip()
    if not texto:
        return None
    if isinstance(anterior, str):
        return texto
    return [parte.strip() for parte in texto.split(",") if parte.strip()]


def somente_mudancas(antes: dict, formulario: dict) -> dict:
    """As chaves do formulario que realmente diferem do arquivo.

    Chave que nao esta no arquivo e que o formulario devolve com o valor
    neutro (`False` para booleano, `None` para o resto) NAO entra: ela chegou
    ao formulario como padrao do widget, nao como escolha de ninguem.

    `0` e `0.0` sao valores, nao ausencias: `corte: 0` e o modo denso, que
    inverte o resultado do enriquecimento."""
    mudancas = {}
    for chave, novo in (formulario or {}).items():
        if chave.startswith("_"):
            continue  # nota: documentacao, nao parametro
        velho = antes.get(chave, AUSENTE) if antes else AUSENTE
        if velho is AUSENTE:
            if novo is None or novo is False:
                continue  # padrao do widget, nao escolha
            mudancas[chave] = novo
            continue
        if novo is None and _vazio(velho):
            continue  # campo vazio sobre valor vazio: nada mudou
        if not _igual(velho, novo):
            mudancas[chave] = novo
    return mudancas


def _vazio(valor) -> bool:
    return valor is None or (isinstance(valor, (list, tuple, str))
                             and len(valor) == 0)


def _igual(a, b) -> bool:
    if isinstance(a, bool) or isinstance(b, bool):
        return a is b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return float(a) == float(b)
    return a == b
