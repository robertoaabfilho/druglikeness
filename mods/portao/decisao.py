# -*- coding: utf-8 -*-
"""
decisao.py

A politica do portao: de veredito por criterio a `aprovado`.

Semantica, uma so para todos os criterios
-----------------------------------------
Cada modulo de criterio devolve um veredito INTRINSECO, que e o julgamento
dele sobre o item, independente de politica:

    passa          o criterio olhou e nao achou problema
    marca          achou algo que merece ser visto
    reprova        achou algo que o proprio modulo considera impeditivo
    indeterminado  nao havia dado para julgar

A politica diz o que fazer com `marca` e `reprova`:

    reprovar   veredito adverso (marca ou reprova) tira o item
    marcar     nada tira o item; o veredito fica na coluna

`indeterminado` NUNCA tira o item, com politica nenhuma. Ausencia de dado nao
e evidencia de problema -- e o caso de quase todo produto natural na checagem
de bula, e de todo alvo nao humano na expressao.

Por que o padrao e quase todo `marcar`
--------------------------------------
Medido (ver o desenho da etapa 2): reprovar por PAINS apagaria a
amitriptilina; por BRENK, 65% da base; por expressao em tecido nervoso, SERT e
NET. O unico criterio com especificidade medida e o da bula (0 de 6 falsos
positivos no conjunto-controle), e e o unico que reprova por padrao junto com
o SMILES ilegivel, que nem molecula tem.
"""

CRITERIOS_LIGANTE = ("smiles", "alerta_bula", "alerta_estrutural", "toxicidade")
CRITERIOS_ALVO = ("expressao",)
CRITERIOS = CRITERIOS_LIGANTE + CRITERIOS_ALVO

REPROVAR = "reprovar"
MARCAR = "marcar"

POLITICA_PADRAO = {
    "smiles": REPROVAR,          # sem molecula nao ha o que avaliar
    "alerta_bula": REPROVAR,     # 0/6 de falso positivo medido no controle
    "alerta_estrutural": MARCAR,  # reprovaria a amitriptilina e 65% da base
    "toxicidade": MARCAR,        # percentil nao validado neste conjunto
    "expressao": MARCAR,         # reprovaria SERT e NET
}

VEREDITOS_ADVERSOS = ("marca", "reprova")
VEREDITOS = ("passa", *VEREDITOS_ADVERSOS, "indeterminado")


class CriterioDesconhecido(ValueError):
    """Nome de criterio que nao existe; erra cedo em vez de ignorar a flag."""


def montar_politica(promovidos=(), rebaixados=()) -> dict[str, str]:
    """Politica em vigor. `promovidos` passam a reprovar, `rebaixados` a so
    marcar. Nome desconhecido levanta: uma flag com erro de digitacao que
    passasse em silencio faria o usuario acreditar numa politica que nao
    esta valendo."""
    politica = dict(POLITICA_PADRAO)
    for nomes, valor in ((promovidos, REPROVAR), (rebaixados, MARCAR)):
        for nome in nomes:
            nome = nome.strip()
            if not nome:
                continue
            if nome not in politica:
                raise CriterioDesconhecido(
                    f"criterio desconhecido: {nome!r}; os que existem sao "
                    + ", ".join(CRITERIOS))
            politica[nome] = valor
    return politica


def texto_da_politica(politica: dict[str, str]) -> str:
    """A politica numa celula, para que a linha diga sob qual regra foi julgada
    -- duas execucoes com politicas diferentes nao se confundem."""
    return ";".join(f"{c}={politica[c]}" for c in CRITERIOS if c in politica)


def aplicar(vereditos: dict[str, tuple[str, str]], politica: dict[str, str]):
    """(aprovado, motivos) a partir de {criterio: (veredito, valor)}.

    `motivos` lista os criterios que tiraram o item, com o valor que decidiu,
    para que a reprovacao seja lida sem consultar outra coluna."""
    motivos = []
    for criterio, (veredito, valor) in sorted(vereditos.items()):
        if veredito not in VEREDITOS:
            raise ValueError(f"veredito invalido para {criterio}: {veredito!r}")
        if veredito in VEREDITOS_ADVERSOS and politica.get(criterio) == REPROVAR:
            motivos.append(f"{criterio}: {valor}" if valor else criterio)
    return (not motivos), "; ".join(motivos)


def linha_de_portao(identidade: dict, vereditos: dict, politica: dict) -> dict:
    aprovado, motivos = aplicar(vereditos, politica)
    linha = {**identidade, "aprovado": aprovado, "motivos": motivos,
             "politica": texto_da_politica(politica)}
    for criterio, (veredito, valor) in vereditos.items():
        linha[f"{criterio}_veredito"] = veredito
        linha[f"{criterio}_valor"] = valor
    return linha


class Contador:
    """Quantos itens cada criterio reprovaria SE fosse promovido, e quantos
    reprova de fato. E o numero que o resumo imprime: o botao e medido antes
    de ser apertado, em vez de o usuario promover um criterio e descobrir
    depois que perdeu a amitriptilina."""

    def __init__(self, criterios):
        self.criterios = tuple(criterios)
        self.adversos = {c: 0 for c in self.criterios}
        self.vereditos = {c: {v: 0 for v in VEREDITOS} for c in self.criterios}
        self.itens = 0
        self.aprovados = 0

    def contar(self, vereditos: dict, aprovado: bool) -> None:
        self.itens += 1
        self.aprovados += bool(aprovado)
        for criterio, (veredito, _) in vereditos.items():
            if criterio in self.adversos:
                self.vereditos[criterio][veredito] += 1
                if veredito in VEREDITOS_ADVERSOS:
                    self.adversos[criterio] += 1

    def resumo(self, politica: dict) -> dict:
        return {
            "itens": self.itens,
            "aprovados": self.aprovados,
            "reprovados": self.itens - self.aprovados,
            "por_criterio": {
                c: {"politica": politica.get(c, ""),
                    # Com politica `marcar`, este e o numero que o resumo
                    # anuncia como "reprovaria N se promovido"; com
                    # `reprovar`, e o numero que ja reprovou.
                    "adversos": self.adversos[c],
                    **self.vereditos[c]}
                for c in self.criterios},
        }
