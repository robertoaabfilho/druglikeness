# -*- coding: utf-8 -*-
"""
negativo.py

O canal negativo DEDUZIDO: atividade em alvo que ninguem pediu.

    negativo = alvo NAO afirmativo
             & organismo preenchido E != "homo sapiens"
             & n_moleculas_treino >= MINIMO_TREINO
             & pref_name != pref_name de qualquer alvo afirmativo

Medido no controle de 12 (adendo de 2026-11-26 aos requisitos da fase 3,
ECFP4, `ecfp4_fold_4.pt`, varrendo os 7.676 alvos do modelo):

| canal negativo | AUC (trata > causa) |
|---|---|
| nenhum -- so os 9 alvos pedidos | 0,750 |
| **nao pedido + NAO HUMANO** (1.169 brutos; 1.129 no canal final) | **0,806** |
| nao pedido + humano (2.376 alvos) | 0,694 -- PIORA |
| contagem de promiscuidade | 0,583 -- PIORA, e zera a amitriptilina |
| painel curado de seguranca | 0,750 -- nada |

0,806 e exatamente o que o alvo do patogeno DECLARADO dava: deduzir empata com
declarar. O lado humano fica fora porque, sozinho e invertido, ele da AUC
0,389 -- atividade humana nao pedida esta associada a TRATAR, e o motivo tem
nome: amitriptilina, 0,910. Punir isso repetiria o erro do PAINS na etapa 2,
apagando o positivo de referencia.

Nao e a regra "nao humano = negativo" que os requisitos derrubaram: o filtro de
PROCEDENCIA vem primeiro. No caso do *C. jejuni*, as enzimas que o usuario
listou sao pedidas, logo saem do canal antes do filtro de organismo.

A exclusao de ortologo custa 9 alvos dos 1.169 e nao muda a AUC: e rede de
seguranca de graca que cobre o CHEMBL304 (NET de rato), que e nao humano e
AFIRMATIVO quando o NET humano foi pedido.

Ressalva de validade: n = 12 e os seis causadores do controle sao todos
antibioticos, logo a atividade nao humana os separa por construcao do
conjunto. Daqui nao sai regra geral, e o canal e desligavel.
"""

import numpy as np

from mods.evidencias.resumo import MINIMO_TREINO

HUMANO = "homo sapiens"      # comparado em casefold: ver `mascara`


def mascara(rotulos: list[dict], idx_afirmativos: list[int], *,
            minimo_treino: int = MINIMO_TREINO) -> np.ndarray:
    """Quais alvos do modelo entram no canal negativo, como vetor de bool.

    `rotulos` e o que `prever_alvos.carregar_rotulos` devolve, indexado por
    `target_idx` -- a mesma indexacao da saida da rede, para que a mascara
    possa ser aplicada direto ao vetor de probabilidades.
    """
    afirmativos = set(idx_afirmativos)
    nomes_afirmativos = {(rotulos[i].get("nome") or "").strip().casefold()
                         for i in afirmativos}
    nomes_afirmativos.discard("")
    return np.array(
        [(i not in afirmativos)
         # Organismo vazio NAO entra: ausencia de dado nao e evidencia de "nao
         # humano". O caso real e `CHEMBL612545 "Unchecked"`, o balde do ChEMBL
         # para ensaios nao curados -- organismo vazio e 217.067 moleculas
         # ativas, o alvo mais treinado do modelo, que liderava o canal de
         # quase todo composto. Medido: excluir os 31 alvos sem organismo
         # mantem a AUC em 0,806 e o maximo passa a ser patogeno de verdade
         # (S. aureus, E. coli, M. tuberculosis).
         and bool((r.get("organismo") or "").strip())
         # casefold, como `predito.py` e `medido.py` ja fazem: um
         # `target_labels.tsv` com "homo sapiens" jogaria alvos HUMANOS para
         # dentro do canal, que e o desenho medido como prejudicial.
         and (r.get("organismo") or "").strip().casefold() != HUMANO
         and int(r.get("n_ativos") or 0) >= minimo_treino
         and (r.get("nome") or "").strip().casefold() not in nomes_afirmativos
         for i, r in enumerate(rotulos)], dtype=bool)


def agregar(probabilidades, mascara_do_canal) -> float | None:
    """`max` das probabilidades dentro da mascara, ou None se ela for vazia.

    None, e nunca 0.0: zero afirmaria "nenhuma atividade indesejada", e
    mascara vazia significa que nao havia o que olhar. Ausencia de dado nao e
    evidencia -- a mesma regra do `indeterminado` do portao da etapa 2."""
    p = np.asarray(probabilidades)
    m = np.asarray(mascara_do_canal, dtype=bool)
    if not m.any():
        return None
    return float(p[m].max())


def alvo_do_maximo(rotulos: list[dict], probabilidades,
                   mascara_do_canal) -> str:
    """Nome do alvo que produziu o maximo do canal, ou "" se a mascara e vazia.

    Existe para que a penalidade seja auditavel: um composto punido sem dizer
    por qual alvo nao e um resultado que alguem possa conferir."""
    p = np.asarray(probabilidades)
    m = np.asarray(mascara_do_canal, dtype=bool)
    if not m.any():
        return ""
    indices = np.flatnonzero(m)
    i = int(indices[int(p[indices].argmax())])
    r = rotulos[i]
    partes = [x for x in ((r.get("chembl_id") or "").strip(),
                          (r.get("nome") or "").strip()) if x]
    organismo = (r.get("organismo") or "").strip()
    texto = " ".join(partes)
    return f"{texto} ({organismo})" if organismo else texto
