# -*- coding: utf-8 -*-
"""
dominio.py

Quando a PREVISAO do modelo tem habilidade medida, e quando ela e azar.

O que foi medido
----------------
`benchmark_externo/`, adendo 8: holdout temporal em 145 alvos do painel
ChEMBL 28, 1.769 ativos publicados a partir de 2022 contra 9.780 decoys
bioativos (`actives-as-decoys`). O eixo que separa acerto de azar NAO e a
probabilidade do modelo -- e a distancia estrutural do composto aos ativos
JA CONHECIDOS do alvo, medida como `sim_top5`: a media das 5 maiores
similaridades Tanimoto ECFP4 aos ativos de referencia daquele alvo.

    sim_top5 <  0,20   AUROC 0,4987   azar exato
    sim_top5 >= 0,20   AUROC 0,6819

Por decil, a transicao e abrupta (`resultado_decis_similaridade.tsv`):
quatro dos cinco primeiros decis ficam indistinguiveis de 0,5 e dois ficam
ABAIXO; do sexto em diante, 0,65 a 0,71.

Comparacao com o que o projeto usava como proxy de dominio, no mesmo dado
(`resultado_proxy_dominio.tsv`, criterio: parte o estrato na mediana do proxy
e compara a AUROC do modelo nas duas metades):

| proxy                            | alta  | baixa | delta  |    z |
|----------------------------------|-------|-------|--------|------|
| **media das 5 maiores (sim_top5)** | 0,682 | 0,499 | +0,183 | +12,1 |
| similaridade maxima              | 0,676 | 0,533 | +0,144 |  +9,6 |
| numero de bits do ECFP4          | 0,640 | 0,514 | +0,127 |  +7,1 |
| probabilidade maxima entre cabecas | 0,654 | 0,590 | +0,065 |  +4,2 |

Por que ha um minimo de referencias
-----------------------------------
`sim_top5` e a media das 5 MAIORES similaridades. Com poucas referencias, as
5 maiores sao simplesmente as unicas que existem, e o valor cai por
construcao -- nao porque o composto seja distante do conhecido, mas porque
pouco e conhecido. Medido nos 9 alvos de dor do projeto contra os 107.123
ligantes (`resultado_posicao_biblioteca.tsv`):

| alvo       | referencias | fracao com sim_top5 < 0,20 |
|------------|-------------|-----------------------------|
| CHEMBL233  | 5.790       | 28,5%  |
| CHEMBL228  | 4.434       | 53,5%  |
| CHEMBL4794 | 3.085       | 45,4%  |
| CHEMBL222  | 2.856       | 61,5%  |
| CHEMBL1904 | 1.549       | 75,7%  |
| CHEMBL1919 |   373       | 96,4%  |
| CHEMBL2015 |   113       | 96,5%  |
| CHEMBL1845 |    67       | 98,8%  |
| CHEMBL3896 |    19       | 99,9%  |

O salto esta entre 1.549 e 373 referencias. Abaixo disso, `sim_top5 < 0,20`
nao quer dizer "fora do dominio", quer dizer "nao da para saber" -- e punir
o composto por isso rebaixaria o alvo inteiro por falta de bibliografia, nao
por falta de sinal. Entao esses alvos recebem `indeterminado`.

O valor de 500 e DECLARADO, nao calibrado: o dado so o constrange ao
intervalo aberto (373, 1.549). Escolhido o numero redondo mais baixo que
separa, para nao descartar alvo por folga.

Limites de validade
-------------------
* A calibracao de 0,20 vem de um painel de alvos MADUROS (>= 50 ativos
  publicados depois de 2022). Nao foi medida em alvo pouco explorado.
* Vem do modelo multitarefa do ChEMBL 28, nao do PP4CKD. Para o PP4CKD nao
  existe numero sem vazamento com dado do ChEMBL (adendo 5).
* `0,20` e `0,25` entram no manifesto da execucao. Nao sao constante de
  codigo a ser esquecida: sao parametro do ranqueamento.
"""

LIMITE_BANDA = 0.20
LIMITE_ANALOGO = 0.25
MIN_REFERENCIAS = 500

ANALOGO = "analogo"
BANDA = "banda"
FORA = "fora"
INDETERMINADO = "indeterminado"

# Ordem de confiabilidade para o desempate do ranking. Menor ordena antes.
#
# `indeterminado` fica ENTRE os dois medidos, e nao junto de `fora`, porque
# nao saber nao e a mesma coisa que saber que e azar. Rebaixar o
# indeterminado ao fundo puniria o alvo por falta de bibliografia.
_ORDEM = {ANALOGO: 0, BANDA: 0, INDETERMINADO: 1, FORA: 2}


def sim_top5(fingerprint, fingerprints_de_referencia, k: int = 5):
    """Media das `k` maiores similaridades Tanimoto, ou None sem referencia."""
    if fingerprint is None or not fingerprints_de_referencia:
        return None
    from rdkit import DataStructs
    sims = DataStructs.BulkTanimotoSimilarity(
        fingerprint, list(fingerprints_de_referencia))
    sims.sort(reverse=True)
    maiores = sims[:k]
    return sum(maiores) / len(maiores)


def zona(similaridade, n_referencias: int) -> str:
    """Em que zona de confiabilidade a previsao cai.

    `indeterminado` quando faltam referencias para o valor significar algo,
    ou quando nao houve similaridade a calcular.
    """
    if n_referencias < MIN_REFERENCIAS or similaridade is None:
        return INDETERMINADO
    if similaridade >= LIMITE_ANALOGO:
        return ANALOGO
    if similaridade >= LIMITE_BANDA:
        return BANDA
    return FORA


def ordem(zona_do_composto: str | None) -> int:
    """Posicao de desempate. Zona ausente ordena como `indeterminado`."""
    if not zona_do_composto:
        return _ORDEM[INDETERMINADO]
    return _ORDEM.get(zona_do_composto, _ORDEM[INDETERMINADO])


def parametros() -> dict:
    """O que vai para o manifesto da execucao."""
    return {"limite_banda": LIMITE_BANDA,
            "limite_analogo": LIMITE_ANALOGO,
            "minimo_referencias": MIN_REFERENCIAS,
            "k_vizinhos": 5}
