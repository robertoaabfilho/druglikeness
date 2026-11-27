"""Metricas de triagem virtual, como a literatura de LIT-PCBA as reporta.

EF e BEDROC sao as que contam em LIT-PCBA; a AUROC fica em 51-73% e quase
nao separa metodos (argumento dos proprios autores do DrugCLIP).
"""
import math


def _ordenar(escores, rotulos):
    """Decrescente por escore. Empate resolvido pelo rotulo MENOR primeiro,
    isto e, de forma pessimista para o metodo -- empates nao viram credito."""
    pares = sorted(zip(escores, rotulos), key=lambda p: (-p[0], p[1]))
    return [r for _, r in pares]


def auroc(escores, rotulos):
    """AUC = P(ativo > inativo), empate = 0,5. Calculo por postos."""
    pares = sorted(zip(escores, rotulos), key=lambda p: p[0])
    n = len(pares)
    postos = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and pares[j + 1][0] == pares[i][0]:
            j += 1
        medio = (i + j) / 2.0 + 1.0          # posto medio, base 1
        for k in range(i, j + 1):
            postos[k] = medio
        i = j + 1
    n_pos = sum(1 for _, r in pares if r == 1)
    n_neg = n - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    soma_pos = sum(p for p, (_, r) in zip(postos, pares) if r == 1)
    return (soma_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def fator_enriquecimento(escores, rotulos, fracao):
    """EF a uma fracao do topo (0,01 = EF1%).

    EF = (ativos no topo / n_topo) / (ativos totais / n_total).
    """
    ordem = _ordenar(escores, rotulos)
    n = len(ordem)
    n_topo = max(1, int(round(fracao * n)))
    ativos_topo = sum(ordem[:n_topo])
    n_ativos = sum(ordem)
    if n_ativos == 0:
        return float("nan")
    return (ativos_topo / n_topo) / (n_ativos / n)


def bedroc(escores, rotulos, alfa=20.5):
    """BEDROC de Truchon & Bayly (2007). alfa=20,5 pesa o topo 8%;
    alfa=0,85 pesa o topo ~80% (os dois aparecem na literatura de LIT-PCBA)."""
    ordem = _ordenar(escores, rotulos)
    n = len(ordem)
    n_ativos = sum(ordem)
    if n_ativos == 0 or n_ativos == n:
        return float("nan")
    ra = n_ativos / n
    # soma sobre as posicoes (base 1) dos ativos
    soma = sum(math.exp(-alfa * (i + 1) / n)
               for i, r in enumerate(ordem) if r == 1)
    rie_num = soma / n_ativos
    rie_den = ((1 - math.exp(-alfa)) /
               (n * (math.exp(alfa / n) - 1)))
    rie = rie_num / rie_den
    return _bedroc_final(rie, ra, alfa)


def _bedroc_final(rie, ra, alfa):
    """BEDROC = RIE * (Ra*senh(a/2) / (cosh(a/2) - cosh(a/2 - a*Ra)))
                + 1/(1 - e^{a(1-Ra)})"""
    termo = (ra * math.sinh(alfa / 2) /
             (math.cosh(alfa / 2) - math.cosh(alfa / 2 - alfa * ra)))
    return rie * termo + 1.0 / (1.0 - math.exp(alfa * (1 - ra)))


def resumo(escores, rotulos):
    return {
        "n": len(rotulos),
        "n_ativos": sum(rotulos),
        "auroc": auroc(escores, rotulos),
        "ef_0.5%": fator_enriquecimento(escores, rotulos, 0.005),
        "ef_1%": fator_enriquecimento(escores, rotulos, 0.01),
        "ef_5%": fator_enriquecimento(escores, rotulos, 0.05),
        "bedroc_20.5": bedroc(escores, rotulos, 20.5),
        "bedroc_0.85": bedroc(escores, rotulos, 0.85),
    }


if __name__ == "__main__":
    # Sanidade: ranking perfeito, aleatorio e invertido.
    import random
    random.seed(0)
    n, na = 10000, 100
    rot = [1] * na + [0] * (n - na)
    perfeito = [1.0] * na + [0.0] * (n - na)
    invertido = [0.0] * na + [1.0] * (n - na)
    aleat = [random.random() for _ in range(n)]
    for nome, esc in [("perfeito", perfeito), ("aleatorio", aleat),
                      ("invertido", invertido)]:
        r = resumo(esc, rot)
        print(f"{nome:<10} auroc={r['auroc']:.3f} ef1%={r['ef_1%']:.2f} "
              f"bedroc(20.5)={r['bedroc_20.5']:.3f}")
    print("esperado:  perfeito auroc~1 ef1%=100 bedroc~1 | "
          "aleatorio auroc~0,5 ef1%~1 bedroc~0,08 | invertido auroc~0 ef1%=0")
