"""Distribuicao exata de Mann-Whitney para n1=n2=6 (o conjunto-controle).
Responde: uma AUC de 0,806 em 6 vs 6 e distinguivel do acaso?"""
from itertools import combinations
from fractions import Fraction

n1 = n2 = 6
N = n1 + n2
# Sob H0 toda atribuicao de 6 postos entre 12 e igualmente provavel.
todas = list(combinations(range(1, N + 1), n1))
total = len(todas)

def auc_de(postos):
    # U = soma dos postos do grupo 'trata' - n1(n1+1)/2 ; AUC = U/(n1*n2)
    U = sum(postos) - n1 * (n1 + 1) // 2
    return Fraction(U, n1 * n2)

aucs = sorted(auc_de(p) for p in todas)
print(f"n1={n1} n2={n2}  -> {total} atribuicoes, AUC assume "
      f"{len(set(aucs))} valores distintos, passo = {Fraction(1, n1*n2)}")

def p_uma_cauda(a):
    a = Fraction(a).limit_denominator(36)
    return sum(1 for x in aucs if x >= a) / total

for alvo, rotulo in [(Fraction(27,36), "0,750  (so alvos humanos)"),
                     (Fraction(29,36), "0,806  (p_humano x (1-p_bact))"),
                     (Fraction(36,36), "1,000  (CHEMBL360 isolado)"),
                     (Fraction(33,36), "0,917"),
                     (Fraction(31,36), "0,861")]:
    print(f"  AUC = {rotulo:<30} p (uma cauda) = {p_uma_cauda(alvo):.4f}")

alfa = 0.05
limiar = min(x for x in sorted(set(aucs)) if p_uma_cauda(x) <= alfa)
print(f"\nMenor AUC que atinge p <= {alfa} com 6 vs 6: {float(limiar):.3f} "
      f"({limiar})")
print(f"AUC = 1,000 tem p = {p_uma_cauda(Fraction(1)):.5f} "
      f"(= 1/{total}, o minimo possivel neste desenho)")
