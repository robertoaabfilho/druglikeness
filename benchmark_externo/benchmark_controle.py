"""Benchmark do conjunto-controle: reproduz as AUCs do README e acrescenta as
linhas de base triviais que faltam.

Roda do diretorio do repositorio:
    python3 benchmark_controle.py
"""
import sys
from itertools import combinations
from pathlib import Path

REPO = Path("/home/user/druglikeness")
sys.path.insert(0, str(REPO))

CKPT = "/home/user/pp4ckd_download/extraido/models/ecfp4_fold_4.pt"
ROTULOS = "/home/user/jixing475/pp4ckd/data/target_labels.tsv"

ALVOS_HUMANOS = ["CHEMBL1919", "CHEMBL3896", "CHEMBL228", "CHEMBL222",
                 "CHEMBL1845", "CHEMBL4794", "CHEMBL233", "CHEMBL2015",
                 "CHEMBL1904"]
MTB = "CHEMBL360"          # M. tuberculosis, celula inteira
SAUREUS = "CHEMBL352"      # S. aureus, celula inteira

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator, DataStructs
RDLogger.DisableLog("rdApp.*")

from mods import modelo as pa


def ler_controle(caminho):
    compostos = []
    with open(caminho) as fh:
        for linha in fh:
            if linha.startswith("#") or not linha.strip():
                continue
            campos = linha.rstrip("\n").split("\t")
            compostos.append({"nome": campos[0], "smiles": campos[1],
                              "grupo": campos[3]})
    return compostos


def auc(pos, neg):
    """AUC = P(score_pos > score_neg), empates contam meio."""
    n = 0.0
    for p in pos:
        for q in neg:
            n += 1.0 if p > q else (0.5 if p == q else 0.0)
    return n / (len(pos) * len(neg))


_NULA = None


def p_exato(valor_auc, n1, n2):
    """p de uma cauda pela distribuicao exata de Mann-Whitney."""
    global _NULA
    chave = (n1, n2)
    if _NULA is None or _NULA[0] != chave:
        todas = list(combinations(range(1, n1 + n2 + 1), n1))
        aucs = sorted((sum(c) - n1 * (n1 + 1) / 2) / (n1 * n2) for c in todas)
        _NULA = (chave, aucs)
    aucs = _NULA[1]
    tol = 1e-9
    return sum(1 for x in aucs if x >= valor_auc - tol) / len(aucs)


def main():
    compostos = ler_controle(REPO / "exemplos/controle.dor_neuropatica.tsv")
    trata = [c for c in compostos if c["grupo"] == "trata"]
    causa = [c for c in compostos if c["grupo"] == "causa"]
    print(f"conjunto-controle: {len(trata)} 'trata' x {len(causa)} 'causa'")
    print("  trata:", ", ".join(c["nome"] for c in trata))
    print("  causa:", ", ".join(c["nome"] for c in causa))

    rede, cfg, _disp, nome_ckpt = pa.carregar_modelo(
        None, fp="ECFP4", checkpoint=CKPT, dispositivo="cpu")
    rotulos = pa.carregar_rotulos(ROTULOS)
    idx = {r["chembl_id"]: i for i, r in enumerate(rotulos) if r["chembl_id"]}
    n_ativos = {r["chembl_id"]: r["n_ativos"] for r in rotulos if r["chembl_id"]}
    print(f"checkpoint: {nome_ckpt}  |  alvos no modelo: {len(rotulos)}")

    import torch
    probs = {}
    for c in compostos:
        fp = pa.fingerprint_de(c["smiles"], cfg)
        with torch.no_grad():
            v = rede(torch.tensor(fp, dtype=torch.float32).unsqueeze(0))
            # A rede devolve LOGITS. O sigmoide e obrigatorio: o termo
            # negativo (1 - p) so faz sentido com p em [0,1].
            v = torch.sigmoid(v)
        probs[c["nome"]] = v.squeeze(0).numpy()

    def p_de(nome, chembl):
        return float(probs[nome][int(idx[chembl])])

    # --- fingerprints para a linha de base de similaridade -----------------
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = {c["nome"]: gen.GetFingerprint(Chem.MolFromSmiles(c["smiles"]))
           for c in compostos}

    def sim_max_aos_trata(nome):
        """Maior Tanimoto aos 'trata'. Leave-one-out quando o proprio e trata."""
        outros = [t["nome"] for t in trata if t["nome"] != nome]
        return max(DataStructs.TanimotoSimilarity(fps[nome], fps[o])
                   for o in outros)

    # --- escores ------------------------------------------------------------
    def agg(nome, alvos, como):
        vals = [p_de(nome, a) for a in alvos]
        return {"max": max, "media": lambda v: sum(v) / len(v),
                "soma": sum}[como](vals)

    escores = {}
    for como in ("max", "media", "soma"):
        escores[f"PP4CKD alvos humanos ({como})"] = (
            lambda n, c=como: agg(n, ALVOS_HUMANOS, c))
    escores["PP4CKD humano(max) x (1 - p_Mtb)"] = (
        lambda n: agg(n, ALVOS_HUMANOS, "max") * (1 - p_de(n, MTB)))
    escores["PP4CKD humano(media) x (1 - p_Mtb)"] = (
        lambda n: agg(n, ALVOS_HUMANOS, "media") * (1 - p_de(n, MTB)))
    escores["PP4CKD humano(max) x (1 - p_Mtb) x (1 - p_Saureus)"] = (
        lambda n: agg(n, ALVOS_HUMANOS, "max") * (1 - p_de(n, MTB))
        * (1 - p_de(n, SAUREUS)))
    escores["-p_Mtb isolado (CHEMBL360)"] = lambda n: -p_de(n, MTB)
    # linhas de base que o repositorio nao tem
    escores["BASE similaridade Tanimoto aos 'trata'"] = sim_max_aos_trata
    escores["BASE promiscuidade (soma dos 7.676)"] = (
        lambda n: float(probs[n].sum()))
    escores["BASE -promiscuidade"] = lambda n: -float(probs[n].sum())

    print(f"\n{'escore':<52} {'AUC':>6} {'p exato':>9}")
    print("-" * 70)
    linhas = []
    for nome_escore, fn in escores.items():
        sp = [fn(c["nome"]) for c in trata]
        sn = [fn(c["nome"]) for c in causa]
        a = auc(sp, sn)
        linhas.append((a, nome_escore, p_exato(a, len(trata), len(causa))))
    for a, nome_escore, p in sorted(linhas, reverse=True):
        marca = " *" if p <= 0.05 else ""
        print(f"{nome_escore:<52} {a:6.3f} {p:9.4f}{marca}")

    print("\n* p <= 0,05 (uma cauda, Mann-Whitney exato, 6 vs 6)")
    print("Nota: com 6 vs 6 a AUC so assume multiplos de 1/36 = 0,0278,")
    print("      e 0,806 e o menor valor que alcanca p <= 0,05.")

    print(f"\n{'composto':<16} {'grupo':<6} {'p_hum(max)':>11} "
          f"{'p_Mtb':>8} {'sim_trata':>10} {'promisc':>9}")
    for c in sorted(compostos, key=lambda x: x["grupo"]):
        n = c["nome"]
        print(f"{n:<16} {c['grupo']:<6} {agg(n, ALVOS_HUMANOS, 'max'):11.4f} "
              f"{p_de(n, MTB):8.4f} {sim_max_aos_trata(n):10.3f} "
              f"{probs[n].sum():9.2f}")


if __name__ == "__main__":
    main()
