"""O canal predito enriquece em compostos DISTANTES dos ativos conhecidos?

E a pergunta do reposicionamento. Se os acertos do PP4CKD estao todos na
vizinhanca quimica dos ativos conhecidos, ele nao acrescenta nada ao 1-NN e
nao serve para reposicionar. Se ele acerta a baixa similaridade, ai esta a
novidade.

Desenho: LIT-PCBA AVE-unbiased. Para cada molecula do conjunto V calcula
`sim_T` = maior Tanimoto aos ativos de TREINO, e estratifica por essa
similaridade. DENTRO de cada estrato -- isto e, comparando ativos e inativos
igualmente distantes do treino -- mede se `p_alvo` ainda separa.

Estratificar por sim_T e o controle de VIES DE ANALOGIA no espirito de
Chen et al. 2019: no estrato de baixa similaridade, um modelo que so memoriza
vizinhanca cai para AUROC 0,5.

Agregacao entre alvos: o escore e convertido em posto percentual DENTRO de
cada alvo antes de juntar, para que alvos com distribuicoes diferentes nao se
somem indevidamente.
"""
import csv
import sys
import time
from pathlib import Path

REPO = Path("/home/user/druglikeness")
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).parent))

AVE = Path("/home/user/litpcba_download/ave/LIT-PCBA_AVE_unbiased")
CKPT = "/home/user/pp4ckd_download/extraido/models/ecfp4_fold_4.pt"
ROTULOS = "/home/user/jixing475/pp4ckd/data/target_labels.tsv"
SAIDA = Path(__file__).parent / "resultado_novidade.tsv"

LOTE = 2048
MIN_ATIVOS_V = 15
FAIXAS = [(0.0, 0.25), (0.25, 0.35), (0.35, 0.45), (0.45, 0.60), (0.60, 1.01)]

MAPA = {
    "ADRB2": "CHEMBL210", "ALDH1": "CHEMBL3577", "ESR1_ago": "CHEMBL206",
    "ESR1_ant": "CHEMBL206", "FEN1": "CHEMBL5027", "GBA": "CHEMBL2179",
    "IDH1": "CHEMBL2007625", "KAT2A": "CHEMBL5501", "MAPK1": "CHEMBL4040",
    "MTORC1": "CHEMBL2842", "OPRK1": "CHEMBL237", "PKM2": "CHEMBL1075189",
    "PPARG": "CHEMBL235", "TP53": "CHEMBL4096", "VDR": "CHEMBL1977",
}

import numpy as np
import torch
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator, DataStructs
RDLogger.DisableLog("rdApp.*")

from mods import modelo as pa
import metricas_vs as mv


def ler_smi(c):
    out = []
    with open(c) as fh:
        for l in fh:
            p = l.split()
            if p:
                out.append(p[0])
    return out


def postos(v):
    v = np.asarray(v, float)
    if len(v) < 2:
        return np.zeros(len(v))
    o = np.argsort(v, kind="mergesort")
    r = np.empty(len(v), float)
    r[o] = np.arange(len(v), dtype=float)
    return r / (len(v) - 1)


def main():
    t0 = time.time()
    rede, cfg, _d, ck = pa.carregar_modelo(None, fp="ECFP4", checkpoint=CKPT,
                                           dispositivo="cpu")
    rot_mod = pa.carregar_rotulos(ROTULOS)
    col = {r["chembl_id"]: i for i, r in enumerate(rot_mod) if r["chembl_id"]}
    torch.set_num_threads(4)
    gen = rdFingerprintGenerator.GetMorganGenerator(
        radius=cfg["raio"], fpSize=cfg["bits"])
    print(f"checkpoint {ck}", flush=True)

    # acumuladores agregados (postos percentuais por alvo)
    agreg = {f: {"p": [], "s": [], "r": []} for f in FAIXAS}
    linhas = []

    for alvo in sorted(MAPA):
        cid = MAPA[alvo]
        d = AVE / alvo
        at_T, at_V, in_V = (ler_smi(d / f"{n}.smi")
                            for n in ("active_T", "active_V", "inactive_V"))
        if len(at_V) < MIN_ATIVOS_V:
            continue
        coluna = col[cid]
        t = time.time()
        fps_T = [gen.GetFingerprint(m) for m in
                 (Chem.MolFromSmiles(s) for s in at_T) if m is not None]

        p_list, fp_list, rotu = [], [], []
        bi, bv = [], []

        def descarregar():
            if not bi:
                return
            x = torch.tensor(np.array(bv, dtype=np.float32))
            with torch.no_grad():
                p = torch.sigmoid(rede(x)[:, coluna]).numpy()
            p_list.extend(p.tolist())
            bi.clear(); bv.clear()

        for lista, marca in ((at_V, 1), (in_V, 0)):
            for s in lista:
                mol = Chem.MolFromSmiles(s)
                if mol is None:
                    continue
                fp_list.append(gen.GetFingerprint(mol))
                rotu.append(marca)
                bi.append(1)
                bv.append(gen.GetFingerprintAsNumPy(mol))
                if len(bi) >= LOTE:
                    descarregar()
        descarregar()

        p_alvo = np.array(p_list)
        rotu = np.array(rotu)
        sim_T = np.array([max(DataStructs.BulkTanimotoSimilarity(f, fps_T))
                          for f in fp_list])
        rp, rs = postos(p_alvo), postos(sim_T)

        for lo, hi in FAIXAS:
            sel = (sim_T >= lo) & (sim_T < hi)
            na, ni = int(rotu[sel].sum()), int((1 - rotu[sel]).sum())
            agreg[(lo, hi)]["p"].extend(rp[sel].tolist())
            agreg[(lo, hi)]["s"].extend(rs[sel].tolist())
            agreg[(lo, hi)]["r"].extend(rotu[sel].tolist())
            if na >= 3 and ni >= 20:
                res = mv.resumo(p_alvo[sel].tolist(), rotu[sel].tolist())
                linhas.append(dict(
                    alvo=alvo, faixa=f"{lo:.2f}-{hi:.2f}", n=na + ni,
                    ativos=na, inativos=ni,
                    taxa_base=round(na / (na + ni), 5),
                    auroc_p=round(res["auroc"], 4),
                    ef1_p=round(res["ef_1%"], 3),
                    auroc_s=round(mv.auroc(sim_T[sel].tolist(),
                                           rotu[sel].tolist()), 4)))
        print(f"=== {alvo}: {len(at_V)} ativos V, {len(fps_T)} ativos T "
              f"({time.time()-t:.0f}s)", flush=True)

    with open(SAIDA, "w", newline="") as fh:
        campos = ["alvo", "faixa", "n", "ativos", "inativos", "taxa_base",
                  "auroc_p", "ef1_p", "auroc_s"]
        w = csv.DictWriter(fh, fieldnames=campos, delimiter="\t")
        w.writeheader(); w.writerows(linhas)

    print("\n=== AGREGADO (postos percentuais dentro de cada alvo) ===")
    print(f"{'faixa sim_T':<14} {'n':>8} {'ativos':>7} {'taxa-base':>10} "
          f"{'AUROC p_alvo':>13} {'AUROC sim_T':>12}")
    print("-" * 70)
    for f in FAIXAS:
        a = agreg[f]
        if sum(a["r"]) < 3:
            continue
        n = len(a["r"]); na = int(sum(a["r"]))
        ap = mv.auroc(a["p"], a["r"])
        as_ = mv.auroc(a["s"], a["r"])
        print(f"{f[0]:.2f}-{f[1]:.2f}      {n:>8} {na:>7} "
              f"{na/n*100:9.3f}% {ap:13.3f} {as_:12.3f}")
    print(f"\n{len(linhas)} linhas em {SAIDA} ({(time.time()-t0)/60:.1f} min)")
    print("\nLeitura: na faixa de BAIXA similaridade, AUROC do p_alvo acima de")
    print("0,5 e novidade real -- o modelo acerta onde a vizinhanca nao ajuda.")
    print("Se cair a 0,5, o canal predito nao serve para reposicionamento.")


if __name__ == "__main__":
    main()
