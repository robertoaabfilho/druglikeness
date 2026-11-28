"""Descritor (QSAR/fisico-quimico) como espaco de similaridade ALTERNATIVO,
dentro do estrato de novidade.

Tres perguntas, todas no estrato `sim_T < 0.25` (ECFP4) do LIT-PCBA
AVE-unbiased:

1. O canal de descritor (`sim_desc`, proximidade aos ativos de treino no
   espaco de descritores padronizado) separa ativo de inativo ali?
2. Ele ACRESCENTA ao `p_alvo` quando fundido? (A fusao foi medida no agregado,
   no Adendo 2, e nao neste estrato, onde a linha de base e outra.)
3. Quanto os dois espacos DISCORDAM sobre o que e novo? Um composto distante
   em ECFP4 e vizinho em descritor nao e novidade no sentido que interessa.

Descritores: conjunto padrao do RDKit, calculado localmente. O preditor de
ADMET por servico externo e API
externa e rejeita 30-35% das moleculas, entao nao entra aqui; o bloco
fisico-quimico e o que o projeto chama de druglikeness e e computavel.
Padronizacao por z-score sobre o proprio conjunto de triagem.
"""
import csv
import sys
import time
from pathlib import Path

import caminhos_dados as cd_

REPO = cd_.REPO
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).parent))

AVE = cd_.de_ambiente("LITPCBA_AVE")
CKPT = str(cd_.de_ambiente("CHECKPOINT_PP4CKD"))
ROTULOS = str(cd_.de_ambiente("ROTULOS_PP4CKD"))
SAIDA = Path(__file__).parent / "resultado_descritor.tsv"

LOTE = 2048
MIN_ATIVOS_V = 15
CORTE_NOVIDADE = 0.25

MAPA = {
    "ALDH1": "CHEMBL3577", "ESR1_ant": "CHEMBL206", "FEN1": "CHEMBL5027",
    "GBA": "CHEMBL2179", "KAT2A": "CHEMBL5501", "MAPK1": "CHEMBL4040",
    "MTORC1": "CHEMBL2842", "PKM2": "CHEMBL1075189", "TP53": "CHEMBL4096",
    "VDR": "CHEMBL1977",
}

import numpy as np
import torch
from rdkit import Chem, RDLogger
from rdkit.Chem import (Crippen, Descriptors, QED, rdFingerprintGenerator,
                        rdMolDescriptors, DataStructs)
RDLogger.DisableLog("rdApp.*")

from mods import modelo as pa
import metricas_vs as mv


def descritores(mol):
    """12 descritores fisico-quimicos -- o bloco de druglikeness."""
    return [
        Descriptors.MolWt(mol),
        Crippen.MolLogP(mol),
        rdMolDescriptors.CalcTPSA(mol),
        rdMolDescriptors.CalcNumHBD(mol),
        rdMolDescriptors.CalcNumHBA(mol),
        rdMolDescriptors.CalcNumRotatableBonds(mol),
        rdMolDescriptors.CalcNumAromaticRings(mol),
        rdMolDescriptors.CalcFractionCSP3(mol),
        mol.GetNumHeavyAtoms(),
        QED.qed(mol),
        Crippen.MolMR(mol),
        rdMolDescriptors.CalcNumRings(mol),
    ]


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
    rm = pa.carregar_rotulos(ROTULOS)
    col = {r["chembl_id"]: i for i, r in enumerate(rm) if r["chembl_id"]}
    torch.set_num_threads(4)
    gen = rdFingerprintGenerator.GetMorganGenerator(
        radius=cfg["raio"], fpSize=cfg["bits"])
    print(f"checkpoint {ck}", flush=True)

    ag = {k: [] for k in ("p", "sd", "r", "fus_media", "fus_prod")}
    discord = {"ecfp_longe_desc_perto": 0, "ecfp_longe_desc_longe": 0}
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

        mols_T = [m for m in (Chem.MolFromSmiles(s) for s in at_T)
                  if m is not None]
        fps_T = [gen.GetFingerprint(m) for m in mols_T]
        desc_T = np.array([descritores(m) for m in mols_T])

        ps, fps, rotu, desc_V = [], [], [], []
        bv = []

        def descarregar():
            if not bv:
                return
            x = torch.tensor(np.array(bv, dtype=np.float32))
            with torch.no_grad():
                ps.extend(torch.sigmoid(rede(x)[:, coluna]).numpy().tolist())
            bv.clear()

        for lista, marca in ((at_V, 1), (in_V, 0)):
            for s in lista:
                mol = Chem.MolFromSmiles(s)
                if mol is None:
                    continue
                fps.append(gen.GetFingerprint(mol))
                desc_V.append(descritores(mol))
                rotu.append(marca)
                bv.append(gen.GetFingerprintAsNumPy(mol))
                if len(bv) >= LOTE:
                    descarregar()
        descarregar()

        p_alvo = np.array(ps)
        rotu = np.array(rotu)
        desc_V = np.array(desc_V)
        sim_T = np.array([max(DataStructs.BulkTanimotoSimilarity(f, fps_T))
                          for f in fps])

        # z-score pelo conjunto de triagem (inativos dominam, e o certo:
        # a escala e a da biblioteca que se esta triando)
        mu, sd = desc_V.mean(axis=0), desc_V.std(axis=0)
        sd[sd == 0] = 1.0
        zV = (desc_V - mu) / sd
        zT = (desc_T - mu) / sd
        # distancia euclidiana ao ativo de treino MAIS PROXIMO -> similaridade
        dist = np.sqrt(((zV[:, None, :] - zT[None, :, :]) ** 2).sum(axis=2)) \
            if len(zT) * len(zV) < 4e7 else None
        if dist is None:
            # em lotes, para nao estourar memoria
            menor = np.empty(len(zV))
            passo = max(1, int(4e7 // max(1, len(zT))))
            for i in range(0, len(zV), passo):
                bloco = zV[i:i + passo]
                dd = np.sqrt(((bloco[:, None, :] - zT[None, :, :]) ** 2)
                             .sum(axis=2))
                menor[i:i + passo] = dd.min(axis=1)
        else:
            menor = dist.min(axis=1)
        sim_desc = 1.0 / (1.0 + menor)

        novo = sim_T < CORTE_NOVIDADE
        # discordancia entre os dois espacos, dentro do estrato novo
        corte_desc = np.median(sim_desc)
        discord["ecfp_longe_desc_perto"] += int(
            (novo & (sim_desc >= corte_desc)).sum())
        discord["ecfp_longe_desc_longe"] += int(
            (novo & (sim_desc < corte_desc)).sum())

        na, ni = int(rotu[novo].sum()), int((1 - rotu[novo]).sum())
        if na < 3 or ni < 20:
            print(f"=== {alvo}: estrato novo com {na} ativos, pulado",
                  flush=True)
            continue
        rp = postos(p_alvo[novo])
        rd = postos(sim_desc[novo])
        ag["p"].extend(rp.tolist())
        ag["sd"].extend(rd.tolist())
        ag["r"].extend(rotu[novo].tolist())
        ag["fus_media"].extend(((rp + rd) / 2).tolist())
        ag["fus_prod"].extend((rp * rd).tolist())

        linhas.append(dict(
            alvo=alvo, n=na + ni, ativos=na,
            taxa_base=round(na / (na + ni), 5),
            auroc_p=round(mv.auroc(rp.tolist(), rotu[novo].tolist()), 4),
            auroc_desc=round(mv.auroc(rd.tolist(), rotu[novo].tolist()), 4),
            auroc_fusao=round(mv.auroc(((rp + rd) / 2).tolist(),
                                       rotu[novo].tolist()), 4)))
        print(f"=== {alvo}: estrato novo n={na+ni}, {na} ativos "
              f"({time.time()-t:.0f}s)", flush=True)

    with open(SAIDA, "w", newline="") as fh:
        campos = ["alvo", "n", "ativos", "taxa_base", "auroc_p",
                  "auroc_desc", "auroc_fusao"]
        w = csv.DictWriter(fh, fieldnames=campos, delimiter="\t")
        w.writeheader(); w.writerows(linhas)

    import math
    print(f"\n=== AGREGADO, estrato de novidade (sim_T < {CORTE_NOVIDADE}) ===")
    n, na = len(ag["r"]), int(sum(ag["r"]))
    print(f"n = {n:,}   ativos = {na}   taxa-base = {na/n*100:.3f}%\n")
    print(f"{'canal':<22} {'AUROC':>8} {'EF 0,1%':>9} {'prec@20':>9}")
    print("-" * 52)
    for nome, chave in (("p_alvo (PP4CKD)", "p"),
                        ("sim_desc (descritor)", "sd"),
                        ("fusao media", "fus_media"),
                        ("fusao produto", "fus_prod")):
        e = ag[chave]
        a = mv.auroc(e, ag["r"])
        ef = mv.fator_enriquecimento(e, ag["r"], 0.001)
        pares = sorted(zip(e, ag["r"]), key=lambda p: (-p[0], p[1]))[:20]
        pk = sum(r for _, r in pares) / 20
        print(f"{nome:<22} {a:8.4f} {ef:9.2f} {pk*100:8.2f}%")
    tot = sum(discord.values())
    print(f"\nDiscordancia entre os dois espacos no estrato novo (n={tot:,}):")
    for k, v in discord.items():
        print(f"  {k:<26} {v:>9,}  ({v/tot*100:5.1f}%)")
    print(f"\n({(time.time()-t0)/60:.1f} min)")


if __name__ == "__main__":
    main()
