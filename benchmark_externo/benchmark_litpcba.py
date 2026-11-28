"""Benchmark do PP4CKD em LIT-PCBA, contra a linha de base 1-NN Tanimoto.

Protocolo
---------
Para cada um dos 15 conjuntos do LIT-PCBA (conjunto 'full'):

  * PP4CKD  -- escore = sigmoide(logit) na cabeca do alvo ChEMBL mapeado.
               Zero-shot: o modelo nunca viu rotulo do LIT-PCBA. (Mas viu o
               ChEMBL, que se sobrepoe -- e a ressalva de vazamento.)
  * 1-NN    -- maior Tanimoto ECFP4 aos ATIVOS do proprio alvo, leave-one-out
               (um ativo nao conta como vizinho de si mesmo). E a linha de base
               que Chen et al. 2019 mostram resolver quase todo o DUD-E, e que
               o PPB2 mostra ter a melhor recall em CV.

Para o cara-a-cara, os inativos sao subamostrados uniformemente (SUBAMOSTRA por
alvo). Subamostra uniforme nao enviesa o EF em esperanca, e e declarada.
O PP4CKD tambem roda sobre o conjunto COMPLETO, que e o comparavel com os
numeros publicados.

Saida: TSV com uma linha por (alvo, metodo, conjunto).
"""
import csv
import random
import sys
import time
from pathlib import Path

import caminhos_dados as cd_

REPO = cd_.REPO
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).parent))

LITPCBA = cd_.de_ambiente("LITPCBA_FULL")
CKPT = str(cd_.de_ambiente("CHECKPOINT_PP4CKD"))
ROTULOS = str(cd_.de_ambiente("ROTULOS_PP4CKD"))
SAIDA = Path(__file__).parent / "resultado_litpcba.tsv"

SUBAMOSTRA = 60_000      # inativos por alvo no cara-a-cara
LOTE = 2048
SEMENTE = 20261127

# LIT-PCBA -> alvo ChEMBL no espaco do modelo. ESR1_ago e ESR1_ant mapeiam
# para O MESMO alvo: o modelo tem uma saida para ESR1 e nao distingue
# agonista de antagonista. Isso e medido, nao contornado.
MAPA = {
    "ADRB2":    "CHEMBL210",
    "ALDH1":    "CHEMBL3577",
    "ESR1_ago": "CHEMBL206",
    "ESR1_ant": "CHEMBL206",
    "FEN1":     "CHEMBL5027",
    "GBA":      "CHEMBL2179",
    "IDH1":     "CHEMBL2007625",
    "KAT2A":    "CHEMBL5501",
    "MAPK1":    "CHEMBL4040",
    "MTORC1":   "CHEMBL2842",
    "OPRK1":    "CHEMBL237",
    "PKM2":     "CHEMBL1075189",
    "PPARG":    "CHEMBL235",
    "TP53":     "CHEMBL4096",
    "VDR":      "CHEMBL1977",
}

import numpy as np
import torch
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator, DataStructs
RDLogger.DisableLog("rdApp.*")

from mods import modelo as pa
import metricas_vs as mv


def ler_smi(caminho):
    """LIT-PCBA .smi: 'SMILES id' por linha."""
    saida = []
    with open(caminho) as fh:
        for linha in fh:
            partes = linha.split()
            if partes:
                saida.append(partes[0])
    return saida


def main():
    random.seed(SEMENTE)
    t0 = time.time()
    rede, cfg, _disp, nome_ckpt = pa.carregar_modelo(
        None, fp="ECFP4", checkpoint=CKPT, dispositivo="cpu")
    rotulos = pa.carregar_rotulos(ROTULOS)
    col = {r["chembl_id"]: i for i, r in enumerate(rotulos) if r["chembl_id"]}
    n_ativos_treino = {r["chembl_id"]: r["n_ativos"] for r in rotulos
                       if r["chembl_id"]}
    torch.set_num_threads(4)
    print(f"checkpoint {nome_ckpt}; {len(rotulos)} alvos", flush=True)

    gen = rdFingerprintGenerator.GetMorganGenerator(
        radius=cfg["raio"], fpSize=cfg["bits"])

    def pp4ckd_escores(smiles_lista, coluna):
        """Sigmoide do logit na coluna do alvo. SMILES ilegivel -> None."""
        out = [None] * len(smiles_lista)
        buf_i, buf_v = [], []

        def descarregar():
            if not buf_i:
                return
            x = torch.tensor(np.array(buf_v, dtype=np.float32))
            with torch.no_grad():
                p = torch.sigmoid(rede(x)[:, coluna]).numpy()
            for k, i in enumerate(buf_i):
                out[i] = float(p[k])
            buf_i.clear()
            buf_v.clear()

        for i, s in enumerate(smiles_lista):
            mol = Chem.MolFromSmiles(s)
            if mol is None:
                continue
            buf_i.append(i)
            buf_v.append(gen.GetFingerprintAsNumPy(mol))
            if len(buf_i) >= LOTE:
                descarregar()
        descarregar()
        return out

    def fps_de(smiles_lista):
        saida = []
        for s in smiles_lista:
            mol = Chem.MolFromSmiles(s)
            saida.append(gen.GetFingerprint(mol) if mol is not None else None)
        return saida

    linhas = []
    for alvo in sorted(MAPA):
        cid = MAPA[alvo]
        coluna = col[cid]
        d = LITPCBA / alvo
        ativos = ler_smi(d / "actives.smi")
        inativos = ler_smi(d / "inactives.smi")
        print(f"\n=== {alvo} -> {cid} ({n_ativos_treino[cid]} ativos de "
              f"treino no modelo) | {len(ativos)} ativos, "
              f"{len(inativos)} inativos", flush=True)

        # ---- PP4CKD no conjunto COMPLETO (comparavel com o publicado) -----
        t = time.time()
        e_at = pp4ckd_escores(ativos, coluna)
        e_in = pp4ckd_escores(inativos, coluna)
        esc, rot = [], []
        for v in e_at:
            if v is not None:
                esc.append(v); rot.append(1)
        for v in e_in:
            if v is not None:
                esc.append(v); rot.append(0)
        r = mv.resumo(esc, rot)
        r.update(alvo=alvo, chembl=cid, metodo="PP4CKD", conjunto="completo",
                 segundos=round(time.time() - t, 1))
        linhas.append(r)
        print(f"  PP4CKD completo: n={r['n']} ativos={r['n_ativos']} "
              f"AUROC={r['auroc']:.3f} EF1%={r['ef_1%']:.2f} "
              f"BEDROC={r['bedroc_20.5']:.3f} ({r['segundos']}s)", flush=True)

        # ---- cara-a-cara com subamostra de inativos -----------------------
        if len(inativos) > SUBAMOSTRA:
            idx = random.sample(range(len(inativos)), SUBAMOSTRA)
        else:
            idx = list(range(len(inativos)))
        sub_in = [inativos[i] for i in idx]
        sub_e_in = [e_in[i] for i in idx]

        esc_s, rot_s = [], []
        for v in e_at:
            if v is not None:
                esc_s.append(v); rot_s.append(1)
        for v in sub_e_in:
            if v is not None:
                esc_s.append(v); rot_s.append(0)
        r = mv.resumo(esc_s, rot_s)
        r.update(alvo=alvo, chembl=cid, metodo="PP4CKD",
                 conjunto=f"sub{SUBAMOSTRA}", segundos=0)
        linhas.append(r)

        # 1-NN Tanimoto leave-one-out
        t = time.time()
        fp_at = fps_de(ativos)
        fp_at_ok = [f for f in fp_at if f is not None]
        fp_sub = fps_de(sub_in)
        esc_n, rot_n = [], []
        for i, f in enumerate(fp_at):
            if f is None:
                continue
            outros = [g for j, g in enumerate(fp_at)
                      if g is not None and j != i]
            if not outros:
                continue
            esc_n.append(max(DataStructs.BulkTanimotoSimilarity(f, outros)))
            rot_n.append(1)
        for f in fp_sub:
            if f is None:
                continue
            esc_n.append(max(DataStructs.BulkTanimotoSimilarity(f, fp_at_ok)))
            rot_n.append(0)
        r = mv.resumo(esc_n, rot_n)
        r.update(alvo=alvo, chembl=cid, metodo="1NN_Tanimoto_LOO",
                 conjunto=f"sub{SUBAMOSTRA}", segundos=round(time.time()-t, 1))
        linhas.append(r)
        print(f"  1-NN  sub:       AUROC={r['auroc']:.3f} "
              f"EF1%={r['ef_1%']:.2f} BEDROC={r['bedroc_20.5']:.3f} "
              f"({r['segundos']}s)", flush=True)

        campos = ["alvo", "chembl", "metodo", "conjunto", "n", "n_ativos",
                  "auroc", "ef_0.5%", "ef_1%", "ef_5%", "bedroc_20.5",
                  "bedroc_0.85", "segundos"]
        with open(SAIDA, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=campos, delimiter="\t")
            w.writeheader()
            w.writerows(linhas)

    print(f"\n{len(linhas)} linhas em {SAIDA}  "
          f"({(time.time()-t0)/60:.1f} min)", flush=True)


if __name__ == "__main__":
    main()
