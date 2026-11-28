"""Holdout TEMPORAL de verdade, nos alvos de dor do projeto.

O problema que este script corrige
----------------------------------
`calibrar_dor.py` dividiu os ativos do ChEMBL por ano e mediu o PP4CKD. Deu
AUROC 1,0 -- memorizacao, nao desempenho: o PP4CKD foi treinado em ChEMBL 36
(2025), entao todo "ativo de teste" publicado depois de 2015 esta dentro do
treino dele. A divisao temporal protegia o canal de similaridade e nao o
modelo.

A correcao
----------
Usa um modelo de ERA ANTERIOR aos ativos de teste: o
`chembl_28_multitask.onnx` oficial do EMBL-EBI, de ChEMBL 28 (jan/2021), e
testa em ativos cujo PRIMEIRO documento e de 2022 ou depois. Esses compostos
nao existiam no treino do modelo. E a pergunta do reposicionamento na sua
forma honesta: o que o modelo de 2021 teria apontado antes de alguem medir?

Cobertura: dos 9 alvos de dor, o painel de 616 alvos do ChEMBL 28 tem 3 --
SERT, NET e receptor mu. Sao os tres mais ricos em dado.

Decoys: ativos do LIT-PCBA (bioativos, drug-like, alvos sem relacao com dor)
que nao tenham medida no alvo em questao -- o controle `actives-as-decoys`.
"""
import csv
import json
import math
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

AQUI = Path(__file__).parent
CACHE = AQUI / "cache_chembl_dor"
LITPCBA = Path("/home/user/litpcba_download/full/LIT-PCBA_full")
MT28 = AQUI / "mt/trained_models/chembl_28_model"
SAIDA = AQUI / "resultado_calibracao_temporal.tsv"

ANO_REF_MAX = 2020      # referencia: ate aqui (ChEMBL 28 e de jan/2021)
ANO_TESTE_MIN = 2022    # teste: daqui em diante, fora do treino do modelo
CORTE_NOVIDADE = 0.25
FP_SIZE, RADIUS = 1024, 2

ALVOS = {
    "CHEMBL228": "Transportador de serotonina (SERT)",
    "CHEMBL222": "Transportador de noradrenalina (NET)",
    "CHEMBL233": "Receptor opioide mu",
}

import numpy as np
import onnxruntime
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors, rdFingerprintGenerator, DataStructs
RDLogger.DisableLog("rdApp.*")

import metricas_vs as mv
from calibrar_dor import buscar_atividades


def main():
    t0 = time.time()
    alvos_mt = json.load(open(MT28 / "targets_28.json"))
    col = {c: i for i, c in enumerate(alvos_mt)}
    opc = onnxruntime.SessionOptions()
    opc.intra_op_num_threads = 4
    ses = onnxruntime.InferenceSession(str(MT28 / "chembl_28_multitask.onnx"),
                                       sess_options=opc,
                                       providers=["CPUExecutionProvider"])
    nin = ses.get_inputs()[0].name
    print(f"modelo: chembl_28_multitask.onnx ({len(alvos_mt)} alvos)")
    print(f"referencia <= {ANO_REF_MAX} | teste >= {ANO_TESTE_MIN}\n",
          flush=True)

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=RADIUS,
                                                    fpSize=FP_SIZE)

    decoys = {}
    for d in sorted(LITPCBA.iterdir()):
        f = d / "actives.smi"
        if f.is_file():
            for linha in open(f):
                p = linha.split()
                if len(p) >= 2:
                    decoys[p[1]] = p[0]
    print(f"decoys disponiveis: {len(decoys)}\n", flush=True)

    def pontuar(smiles_lista, cid):
        """(probabilidades, fingerprints) alinhados; pula SMILES ilegivel."""
        ps, fps = [], []
        for s in smiles_lista:
            mol = Chem.MolFromSmiles(s)
            if mol is None:
                continue
            f = np.array(rdMolDescriptors.GetMorganFingerprintAsBitVect(
                mol, RADIUS, nBits=FP_SIZE), dtype=np.float32)
            ps.append(float(ses.run([cid], {nin: f})[0][0]))
            fps.append(gen.GetFingerprint(mol))
        return np.array(ps), fps

    linhas, ag = [], {"p": [], "r": []}
    for alvo, nome in ALVOS.items():
        print(f"=== {alvo}  {nome}", flush=True)
        ativ = buscar_atividades(alvo)
        por_mol = {}
        for a in ativ:
            if a.get("assay_type") not in ("B", "F"):
                continue
            sm, mid, ano = (a.get("canonical_smiles"),
                            a.get("molecule_chembl_id"), a.get("document_year"))
            if not (sm and mid and ano):
                continue
            ano = int(ano)
            if mid not in por_mol or ano < por_mol[mid][1]:
                por_mol[mid] = (sm, ano)
        ref = {m: v[0] for m, v in por_mol.items() if v[1] <= ANO_REF_MAX}
        tes = {m: v[0] for m, v in por_mol.items() if v[1] >= ANO_TESTE_MIN}
        dec = {m: s for m, s in decoys.items() if m not in por_mol}
        print(f"    referencia={len(ref)}  teste={len(tes)}  decoys={len(dec)}",
              flush=True)
        if len(tes) < 10:
            print("    teste pequeno demais, pulado\n", flush=True)
            continue

        cid = alvo
        if cid not in col:
            print("    fora do painel de 616, pulado\n", flush=True)
            continue
        _p, fps_ref = pontuar(list(ref.values()), cid)
        p_tes, fps_tes = pontuar(list(tes.values()), cid)
        p_dec, fps_dec = pontuar(list(dec.values()), cid)

        p_all = np.concatenate([p_tes, p_dec])
        rot = np.array([1] * len(p_tes) + [0] * len(p_dec))
        sim = np.array([max(DataStructs.BulkTanimotoSimilarity(f, fps_ref))
                        for f in (fps_tes + fps_dec)])

        for marca, sel in (("tudo", np.ones(len(rot), bool)),
                           ("novidade", sim < CORTE_NOVIDADE)):
            na, ni = int(rot[sel].sum()), int((1 - rot[sel]).sum())
            if na < 5 or ni < 50:
                print(f"    [{marca}] {na} ativos / {ni} decoys -- "
                      f"insuficiente", flush=True)
                continue
            res = mv.resumo(p_all[sel].tolist(), rot[sel].tolist())
            pares = sorted(zip(p_all[sel].tolist(), rot[sel].tolist()),
                           key=lambda x: (-x[0], x[1]))
            p20 = sum(r for _, r in pares[:20]) / 20
            base = na / (na + ni)
            linhas.append(dict(
                alvo=alvo, nome=nome, estrato=marca, n=na + ni, ativos=na,
                taxa_base=round(base, 5), auroc=round(res["auroc"], 4),
                ef_1=round(res["ef_1%"], 3), prec_20=round(p20, 4),
                ganho_20=round(p20 / base, 1) if base else ""))
            print(f"    [{marca}] n={na+ni} ativos={na} "
                  f"base={base*100:.2f}% AUROC={res['auroc']:.3f} "
                  f"EF1%={res['ef_1%']:.2f} prec@20={p20*100:.0f}% "
                  f"({p20/base:.0f}x)", flush=True)
            if marca == "novidade":
                v = p_all[sel]
                o = np.argsort(v, kind="mergesort")
                r = np.empty(len(v)); r[o] = np.arange(len(v), dtype=float)
                ag["p"].extend((r / max(1, len(v) - 1)).tolist())
                ag["r"].extend(rot[sel].tolist())
        print(flush=True)

        with open(SAIDA, "w", newline="") as fh:
            campos = ["alvo", "nome", "estrato", "n", "ativos", "taxa_base",
                      "auroc", "ef_1", "prec_20", "ganho_20"]
            w = csv.DictWriter(fh, fieldnames=campos, delimiter="\t")
            w.writeheader(); w.writerows(linhas)

    if sum(ag["r"]) >= 5:
        n, na = len(ag["r"]), int(sum(ag["r"]))
        a = mv.auroc(ag["p"], ag["r"])
        q1, q2 = a / (2 - a), 2 * a * a / (1 + a)
        nn = n - na
        se = math.sqrt((a * (1 - a) + (na - 1) * (q1 - a**2)
                        + (nn - 1) * (q2 - a**2)) / (na * nn))
        pares = sorted(zip(ag["p"], ag["r"]), key=lambda x: (-x[0], x[1]))
        print("=== AGREGADO, holdout temporal, estrato de novidade ===")
        print(f"n = {n:,}  ativos = {na}  taxa-base = {na/n*100:.3f}%")
        print(f"AUROC = {a:.4f}  EP {se:.4f}  "
              f"IC95% [{a-1.96*se:.4f}, {a+1.96*se:.4f}]  z = {(a-0.5)/se:.2f}")
        print(f"p (uma cauda) = {math.erfc((a-0.5)/se/math.sqrt(2))/2:.3e}")
        for k in (20, 50, 100):
            ac = sum(r for _, r in pares[:k])
            print(f"  precisao@{k:<4} = {ac}/{k} = {ac/k*100:5.2f}%  "
                  f"(ganho {ac/k/(na/n):.0f}x)")
    print(f"\n{len(linhas)} linhas em {SAIDA} ({(time.time()-t0)/60:.1f} min)")


if __name__ == "__main__":
    main()
