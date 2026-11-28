"""Holdout temporal ESCALADO: o painel ChEMBL 28 inteiro, nao 3 alvos.

Objetivo: apertar o intervalo de confianca do Adendo 5, que tinha 29 ativos no
estrato de novidade. Mesma receita, muitos mais alvos.

Desenho (identico ao Adendo 5, so maior)
----------------------------------------
modelo    chembl_28_multitask.onnx (ChEMBL 28, jan/2021, 616 alvos)
ativos    ChEMBL, pchembl >= 6, ensaio B ou F
          referencia = primeiro documento <= 2020
          teste      = primeiro documento >= 2022  (fora do treino do modelo)
decoys    os ativos do LIT-PCBA sem medida no alvo em questao
novidade  maior Tanimoto ECFP4 aos ativos de referencia < 0,25

Otimizacao que torna isso barato: o ONNX devolve as 616 saidas por chamada, e
o conjunto de decoys e o mesmo para todo alvo. Entao cada molecula e pontuada
UMA vez e a matriz (molecula x 616) fica em memoria; cada alvo so indexa a
sua coluna. Os ativos de referencia nao precisam do modelo -- so de
fingerprint.
"""
import csv
import json
import math
import sys
import time
from pathlib import Path

import caminhos_dados as cd_

AQUI = Path(__file__).parent
sys.path.insert(0, str(AQUI))

MT28 = AQUI / "mt/trained_models/chembl_28_model"
LITPCBA = cd_.de_ambiente("LITPCBA_FULL")
SONDAGEM = AQUI / "sondagem_alvos.json"
SAIDA = AQUI / "resultado_holdout_escalado.tsv"

ANO_REF_MAX, ANO_TESTE_MIN = 2020, 2022
CORTE_NOVIDADE = 0.25
PCHEMBL_MIN = 6.0
FP_SIZE, RADIUS = 1024, 2
MIN_RECENTES = 50        # ativos pos-2022 para o alvo valer a coleta
MAX_ALVOS = 150          # limita o tempo de coleta
MAX_ATIVIDADES = 15000   # pula alvo gigante: a coleta dele nao cabe no tempo
SEMENTE = 20261128       # amostra reproduzivel
MIN_ATIVOS_ESTRATO = 3   # minimo para o alvo entrar numa linha de resultado

import numpy as np
import onnxruntime
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors, rdFingerprintGenerator, DataStructs
RDLogger.DisableLog("rdApp.*")

import metricas_vs as mv
from calibrar_dor import buscar_atividades


def main():
    t0 = time.time()
    painel = json.load(open(MT28 / "targets_28.json"))
    col = {c: i for i, c in enumerate(painel)}
    sond = json.loads(SONDAGEM.read_text())
    elegiveis = sorted(a for a, v in sond.items()
                       if a in col and v.get("recentes")
                       and v["recentes"] >= MIN_RECENTES)
    # AMOSTRA ALEATORIA, nao os maiores. Pegar os maiores enviesaria para
    # alvo muito estudado e deixaria o agregado ser puxado por poucas
    # familias -- o defeito do Adendo 5, em que tudo vinha do receptor mu.
    import random
    rnd = random.Random(SEMENTE)
    escolhidos = (sorted(rnd.sample(elegiveis, MAX_ALVOS))
                  if len(elegiveis) > MAX_ALVOS else elegiveis)
    cand = [(a, sond[a]["recentes"]) for a in escolhidos]
    print(f"{len(elegiveis)} elegiveis -> amostra de {len(cand)} "
          f"(semente {SEMENTE})", flush=True)
    print(f"painel {len(painel)} alvos | {len(cand)} selecionados "
          f"(>= {MIN_RECENTES} ativos pos-{ANO_TESTE_MIN})", flush=True)
    print(f"soma de ativos recentes: {sum(c for _, c in cand):,}\n", flush=True)

    opc = onnxruntime.SessionOptions()
    opc.intra_op_num_threads = 4
    ses = onnxruntime.InferenceSession(str(MT28 / "chembl_28_multitask.onnx"),
                                       sess_options=opc,
                                       providers=["CPUExecutionProvider"])
    nin = ses.get_inputs()[0].name
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=RADIUS,
                                                    fpSize=FP_SIZE)

    # ---- decoys -----------------------------------------------------------
    decoys = {}
    for d in sorted(LITPCBA.iterdir()):
        f = d / "actives.smi"
        if f.is_file():
            for linha in open(f):
                p = linha.split()
                if len(p) >= 2:
                    decoys[p[1]] = p[0]
    print(f"decoys: {len(decoys)}", flush=True)

    # ---- coleta -----------------------------------------------------------
    from sondar_alvos import contar
    por_alvo = {}
    pulados_grandes = 0
    for i, (alvo, _n) in enumerate(cand):
        total = contar(alvo)
        if total and total > MAX_ATIVIDADES:
            pulados_grandes += 1
            continue
        try:
            ativ = buscar_atividades(alvo)
        except Exception as erro:
            print(f"  [{i+1}/{len(cand)}] {alvo}: coleta falhou "
                  f"({type(erro).__name__}), pulado", flush=True)
            continue
        pm = {}
        for a in ativ:
            if a.get("assay_type") not in ("B", "F"):
                continue
            sm, mid, ano = (a.get("canonical_smiles"),
                            a.get("molecule_chembl_id"),
                            a.get("document_year"))
            if not (sm and mid and ano):
                continue
            ano = int(ano)
            if mid not in pm or ano < pm[mid][1]:
                pm[mid] = (sm, ano)
        ref = {m: v[0] for m, v in pm.items() if v[1] <= ANO_REF_MAX}
        tes = {m: v[0] for m, v in pm.items() if v[1] >= ANO_TESTE_MIN}
        if len(tes) < 10 or len(ref) < 10:
            continue
        por_alvo[alvo] = {"ref": ref, "tes": tes, "medidos": set(pm)}
        if (i + 1) % 10 == 0:
            print(f"  [{i+1}/{len(cand)}] coletados {len(por_alvo)} alvos "
                  f"({(time.time()-t0)/60:.1f} min)", flush=True)
    print(f"\n{len(por_alvo)} alvos com dado suficiente; "
          f"{pulados_grandes} pulados por tamanho (> {MAX_ATIVIDADES} "
          f"atividades) ({(time.time()-t0)/60:.1f} min)\n", flush=True)

    # ---- universo de moleculas: pontuar UMA vez ---------------------------
    universo = dict(decoys)
    for v in por_alvo.values():
        universo.update(v["tes"])
    refs_todas = {}
    for v in por_alvo.values():
        refs_todas.update(v["ref"])
    print(f"universo a pontuar: {len(universo):,} moleculas; "
          f"referencias (so fingerprint): {len(refs_todas):,}", flush=True)

    fp_cache, probs = {}, {}
    t = time.time()
    for k, (mid, sm) in enumerate(universo.items()):
        mol = Chem.MolFromSmiles(sm)
        if mol is None:
            continue
        fp_cache[mid] = gen.GetFingerprint(mol)
        v = np.array(rdMolDescriptors.GetMorganFingerprintAsBitVect(
            mol, RADIUS, nBits=FP_SIZE), dtype=np.float32)
        r = ses.run(None, {nin: v})
        probs[mid] = np.array([x[0] for x in r], dtype=np.float32)
        if (k + 1) % 5000 == 0:
            print(f"    pontuadas {k+1:,}/{len(universo):,} "
                  f"({(time.time()-t)/60:.1f} min)", flush=True)
    for mid, sm in refs_todas.items():
        if mid in fp_cache:
            continue
        mol = Chem.MolFromSmiles(sm)
        if mol is not None:
            fp_cache[mid] = gen.GetFingerprint(mol)
    print(f"  pontuacao concluida ({(time.time()-t)/60:.1f} min)\n", flush=True)

    # ---- por alvo ---------------------------------------------------------
    linhas, ag = [], {"p": [], "r": []}
    for j, (alvo, v) in enumerate(sorted(por_alvo.items())):
        c = col[alvo]
        fps_ref = [fp_cache[m] for m in v["ref"] if m in fp_cache]
        if len(fps_ref) < 10:
            continue
        ids_t = [m for m in v["tes"] if m in probs]
        ids_d = [m for m in decoys
                 if m in probs and m not in v["medidos"]]
        if len(ids_t) < 10 or len(ids_d) < 100:
            continue
        ids = ids_t + ids_d
        p_all = np.array([probs[m][c] for m in ids])
        rot = np.array([1] * len(ids_t) + [0] * len(ids_d))
        sim = np.array([max(DataStructs.BulkTanimotoSimilarity(
            fp_cache[m], fps_ref)) for m in ids])

        for marca, sel in (("tudo", np.ones(len(rot), bool)),
                           ("novidade", sim < CORTE_NOVIDADE)):
            na, ni = int(rot[sel].sum()), int((1 - rot[sel]).sum())
            if na < MIN_ATIVOS_ESTRATO or ni < 50:
                continue
            res = mv.resumo(p_all[sel].tolist(), rot[sel].tolist())
            base = na / (na + ni)
            # Precisao@k DENTRO do alvo. O agregado por posto percentual
            # empata os topos dos 145 alvos em 1,0 e o desempate pessimista
            # joga todo decoy empatado para a frente, o que zera prec@k para
            # k < numero de alvos. Por alvo nao ha esse artefato, e e a
            # pergunta real: dos k que eu ensaiaria neste alvo, quantos sao
            # ativos de verdade?
            ordem = sorted(zip(p_all[sel].tolist(), rot[sel].tolist()),
                           key=lambda q: (-q[0], q[1]))
            acertos = {k: sum(r for _, r in ordem[:k]) for k in (5, 10, 20)}
            linhas.append(dict(
                alvo=alvo, estrato=marca, n=na + ni, ativos=na,
                ref=len(fps_ref), taxa_base=round(base, 5),
                auroc=round(res["auroc"], 4), ef_1=round(res["ef_1%"], 3),
                bedroc=round(res["bedroc_20.5"], 4),
                acertos_5=acertos[5], acertos_10=acertos[10],
                acertos_20=acertos[20]))
            if marca == "novidade":
                x = p_all[sel]
                o = np.argsort(x, kind="mergesort")
                r = np.empty(len(x)); r[o] = np.arange(len(x), dtype=float)
                ag["p"].extend((r / max(1, len(x) - 1)).tolist())
                ag["r"].extend(rot[sel].tolist())
        if (j + 1) % 20 == 0:
            print(f"  medidos {j+1}/{len(por_alvo)} alvos", flush=True)

    with open(SAIDA, "w", newline="") as fh:
        campos = ["alvo", "estrato", "n", "ativos", "ref", "taxa_base",
                  "auroc", "ef_1", "bedroc", "acertos_5", "acertos_10",
                  "acertos_20"]
        w = csv.DictWriter(fh, fieldnames=campos, delimiter="\t")
        w.writeheader(); w.writerows(linhas)

    import statistics as st
    for marca in ("tudo", "novidade"):
        sel = [l for l in linhas if l["estrato"] == marca]
        if not sel:
            continue
        print(f"\n=== {marca.upper()}: {len(sel)} alvos ===")
        print(f"  AUROC  mediana {st.median(l['auroc'] for l in sel):.3f}  "
              f"media {st.mean(l['auroc'] for l in sel):.3f}")
        print(f"  EF1%   mediana {st.median(l['ef_1'] for l in sel):.2f}  "
              f"media {st.mean(l['ef_1'] for l in sel):.2f}")
        print(f"  alvos com AUROC > 0,5: "
              f"{sum(1 for l in sel if l['auroc'] > 0.5)}/{len(sel)}")
        print(f"  ativos somados: {sum(l['ativos'] for l in sel):,}")
        for k in (5, 10, 20):
            ac = sum(l[f"acertos_{k}"] for l in sel)
            tot = k * len(sel)
            bm = st.mean(l["taxa_base"] for l in sel)
            print(f"  precisao@{k:<3} por alvo, somada = {ac}/{tot} = "
                  f"{ac/tot*100:5.2f}%  (taxa-base media {bm*100:.3f}%, "
                  f"ganho {ac/tot/bm:.0f}x)")
            q = sum(1 for l in sel if l[f"acertos_{k}"] > 0)
            print(f"      alvos com >= 1 acerto no top-{k}: {q}/{len(sel)}")

    if sum(ag["r"]) >= 5:
        n, na = len(ag["r"]), int(sum(ag["r"]))
        a = mv.auroc(ag["p"], ag["r"])
        q1, q2 = a / (2 - a), 2 * a * a / (1 + a)
        nn = n - na
        se = math.sqrt((a * (1 - a) + (na - 1) * (q1 - a**2)
                        + (nn - 1) * (q2 - a**2)) / (na * nn))
        pares = sorted(zip(ag["p"], ag["r"]), key=lambda x: (-x[0], x[1]))
        print("\n=== AGREGADO, estrato de novidade, todos os alvos ===")
        print(f"n = {n:,}  ativos = {na:,}  taxa-base = {na/n*100:.3f}%")
        print(f"AUROC = {a:.4f}  EP {se:.4f}  "
              f"IC95% [{a-1.96*se:.4f}, {a+1.96*se:.4f}]  z = {(a-0.5)/se:.2f}")
        print(f"p (uma cauda) = {math.erfc((a-0.5)/se/math.sqrt(2))/2:.3e}")
        for k in (20, 50, 100, 200):
            ac = sum(r for _, r in pares[:k])
            print(f"  precisao@{k:<4} = {ac}/{k} = {ac/k*100:5.2f}%  "
                  f"(limite pessimista: empates no topo)")
    print(f"\n{len(linhas)} linhas em {SAIDA} ({(time.time()-t0)/60:.1f} min)")


if __name__ == "__main__":
    main()
