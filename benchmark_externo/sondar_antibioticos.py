"""O que o pipeline faz com uma lista de antibioticos.

Nao mede acerto -- nao existe ativo conhecido para "mimetismo molecular na
GBS". Mede o que DA para medir antes de confiar na saida:

1. Dominio de aplicabilidade: Tanimoto ECFP4 maximo de cada antibiotico aos
   193 mil ativos de referencia dos 145 alvos do Adendo 6. Diz em qual
   estrato do Adendo 6 a lista cai -- e portanto qual precisao esperar.
2. Familiaridade: probabilidade MAXIMA entre as 616 cabecas. O proxy de
   dominio que o projeto usa.
3. Qual propriedade fisico-quimica explica o ranque. Se o ranque for
   previsivel por logP/peso, o modelo nao esta acrescentando nada.
4. Ranque nos alvos imunes do painel que teriam alguma plausibilidade em
   GBS, comparado com os ativos do LIT-PCBA como referencia de "composto
   bioativo qualquer".

Lista: classe ATC J01 (antibacterianos sistemicos) do ChEMBL.
"""
import csv, json, sys, time, urllib.parse, urllib.request
from pathlib import Path

import caminhos_dados as cd_

AQUI = Path(__file__).parent
sys.path.insert(0, str(AQUI))
MT28 = AQUI / "mt/trained_models/chembl_28_model"
CACHE = AQUI / "cache_chembl_dor"
LITPCBA = cd_.de_ambiente("LITPCBA_FULL")
SAIDA = AQUI / "resultado_antibioticos.tsv"
LISTA = AQUI / "lista_j01.json"
AGENTE = "druglikeness-benchmark/1.0 (research script)"
FP_SIZE, RADIUS = 1024, 2
ANO_REF_MAX = 2020
PCHEMBL_MIN = 6.0
CORTE_NOVIDADE = 0.25

IMUNES = {"CHEMBL2413": "C5aR1 (complemento)", "CHEMBL4296": "S1PR1",
          "CHEMBL321": "MMP-9", "CHEMBL5251": "BTK", "CHEMBL2835": "JAK1",
          "CHEMBL2971": "JAK2", "CHEMBL1829": "TNF-alfa"}

import numpy as np
import onnxruntime
from rdkit import Chem, RDLogger
from rdkit.Chem import (rdMolDescriptors, rdFingerprintGenerator,
                        DataStructs, Descriptors, Crippen)
RDLogger.DisableLog("rdApp.*")


def baixar_j01():
    if LISTA.is_file():
        return json.loads(LISTA.read_text())
    tudo, deslo = [], 0
    while True:
        q = urllib.parse.urlencode({
            "atc_classifications__level2": "J01",
            "only": "molecule_chembl_id,pref_name,molecule_structures",
            "limit": 200, "offset": deslo})
        url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?{q}"
        req = urllib.request.Request(url, headers={"User-Agent": AGENTE})
        for t in range(5):
            try:
                with urllib.request.urlopen(req, timeout=120) as r:
                    d = json.loads(r.read())
                break
            except Exception:
                if t == 4:
                    raise
                time.sleep(3 * (t + 1))
        lote = d.get("molecules", [])
        tudo.extend(lote)
        deslo += len(lote)
        print(f"  J01: {deslo}/{d['page_meta']['total_count']}", flush=True)
        if deslo >= d["page_meta"]["total_count"] or not lote:
            break
        time.sleep(0.3)
    LISTA.write_text(json.dumps(tudo))
    return tudo


def main():
    t0 = time.time()
    painel = json.load(open(MT28 / "targets_28.json"))
    col = {c: i for i, c in enumerate(painel)}
    ses = onnxruntime.InferenceSession(
        str(MT28 / "chembl_28_multitask.onnx"),
        providers=["CPUExecutionProvider"])
    nin = ses.get_inputs()[0].name
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=RADIUS,
                                                    fpSize=FP_SIZE)

    def fp_e_prob(sm):
        mol = Chem.MolFromSmiles(sm)
        if mol is None:
            return None
        v = np.array(rdMolDescriptors.GetMorganFingerprintAsBitVect(
            mol, RADIUS, nBits=FP_SIZE), dtype=np.float32)
        p = np.array([x[0] for x in ses.run(None, {nin: v})],
                     dtype=np.float32)
        return mol, gen.GetFingerprint(mol), p

    # --- lista J01
    mols = baixar_j01()
    abx = {}
    for m in mols:
        # O cache gravado neste repositorio e a forma enxuta
        # (id/nome/smiles); a resposta crua da API traz o SMILES dentro de
        # molecule_structures. Aceita as duas.
        s = (m.get("canonical_smiles")
             or (m.get("molecule_structures") or {}).get("canonical_smiles"))
        if s:
            abx[m["molecule_chembl_id"]] = (m.get("pref_name") or
                                            m["molecule_chembl_id"], s)
    print(f"\n{len(abx)} antibioticos J01 com estrutura", flush=True)

    # --- referencias: ativos <= 2020 dos alvos coletados no Adendo 6
    refs = {}
    for arq in sorted(CACHE.glob("CHEMBL*.json")):
        for a in json.loads(arq.read_text()):
            if a.get("assay_type") not in ("B", "F"):
                continue
            sm, mid, ano = (a.get("canonical_smiles"),
                            a.get("molecule_chembl_id"),
                            a.get("document_year"))
            if sm and mid and ano and int(ano) <= ANO_REF_MAX:
                refs[mid] = sm
    print(f"{len(refs):,} ativos de referencia (<= {ANO_REF_MAX})", flush=True)
    fps_ref = []
    for sm in refs.values():
        mol = Chem.MolFromSmiles(sm)
        if mol is not None:
            fps_ref.append(gen.GetFingerprint(mol))
    print(f"{len(fps_ref):,} fingerprints de referencia "
          f"({(time.time()-t0)/60:.1f} min)", flush=True)

    # --- decoys LIT-PCBA como comparacao de "bioativo qualquer"
    dec = {}
    for d in sorted(LITPCBA.iterdir()):
        f = d / "actives.smi"
        if f.is_file():
            for linha in open(f):
                p = linha.split()
                if len(p) >= 2:
                    dec[p[1]] = p[0]
    print(f"{len(dec):,} ativos LIT-PCBA (comparacao)\n", flush=True)

    def medir(conj, marca):
        linhas = []
        for i, (mid, par) in enumerate(conj.items()):
            nome, sm = par if isinstance(par, tuple) else (mid, par)
            r = fp_e_prob(sm)
            if r is None:
                continue
            mol, fp, p = r
            sim = max(DataStructs.BulkTanimotoSimilarity(fp, fps_ref))
            linhas.append(dict(
                grupo=marca, id=mid, nome=nome,
                pm=round(Descriptors.MolWt(mol), 1),
                clogp=round(Crippen.MolLogP(mol), 2),
                tpsa=round(rdMolDescriptors.CalcTPSA(mol), 1),
                sim_ref=round(sim, 4),
                p_max=round(float(p.max()), 4),
                alvo_max=painel[int(p.argmax())],
                **{f"p_{v}": round(float(p[col[k]]), 4)
                   for k, v in IMUNES.items()}))
            if (i + 1) % 100 == 0:
                print(f"  {marca} {i+1}/{len(conj)} "
                      f"({(time.time()-t0)/60:.1f} min)", flush=True)
        return linhas

    L = medir(abx, "antibiotico")
    import random
    rnd = random.Random(20261128)
    amostra = dict(rnd.sample(sorted(dec.items()), 600))
    L += medir(amostra, "bioativo_litpcba")

    campos = list(L[0].keys())
    with open(SAIDA, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=campos, delimiter="\t")
        w.writeheader(); w.writerows(L)

    import statistics as st
    for g in ("antibiotico", "bioativo_litpcba"):
        S = [l for l in L if l["grupo"] == g]
        print(f"\n=== {g}: {len(S)} compostos ===")
        for c in ("sim_ref", "p_max", "pm", "clogp", "tpsa"):
            v = [l[c] for l in S]
            print(f"  {c:<8} mediana {st.median(v):8.3f}  "
                  f"q1 {st.quantiles(v)[0]:8.3f}  q3 {st.quantiles(v)[2]:8.3f}")
        nov = sum(1 for l in S if l["sim_ref"] < CORTE_NOVIDADE)
        print(f"  no estrato de novidade (sim < {CORTE_NOVIDADE}): "
              f"{nov}/{len(S)} = {nov/len(S)*100:.0f}%")

    S = [l for l in L if l["grupo"] == "antibiotico"]
    print("\n=== top 15 antibioticos por familiaridade (p_max) ===")
    for l in sorted(S, key=lambda x: -x["p_max"])[:15]:
        print(f"  {l['p_max']:.3f}  sim={l['sim_ref']:.2f}  "
              f"logP={l['clogp']:+5.2f}  {l['nome'][:34]:<34} "
              f"-> {l['alvo_max']}")
    for k, v in IMUNES.items():
        print(f"\n=== top 8 em {v} ===")
        for l in sorted(S, key=lambda x: -x[f"p_{v}"])[:8]:
            print(f"  {l[f'p_{v}']:.4f}  sim={l['sim_ref']:.2f}  "
                  f"{l['nome'][:40]}")
        d = [l[f"p_{v}"] for l in L if l["grupo"] == "bioativo_litpcba"]
        a = [l[f"p_{v}"] for l in S]
        print(f"  mediana antibiotico {st.median(a):.4f} | "
              f"bioativo qualquer {st.median(d):.4f}")

    # correlacao do ranque com fisico-quimica
    import math
    def spearman(x, y):
        def postos(v):
            o = sorted(range(len(v)), key=lambda i: v[i])
            r = [0.0] * len(v)
            for p, i in enumerate(o):
                r[i] = p
            return r
        rx, ry = postos(x), postos(y)
        n = len(x)
        mx, my = sum(rx) / n, sum(ry) / n
        num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
        den = math.sqrt(sum((a - mx) ** 2 for a in rx)
                        * sum((b - my) ** 2 for b in ry))
        return num / den if den else float("nan")
    print("\n=== o ranque e explicado por fisico-quimica? "
          "(Spearman, antibioticos) ===")
    for c in ("clogp", "pm", "tpsa", "sim_ref"):
        print(f"  p_max vs {c:<8} rho = "
              f"{spearman([l['p_max'] for l in S], [l[c] for l in S]):+.3f}")
    print(f"\n{len(L)} linhas em {SAIDA} ({(time.time()-t0)/60:.1f} min)")


if __name__ == "__main__":
    main()
