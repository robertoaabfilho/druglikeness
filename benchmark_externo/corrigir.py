"""As tres correcoes que os Adendos 6 e 7 pediram, medidas no mesmo dado.

1. PROXY DE DOMINIO (Adendo 7, defeito 1). `p_max` esta saturado: mediana
   0,968 em antibiotico e 0,968 em bioativo qualquer. Testa substitutos pelo
   criterio certo -- um proxy de dominio nao deve prever atividade, deve
   prever CONFIANCA. Logo: dentro do estrato de novidade, parte as moleculas
   na mediana do proxy e compara a AUROC do modelo nas duas metades. Proxy
   que serve tem AUROC maior na metade alta.

2. VIES DE TAMANHO (Adendo 7, defeito 2). O ranque e explicado por peso
   molecular (rho +0,407 com p_max) e isso passa por baixo do filtro de
   novidade (rho +0,051). Mede a inclinacao de cada uma das 616 cabecas
   contra o numero de bits ligados do ECFP4 e usa o residuo como escore
   corrigido. Reavalia o holdout do Adendo 6 com o escore corrigido.

3. QUEM SAO OS 19 ALVOS QUE ACERTAM (Adendo 6). Os preditores obvios --
   AUROC no estrato analogo, EF1%, BEDROC -- sao RETROSPECTIVOS: medem-se no
   mesmo conjunto de teste pos-2022. Nao servem de regra de decisao. Aqui so
   entram preditores disponiveis ANTES de conhecer o desfecho: contagem e
   diversidade das referencias <= 2020, classe de proteina, e a forma da
   distribuicao de escore do modelo sobre a propria biblioteca.

Desvio declarado em relacao a holdout_escalado.py: nao refaz a chamada
`contar()` de guarda de tamanho (no Adendo 6 ela nao excluiu nenhum alvo), de
modo que esta rodada usa so o cache de atividades.
"""
import csv
import json
import math
import random
import sys
import time
import urllib.request
from pathlib import Path

import caminhos_dados as cd_

AQUI = Path(__file__).parent
sys.path.insert(0, str(AQUI))

# O modelo clonado e o cache de atividades sao dado baixado (centenas de MB),
# ficam fora do repositorio. Procura no diretorio do script e no rascunho da
# sessao; DADOS_BENCHMARK sobrescreve.
import os
def _base():
    cand = [Path(os.environ["DADOS_BENCHMARK"])] if os.environ.get(
        "DADOS_BENCHMARK") else []
    cand.append(AQUI)
    cand += sorted(Path("/tmp").glob(
        "claude-*/*/*/scratchpad"), reverse=True)
    for d in cand:
        if (d / "mt/trained_models/chembl_28_model").is_dir():
            return d
    raise SystemExit("nao achei mt/trained_models/chembl_28_model; "
                     "defina DADOS_BENCHMARK")
DADOS = _base()
MT28 = DADOS / "mt/trained_models/chembl_28_model"
LITPCBA = cd_.de_ambiente("LITPCBA_FULL")
SONDAGEM = AQUI / "sondagem_alvos.json"
CACHE_ATIV = DADOS / "cache_chembl_dor"
CLASSES = AQUI / "classes_alvo.json"
SAIDA_PROXY = AQUI / "resultado_proxy_dominio.tsv"
SAIDA_VIES = AQUI / "resultado_vies_tamanho.tsv"
SAIDA_ALVO = AQUI / "resultado_preditores_alvo.tsv"
AGENTE = "druglikeness-benchmark/1.0 (research script)"

ANO_REF_MAX, ANO_TESTE_MIN = 2020, 2022
CORTE_NOVIDADE = 0.25
FP_SIZE, RADIUS = 1024, 2
MIN_RECENTES, MAX_ALVOS, SEMENTE = 50, 150, 20261128
MIN_ATIVOS_ESTRATO = 3
AMOSTRA_DIVERSIDADE = 300

import numpy as np
import onnxruntime
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors, rdFingerprintGenerator, DataStructs
RDLogger.DisableLog("rdApp.*")

import metricas_vs as mv
from calibrar_dor import buscar_atividades


def postos(v):
    """Postos medios (trata empate), base 0."""
    o = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
            j += 1
        m = (i + j) / 2.0
        for k in range(i, j + 1):
            r[o[k]] = m
        i = j + 1
    return r


def spearman(x, y):
    rx, ry = postos(x), postos(y)
    n = len(x)
    if n < 3:
        return float("nan")
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx)
                    * sum((b - my) ** 2 for b in ry))
    return num / den if den else float("nan")


def auroc_ep(p, r):
    """AUROC com erro padrao de Hanley-McNeil."""
    n, na = len(r), int(sum(r))
    nn = n - na
    if na < 2 or nn < 2:
        return float("nan"), float("nan"), na
    a = mv.auroc(p, r)
    q1, q2 = a / (2 - a), 2 * a * a / (1 + a)
    ep = math.sqrt((a * (1 - a) + (na - 1) * (q1 - a ** 2)
                    + (nn - 1) * (q2 - a ** 2)) / (na * nn))
    return a, ep, na


def pct(v):
    """Posto percentual dentro do vetor, para empilhar alvos diferentes."""
    r = postos(list(v))
    return [x / max(1, len(v) - 1) for x in r]


def classe_do_alvo(alvo, cache):
    """Classe de proteina de nivel 1 do ChEMBL (enzyme, membrane receptor...)."""
    if alvo in cache:
        return cache[alvo]
    def pega(url):
        req = urllib.request.Request(url, headers={"User-Agent": AGENTE})
        for t in range(4):
            try:
                with urllib.request.urlopen(req, timeout=60) as r:
                    return json.loads(r.read())
            except Exception:
                if t == 3:
                    return None
                time.sleep(2 * (t + 1))
    base = "https://www.ebi.ac.uk/chembl/api/data"
    d = pega(f"{base}/target/{alvo}.json")
    r = {"tipo": None, "classe": None}
    if d:
        r["tipo"] = d.get("target_type")
        comps = d.get("target_components") or []
        if comps:
            c = pega(f"{base}/target_component/{comps[0]['component_id']}.json")
            pcs = (c or {}).get("protein_classifications") or []
            if pcs:
                pc = pega(f"{base}/protein_classification/"
                          f"{pcs[0]['protein_classification_id']}.json")
                desc = (pc or {}).get("protein_class_desc") or ""
                if desc:
                    r["classe"] = desc.split()[0]
    cache[alvo] = r
    CLASSES.write_text(json.dumps(cache))
    return r


def main():
    t0 = time.time()
    painel = json.load(open(MT28 / "targets_28.json"))
    col = {c: i for i, c in enumerate(painel)}
    sond = json.loads(SONDAGEM.read_text())
    elegiveis = sorted(a for a, v in sond.items()
                       if a in col and v.get("recentes")
                       and v["recentes"] >= MIN_RECENTES)
    rnd = random.Random(SEMENTE)
    escolhidos = (sorted(rnd.sample(elegiveis, MAX_ALVOS))
                  if len(elegiveis) > MAX_ALVOS else elegiveis)
    print(f"{len(elegiveis)} elegiveis -> {len(escolhidos)} "
          f"(semente {SEMENTE}, igual ao Adendo 6)", flush=True)

    opc = onnxruntime.SessionOptions()
    opc.intra_op_num_threads = 4
    ses = onnxruntime.InferenceSession(str(MT28 / "chembl_28_multitask.onnx"),
                                       sess_options=opc,
                                       providers=["CPUExecutionProvider"])
    nin = ses.get_inputs()[0].name
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

    # ---- coleta (so cache) ------------------------------------------------
    por_alvo = {}
    for alvo in escolhidos:
        arq = CACHE_ATIV / f"{alvo}.json"
        if not arq.is_file():
            continue
        pm = {}
        for a in json.loads(arq.read_text()):
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
    print(f"{len(por_alvo)} alvos do cache\n", flush=True)

    # ---- pontuar o universo uma vez, guardando nbits ----------------------
    universo = dict(decoys)
    for v in por_alvo.values():
        universo.update(v["tes"])
    refs_todas = {}
    for v in por_alvo.values():
        refs_todas.update(v["ref"])
    print(f"universo {len(universo):,} | referencias {len(refs_todas):,}",
          flush=True)

    fp_cache, probs, nbits = {}, {}, {}
    t = time.time()
    for k, (mid, sm) in enumerate(universo.items()):
        mol = Chem.MolFromSmiles(sm)
        if mol is None:
            continue
        fp_cache[mid] = gen.GetFingerprint(mol)
        v = np.array(rdMolDescriptors.GetMorganFingerprintAsBitVect(
            mol, RADIUS, nBits=FP_SIZE), dtype=np.float32)
        nbits[mid] = float(v.sum())
        probs[mid] = np.array([x[0] for x in ses.run(None, {nin: v})],
                              dtype=np.float32)
        if (k + 1) % 10000 == 0:
            print(f"  {k+1:,}/{len(universo):,} "
                  f"({(time.time()-t)/60:.1f} min)", flush=True)
    for mid, sm in refs_todas.items():
        if mid not in fp_cache:
            mol = Chem.MolFromSmiles(sm)
            if mol is not None:
                fp_cache[mid] = gen.GetFingerprint(mol)
    ids = sorted(probs)
    P = np.stack([probs[m] for m in ids])          # (mol, 616)
    B = np.array([nbits[m] for m in ids])
    print(f"matriz {P.shape} pontuada ({(time.time()-t0)/60:.1f} min)\n",
          flush=True)

    # =====================================================================
    # CORRECAO 2: inclinacao de cada cabeca contra o numero de bits
    # =====================================================================
    eps = 1e-6
    Lg = np.log(np.clip(P, eps, 1 - eps) / (1 - np.clip(P, eps, 1 - eps)))
    bz = (B - B.mean()) / B.std()
    coef = (Lg * bz[:, None]).mean(axis=0)         # inclinacao por cabeca
    Res = Lg - coef[None, :] * bz[:, None]         # escore corrigido
    rho_cab = np.array([spearman(B.tolist(), P[:, j].tolist())
                        for j in range(P.shape[1])])
    var_expl = (coef ** 2) / np.maximum(Lg.var(axis=0), 1e-12)
    print("=== CORRECAO 2: vies de tamanho por cabeca (616 cabecas) ===")
    for q in (10, 25, 50, 75, 90):
        print(f"  rho(p, n_bits) percentil {q:>2}: "
              f"{np.percentile(rho_cab, q):+.3f}")
    print(f"  cabecas com rho > +0,3: "
          f"{int((rho_cab > 0.3).sum())}/{len(rho_cab)}")
    print(f"  cabecas com rho > +0,5: "
          f"{int((rho_cab > 0.5).sum())}/{len(rho_cab)}")
    print(f"  fracao da variancia do logit explicada por tamanho: "
          f"mediana {np.median(var_expl)*100:.1f}%  "
          f"p90 {np.percentile(var_expl, 90)*100:.1f}%", flush=True)
    idx = {m: i for i, m in enumerate(ids)}
    with open(SAIDA_VIES, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["alvo_idx", "alvo", "rho_p_nbits", "coef_logit",
                    "var_explicada"])
        for j, a in enumerate(painel):
            w.writerow([j, a, round(float(rho_cab[j]), 4),
                        round(float(coef[j]), 4),
                        round(float(var_expl[j]), 5)])

    # =====================================================================
    # por alvo: estrato de novidade, escores cru e corrigido, proxies
    # =====================================================================
    PROXIES = ["p_max", "p_max_corrigido", "sim_max", "sim_top5",
               "n_cabecas_90", "n_bits"]
    emp = {k: {"alta": {"p": [], "r": []}, "baixa": {"p": [], "r": []}}
           for k in PROXIES}
    agr = {"cru": {"p": [], "r": []}, "corrigido": {"p": [], "r": []}}
    dec = {q: {"p": [], "r": [], "v": []} for q in range(10)}
    linhas_alvo = []
    classes = json.loads(CLASSES.read_text()) if CLASSES.is_file() else {}

    for j, (alvo, v) in enumerate(sorted(por_alvo.items())):
        c = col[alvo]
        fps_ref = [fp_cache[m] for m in v["ref"] if m in fp_cache]
        if len(fps_ref) < 10:
            continue
        it = [m for m in v["tes"] if m in probs]
        idd = [m for m in decoys if m in probs and m not in v["medidos"]]
        if len(it) < 10 or len(idd) < 100:
            continue
        mids = it + idd
        li = [idx[m] for m in mids]
        rot = np.array([1] * len(it) + [0] * len(idd))
        sims = [sorted(DataStructs.BulkTanimotoSimilarity(
            fp_cache[m], fps_ref), reverse=True) for m in mids]
        sim_max = np.array([s[0] for s in sims])
        sim5 = np.array([float(np.mean(s[:5])) for s in sims])
        sel = sim_max < CORTE_NOVIDADE
        na = int(rot[sel].sum())
        if na < MIN_ATIVOS_ESTRATO or int((1 - rot[sel]).sum()) < 50:
            continue

        cru = P[li, c]
        corr = Res[li, c]
        rs = rot[sel].tolist()
        agr["cru"]["p"].extend(pct(cru[sel])); agr["cru"]["r"].extend(rs)
        agr["corrigido"]["p"].extend(pct(corr[sel]))
        agr["corrigido"]["r"].extend(rs)

        val = {
            "p_max": P[li].max(axis=1),
            "p_max_corrigido": Res[li].max(axis=1),
            "sim_max": sim_max,
            "sim_top5": sim5,
            "n_cabecas_90": (P[li] > 0.9).sum(axis=1).astype(float),
            "n_bits": B[li],
        }
        # Decis de sim_top5 dentro do estrato: onde o sinal morre?
        xs5 = sim5[sel]
        for q in range(10):
            lo, hi = np.percentile(xs5, q * 10), np.percentile(xs5,
                                                               (q + 1) * 10)
            m2 = (xs5 >= lo) & (xs5 <= hi) if q == 9 else (
                (xs5 >= lo) & (xs5 < hi))
            rr = rot[sel][m2]
            if rr.sum() < 1 or len(rr) < 20:
                continue
            dec[q]["p"].extend(pct(cru[sel][m2]))
            dec[q]["r"].extend(rr.tolist())
            dec[q]["v"].extend(xs5[m2].tolist())
        for nome, x in val.items():
            xs = x[sel]
            med = float(np.median(xs))
            for marca, m2 in (("alta", xs >= med), ("baixa", xs < med)):
                rr = rot[sel][m2]
                if rr.sum() < 1 or len(rr) < 20:
                    continue
                emp[nome][marca]["p"].extend(pct(cru[sel][m2]))
                emp[nome][marca]["r"].extend(rr.tolist())

        # ---- preditores PROSPECTIVOS do alvo
        amostra = (random.Random(SEMENTE).sample(fps_ref, AMOSTRA_DIVERSIDADE)
                   if len(fps_ref) > AMOSTRA_DIVERSIDADE else fps_ref)
        soma, cont = 0.0, 0
        for i2 in range(len(amostra) - 1):
            s = DataStructs.BulkTanimotoSimilarity(amostra[i2],
                                                   amostra[i2 + 1:])
            soma += float(sum(s)); cont += len(s)
        div = 1.0 - (soma / cont if cont else 0.0)
        cl = classe_do_alvo(alvo, classes)
        colx = cru[sel]
        a_cru, _e, _n = auroc_ep(pct(colx).copy(), rs)
        a_cor, _e2, _n2 = auroc_ep(pct(corr[sel]).copy(), rs)
        ordem = sorted(zip(colx.tolist(), rs), key=lambda q: (-q[0], q[1]))
        linhas_alvo.append(dict(
            alvo=alvo, classe=cl.get("classe"), tipo=cl.get("tipo"),
            n_ref=len(fps_ref), diversidade_ref=round(div, 4),
            n_novidade=int(sel.sum()), ativos_novos=na,
            auroc_cru=round(a_cru, 4), auroc_corrigido=round(a_cor, 4),
            acertos_20=sum(r for _, r in ordem[:20]),
            p50_coluna=round(float(np.percentile(colx, 50)), 4),
            p90_coluna=round(float(np.percentile(colx, 90)), 4),
            razao_p90_p50=round(float(np.percentile(colx, 90)
                                      / max(np.percentile(colx, 50), 1e-6)), 3),
            max_coluna=round(float(colx.max()), 4)))
        if (j + 1) % 20 == 0:
            print(f"  medidos {j+1}/{len(por_alvo)} "
                  f"({(time.time()-t0)/60:.1f} min)", flush=True)

    with open(SAIDA_ALVO, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(linhas_alvo[0].keys()),
                           delimiter="\t")
        w.writeheader(); w.writerows(linhas_alvo)

    # =====================================================================
    # resultado da CORRECAO 2 no holdout
    # =====================================================================
    print("\n=== CORRECAO 2: o escore corrigido melhora o holdout? ===")
    for nome in ("cru", "corrigido"):
        a, ep, na = auroc_ep(agr[nome]["p"], agr[nome]["r"])
        print(f"  {nome:<10} AUROC {a:.4f}  EP {ep:.4f}  "
              f"IC95% [{a-1.96*ep:.4f}, {a+1.96*ep:.4f}]  (n_ativos {na})")
    ac = sum(l["acertos_20"] for l in linhas_alvo)
    print(f"  precisao@20 por alvo, escore cru: {ac}/{20*len(linhas_alvo)} = "
          f"{ac/(20*len(linhas_alvo))*100:.2f}%")
    import statistics as st
    print(f"  AUROC mediana por alvo: cru "
          f"{st.median(l['auroc_cru'] for l in linhas_alvo):.4f}  "
          f"corrigido "
          f"{st.median(l['auroc_corrigido'] for l in linhas_alvo):.4f}")
    melhorou = sum(1 for l in linhas_alvo
                   if l["auroc_corrigido"] > l["auroc_cru"])
    print(f"  alvos em que corrigir melhorou: {melhorou}/{len(linhas_alvo)}")

    # =====================================================================
    # resultado da CORRECAO 1
    # =====================================================================
    print("\n=== CORRECAO 1: qual proxy de dominio preve confianca? ===")
    print(f"  {'proxy':<18} {'AUROC alta':>11} {'AUROC baixa':>12} "
          f"{'delta':>8}  {'n_ativos a/b':>13}")
    linhas_proxy = []
    for nome in PROXIES:
        aa, ea, nna = auroc_ep(emp[nome]["alta"]["p"], emp[nome]["alta"]["r"])
        ab, eb, nnb = auroc_ep(emp[nome]["baixa"]["p"],
                               emp[nome]["baixa"]["r"])
        d = aa - ab
        epd = math.sqrt(ea ** 2 + eb ** 2)
        print(f"  {nome:<18} {aa:>11.4f} {ab:>12.4f} {d:>+8.4f}  "
              f"{nna:>6}/{nnb:<6}  (EP do delta {epd:.4f}, "
              f"z {d/epd:+.2f})")
        linhas_proxy.append(dict(proxy=nome, auroc_alta=round(aa, 4),
                                 ep_alta=round(ea, 4),
                                 auroc_baixa=round(ab, 4),
                                 ep_baixa=round(eb, 4), delta=round(d, 4),
                                 ep_delta=round(epd, 4),
                                 z=round(d / epd, 2),
                                 ativos_alta=nna, ativos_baixa=nnb))
    with open(SAIDA_PROXY, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(linhas_proxy[0].keys()),
                           delimiter="\t")
        w.writeheader(); w.writerows(linhas_proxy)

    print("\n=== CORRECAO 1b: onde o sinal morre, por decil de sim_top5 ===")
    print(f"  {'decil':<7} {'sim_top5 medio':>15} {'AUROC':>8} {'EP':>7} "
          f"{'ativos':>7}")
    linhas_dec = []
    for q in range(10):
        if len(dec[q]["r"]) < 50:
            continue
        a, ep, na = auroc_ep(dec[q]["p"], dec[q]["r"])
        mv5 = sum(dec[q]["v"]) / len(dec[q]["v"])
        print(f"  {q+1:<7} {mv5:>15.4f} {a:>8.4f} {ep:>7.4f} {na:>7}")
        linhas_dec.append(dict(decil=q + 1, sim_top5_medio=round(mv5, 4),
                               auroc=round(a, 4), ep=round(ep, 4),
                               ativos=na, n=len(dec[q]["r"])))
    with open(AQUI / "resultado_decis_similaridade.tsv", "w",
              newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(linhas_dec[0].keys()),
                           delimiter="\t")
        w.writeheader(); w.writerows(linhas_dec)

    # =====================================================================
    # resultado da CORRECAO 3
    # =====================================================================
    print("\n=== CORRECAO 3: preditores PROSPECTIVOS do sucesso no alvo ===")
    acerta = [1.0 if l["acertos_20"] > 0 else 0.0 for l in linhas_alvo]
    aur = [l["auroc_cru"] for l in linhas_alvo]
    print(f"  {len(linhas_alvo)} alvos, {int(sum(acerta))} acertam no top-20")
    print(f"  {'preditor':<24} {'rho vs acerta':>14} {'rho vs AUROC':>13}")
    for k in ("n_ref", "diversidade_ref", "n_novidade", "p50_coluna",
              "p90_coluna", "razao_p90_p50", "max_coluna"):
        x = [float(l[k]) for l in linhas_alvo]
        print(f"  {k:<24} {spearman(x, acerta):>+14.3f} "
              f"{spearman(x, aur):>+13.3f}")
    print("\n  por classe de proteina:")
    pc = {}
    for l, a in zip(linhas_alvo, acerta):
        pc.setdefault(l["classe"] or "sem classe", []).append((a, l["auroc_cru"]))
    for k, v in sorted(pc.items(), key=lambda x: -len(x[1])):
        print(f"    {k:<20} {len(v):>3} alvos  acertam "
              f"{int(sum(a for a, _ in v)):>2} ({sum(a for a,_ in v)/len(v)*100:>3.0f}%)"
              f"  AUROC mediana {st.median(b for _, b in v):.3f}")
    print(f"\n{(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
