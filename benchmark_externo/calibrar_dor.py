"""Calibracao do canal predito NOS ALVOS DE DOR NEUROPATICA do projeto.

Troca o benchmark emprestado (LIT-PCBA) pelos alvos reais. Responde: qual a
precisao esperada do canal predito, no estrato de novidade, para ESTES alvos?

Desenho
-------
Ativos: ChEMBL, `pchembl >= 6`, ensaio de ligacao (B) ou funcional (F) --
  triagem HTS fica fora, porque o proprio projeto mediu que ela e ruido
  (recall@10 0,111 contra 0,272).

Divisao TEMPORAL por ano do primeiro documento que mediu o par
  composto/alvo: referencia = ate ANO_CORTE, teste = depois. Imita
  reposicionamento -- o que o modelo teria achado antes de se saber.

Decoys: os 10.033 ATIVOS do LIT-PCBA, que sao compostos bioativos e drug-like
  ativos em 15 alvos sem relacao com dor. E o controle `actives-as-decoys` de
  Chen et al. 2019, o mais severo: um modelo que so reconhece "parece
  farmaco" cai para AUROC 0,5 contra ele. Qualquer molecula que seja ativa no
  proprio alvo de dor e removida dos decoys.

Estrato de novidade: maior Tanimoto ECFP4 aos ativos de REFERENCIA < CORTE.
"""
import csv
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

REPO = Path("/home/user/druglikeness")
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).parent))

CACHE = Path(__file__).parent / "cache_chembl_dor"
LITPCBA = Path("/home/user/litpcba_download/full/LIT-PCBA_full")
CKPT = "/home/user/pp4ckd_download/extraido/models/ecfp4_fold_4.pt"
ROTULOS = "/home/user/jixing475/pp4ckd/data/target_labels.tsv"
SAIDA = Path(__file__).parent / "resultado_calibracao_dor.tsv"

ANO_CORTE = 2015
CORTE_NOVIDADE = 0.25
PCHEMBL_MIN = 6.0
LOTE = 2048
# Cliente generico: nunca identificar o usuario a um servico externo.
AGENTE = "druglikeness-benchmark/1.0 (research script)"

ALVOS = {
    "CHEMBL1919": "Canal de calcio alpha-2/delta-1",
    "CHEMBL3896": "Canal de calcio alpha-2/delta-2",
    "CHEMBL228":  "Transportador de serotonina (SERT)",
    "CHEMBL222":  "Transportador de noradrenalina (NET)",
    "CHEMBL1845": "Canal de sodio tipo 1 alpha",
    "CHEMBL4794": "TRPV1",
    "CHEMBL233":  "Receptor opioide mu",
    "CHEMBL2015": "NMDA 1",
    "CHEMBL1904": "NMDA 2B",
}

import numpy as np
import torch
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator, DataStructs
RDLogger.DisableLog("rdApp.*")

from mods import prever_alvos as pa
import metricas_vs as mv


def buscar_atividades(alvo):
    """Todas as atividades do alvo com pchembl >= PCHEMBL_MIN. Em cache."""
    CACHE.mkdir(exist_ok=True)
    arq = CACHE / f"{alvo}.json"
    if arq.is_file():
        return json.loads(arq.read_text())
    campos = ("molecule_chembl_id,canonical_smiles,document_year,"
              "assay_type,pchembl_value")
    tudo, deslo = [], 0
    while True:
        q = urllib.parse.urlencode({
            "target_chembl_id": alvo,
            "pchembl_value__gte": PCHEMBL_MIN,
            "only": campos, "limit": 1000, "offset": deslo})
        url = f"https://www.ebi.ac.uk/chembl/api/data/activity.json?{q}"
        req = urllib.request.Request(url, headers={"User-Agent": AGENTE})
        # A API cai de vez em quando (ConnectionResetError); tenta de novo com
        # espera crescente em vez de perder a coleta inteira.
        for tentativa in range(5):
            try:
                with urllib.request.urlopen(req, timeout=120) as r:
                    d = json.loads(r.read())
                break
            except Exception as erro:
                if tentativa == 4:
                    raise
                print(f"    rede falhou ({type(erro).__name__}), "
                      f"tentativa {tentativa + 2}/5", flush=True)
                time.sleep(3 * (tentativa + 1))
        lote = d.get("activities", [])
        tudo.extend(lote)
        total = d["page_meta"]["total_count"]
        deslo += len(lote)
        print(f"    {alvo}: {deslo}/{total}", flush=True)
        if deslo >= total or not lote:
            break
        time.sleep(0.4)
    arq.write_text(json.dumps(tudo))
    return tudo


def main():
    t0 = time.time()
    rede, cfg, _d, ck = pa.carregar_modelo(None, fp="ECFP4", checkpoint=CKPT,
                                           dispositivo="cpu")
    rm = pa.carregar_rotulos(ROTULOS)
    col = {r["chembl_id"]: i for i, r in enumerate(rm) if r["chembl_id"]}
    treino_mod = {r["chembl_id"]: r["n_ativos"] for r in rm if r["chembl_id"]}
    torch.set_num_threads(4)
    gen = rdFingerprintGenerator.GetMorganGenerator(
        radius=cfg["raio"], fpSize=cfg["bits"])
    print(f"checkpoint {ck}\n", flush=True)

    # --- decoys: ativos do LIT-PCBA (bioativos, drug-like, alvo sem relacao)
    decoys = {}
    for d in sorted(LITPCBA.iterdir()):
        f = d / "actives.smi"
        if not f.is_file():
            continue
        for linha in open(f):
            p = linha.split()
            if len(p) >= 2:
                decoys[p[1]] = p[0]          # id -> smiles
    print(f"decoys (ativos do LIT-PCBA): {len(decoys)} compostos\n", flush=True)

    def p_de(smiles_lista, coluna):
        """Probabilidade na cabeca do alvo; devolve (probs, fps) alinhados,
        pulando SMILES ilegivel."""
        probs, fps, bv = [], [], []

        def flush():
            if not bv:
                return
            x = torch.tensor(np.array(bv, dtype=np.float32))
            with torch.no_grad():
                probs.extend(torch.sigmoid(rede(x)[:, coluna]).numpy().tolist())
            bv.clear()

        for s in smiles_lista:
            mol = Chem.MolFromSmiles(s)
            if mol is None:
                continue
            fps.append(gen.GetFingerprint(mol))
            bv.append(gen.GetFingerprintAsNumPy(mol))
            if len(bv) >= LOTE:
                flush()
        flush()
        return np.array(probs), fps

    linhas = []
    ag = {"p": [], "r": []}
    for alvo, nome in ALVOS.items():
        print(f"=== {alvo}  {nome}", flush=True)
        ativ = buscar_atividades(alvo)
        # primeiro ano por molecula, so ensaio B/F
        por_mol = {}
        for a in ativ:
            if a.get("assay_type") not in ("B", "F"):
                continue
            sm = a.get("canonical_smiles")
            mid = a.get("molecule_chembl_id")
            ano = a.get("document_year")
            if not sm or not mid or not ano:
                continue
            ano = int(ano)
            if mid not in por_mol or ano < por_mol[mid][1]:
                por_mol[mid] = (sm, ano)
        ref = {m: v[0] for m, v in por_mol.items() if v[1] <= ANO_CORTE}
        tes = {m: v[0] for m, v in por_mol.items() if v[1] > ANO_CORTE}
        dec = {m: s for m, s in decoys.items() if m not in por_mol}
        print(f"    {len(por_mol)} ativos B/F pchembl>={PCHEMBL_MIN}; "
              f"referencia(<={ANO_CORTE})={len(ref)}, teste(>{ANO_CORTE})="
              f"{len(tes)}, decoys={len(dec)}", flush=True)
        if len(tes) < 10 or len(ref) < 10:
            print("    n insuficiente, pulado\n", flush=True)
            continue

        coluna = col[alvo]
        _pr, fps_ref = p_de(list(ref.values()), coluna)
        p_tes, fps_tes = p_de(list(tes.values()), coluna)
        p_dec, fps_dec = p_de(list(dec.values()), coluna)

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
                           key=lambda x: (-x[0], x[1]))[:20]
            p20 = sum(r for _, r in pares) / 20
            base = na / (na + ni)
            linhas.append(dict(
                alvo=alvo, nome=nome, estrato=marca,
                treino_modelo=treino_mod[alvo], n=na + ni, ativos=na,
                taxa_base=round(base, 5), auroc=round(res["auroc"], 4),
                ef_1=round(res["ef_1%"], 3), prec_20=round(p20, 4),
                ganho_20=round(p20 / base, 1) if base else ""))
            print(f"    [{marca}] n={na+ni} ativos={na} base={base*100:.2f}% "
                  f"AUROC={res['auroc']:.3f} EF1%={res['ef_1%']:.2f} "
                  f"prec@20={p20*100:.0f}% ({p20/base:.0f}x)", flush=True)
            if marca == "novidade":
                # agregado por posto percentual dentro do alvo
                v = p_all[sel]
                o = np.argsort(v, kind="mergesort")
                r = np.empty(len(v)); r[o] = np.arange(len(v), dtype=float)
                ag["p"].extend((r / max(1, len(v) - 1)).tolist())
                ag["r"].extend(rot[sel].tolist())
        print(flush=True)

        with open(SAIDA, "w", newline="") as fh:
            campos = ["alvo", "nome", "estrato", "treino_modelo", "n",
                      "ativos", "taxa_base", "auroc", "ef_1", "prec_20",
                      "ganho_20"]
            w = csv.DictWriter(fh, fieldnames=campos, delimiter="\t")
            w.writeheader(); w.writerows(linhas)

    if sum(ag["r"]) >= 5:
        import math
        n, na = len(ag["r"]), int(sum(ag["r"]))
        a = mv.auroc(ag["p"], ag["r"])
        q1, q2 = a / (2 - a), 2 * a * a / (1 + a)
        nn = n - na
        se = math.sqrt((a * (1 - a) + (na - 1) * (q1 - a**2)
                        + (nn - 1) * (q2 - a**2)) / (na * nn))
        pares = sorted(zip(ag["p"], ag["r"]), key=lambda x: (-x[0], x[1]))
        print("=== AGREGADO, estrato de novidade, alvos de dor ===")
        print(f"n = {n:,}  ativos = {na}  taxa-base = {na/n*100:.3f}%")
        print(f"AUROC = {a:.4f}  EP {se:.4f}  "
              f"IC95% [{a-1.96*se:.4f}, {a+1.96*se:.4f}]  "
              f"z = {(a-0.5)/se:.2f}")
        for k in (20, 50, 100):
            ac = sum(r for _, r in pares[:k])
            print(f"  precisao@{k:<4} = {ac}/{k} = {ac/k*100:5.2f}%  "
                  f"(ganho {ac/k/(na/n):.0f}x)")
    print(f"\n{len(linhas)} linhas em {SAIDA} ({(time.time()-t0)/60:.1f} min)")


if __name__ == "__main__":
    main()
