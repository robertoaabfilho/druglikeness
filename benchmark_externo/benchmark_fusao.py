"""Fusao de canais e dominio de aplicabilidade, medidos por PRECISAO NO TOPO.

Objetivo: o que economiza bancada e precisao@k para k pequeno, nao AUC sobre
centenas de milhares de compostos.

Divisao: LIT-PCBA AVE-unbiased. O canal de similaridade usa `active_T` como
referencia e e avaliado em `active_V` + `inactive_V` -- sem o leave-one-out
otimista do benchmark anterior, em que um ativo de teste via outros ativos de
teste.

Canais (todos sobre o conjunto V):
  p_alvo  PP4CKD, probabilidade na cabeca do alvo. Zero-shot.
  sim_T   maior Tanimoto ECFP4 aos ativos de TREINO. Honesto.
  famil   maior probabilidade do PP4CKD entre os 7.676 alvos -- proxy de
          "o modelo reconhece esta molecula", candidato a dominio de
          aplicabilidade.

Fusoes:
  media_postos   media dos postos percentuais de p_alvo e sim_T
  rrf            reciprocal rank fusion, k = 60 (padrao da area)
  produto        produto dos postos percentuais
  maximo         maior dos dois postos percentuais

Filtro de dominio: descarta a metade de menor `famil` e remede. Referencia:
omitir a metade de menor aplicabilidade melhora a triagem (Fechner et al.).
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
SAIDA = Path(__file__).parent / "resultado_fusao.tsv"

LOTE = 2048
MIN_ATIVOS_V = 15          # alvo com menos que isso nao sustenta precisao@20

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


def ler_smi(caminho):
    out = []
    with open(caminho) as fh:
        for linha in fh:
            p = linha.split()
            if p:
                out.append(p[0])
    return out


def postos_percentuais(v):
    """Posto percentual (0 = pior, 1 = melhor), empates com posto medio."""
    v = np.asarray(v, dtype=np.float64)
    ordem = np.argsort(v, kind="mergesort")
    postos = np.empty(len(v), dtype=np.float64)
    postos[ordem] = np.arange(len(v), dtype=np.float64)
    # media para empates
    i = 0
    vs = v[ordem]
    while i < len(v):
        j = i
        while j + 1 < len(v) and vs[j + 1] == vs[i]:
            j += 1
        if j > i:
            postos[ordem[i:j + 1]] = (i + j) / 2.0
        i = j + 1
    return postos / max(1, len(v) - 1)


def precisao_em_k(escores, rotulos, k):
    pares = sorted(zip(escores, rotulos), key=lambda p: (-p[0], p[1]))
    topo = pares[:k]
    return sum(r for _, r in topo) / len(topo) if topo else float("nan")


def main():
    t0 = time.time()
    rede, cfg, _d, nome_ckpt = pa.carregar_modelo(
        None, fp="ECFP4", checkpoint=CKPT, dispositivo="cpu")
    rotulos = pa.carregar_rotulos(ROTULOS)
    col = {r["chembl_id"]: i for i, r in enumerate(rotulos) if r["chembl_id"]}
    torch.set_num_threads(4)
    gen = rdFingerprintGenerator.GetMorganGenerator(
        radius=cfg["raio"], fpSize=cfg["bits"])
    print(f"checkpoint {nome_ckpt}", flush=True)

    def pontuar(smiles_lista, coluna):
        """Devolve (p_alvo, famil, fingerprints) alinhados, pulando ilegiveis.
        Tambem devolve os indices mantidos."""
        p_alvo, famil, fps, mantidos = [], [], [], []
        bi, bv = [], []

        def descarregar():
            if not bi:
                return
            x = torch.tensor(np.array(bv, dtype=np.float32))
            with torch.no_grad():
                p = torch.sigmoid(rede(x)).numpy()
            p_alvo.extend(p[:, coluna].tolist())
            famil.extend(p.max(axis=1).tolist())
            bi.clear(); bv.clear()

        for i, s in enumerate(smiles_lista):
            mol = Chem.MolFromSmiles(s)
            if mol is None:
                continue
            mantidos.append(i)
            fps.append(gen.GetFingerprint(mol))
            bi.append(i)
            bv.append(gen.GetFingerprintAsNumPy(mol))
            if len(bi) >= LOTE:
                descarregar()
        descarregar()
        return np.array(p_alvo), np.array(famil), fps, mantidos

    linhas = []
    for alvo in sorted(MAPA):
        cid = MAPA[alvo]
        d = AVE / alvo
        at_T = ler_smi(d / "active_T.smi")
        at_V = ler_smi(d / "active_V.smi")
        in_V = ler_smi(d / "inactive_V.smi")
        if len(at_V) < MIN_ATIVOS_V:
            print(f"=== {alvo}: {len(at_V)} ativos em V (< {MIN_ATIVOS_V}), "
                  f"pulado", flush=True)
            continue
        coluna = col[cid]
        t = time.time()

        fps_T = [gen.GetFingerprint(m) for m in
                 (Chem.MolFromSmiles(s) for s in at_T) if m is not None]
        p_at, f_at, fp_at, _ = pontuar(at_V, coluna)
        p_in, f_in, fp_in, _ = pontuar(in_V, coluna)

        p_alvo = np.concatenate([p_at, p_in])
        famil = np.concatenate([f_at, f_in])
        rot = np.array([1] * len(p_at) + [0] * len(p_in))
        sim_T = np.array([max(DataStructs.BulkTanimotoSimilarity(f, fps_T))
                          for f in (fp_at + fp_in)])

        r_p, r_s = postos_percentuais(p_alvo), postos_percentuais(sim_T)
        n = len(rot)
        # RRF sobre os postos decrescentes (posicao 1 = melhor)
        pos_p = n - np.argsort(np.argsort(p_alvo, kind="mergesort"))
        pos_s = n - np.argsort(np.argsort(sim_T, kind="mergesort"))
        canais = {
            "p_alvo": p_alvo,
            "sim_T": sim_T,
            "famil": famil,
            "media_postos": (r_p + r_s) / 2.0,
            "rrf": 1.0 / (60 + pos_p) + 1.0 / (60 + pos_s),
            "produto": r_p * r_s,
            "maximo": np.maximum(r_p, r_s),
        }

        corte = np.median(famil)
        dentro = famil >= corte
        for nome, esc in canais.items():
            for marca, sel in (("tudo", np.ones(n, bool)), ("dominio", dentro)):
                if marca == "dominio" and nome == "famil":
                    continue
                e, rr = esc[sel].tolist(), rot[sel].tolist()
                if sum(rr) == 0:
                    continue
                res = mv.resumo(e, rr)
                res.update(
                    alvo=alvo, chembl=cid, canal=nome, filtro=marca,
                    ativos_T=len(fps_T),
                    **{f"prec_{k}": precisao_em_k(e, rr, k)
                       for k in (20, 50, 100)})
                linhas.append(res)
        print(f"=== {alvo}: {len(at_V)} ativos V, {len(in_V)} inativos V, "
              f"{len(fps_T)} ativos T ({time.time()-t:.0f}s)", flush=True)

        campos = ["alvo", "chembl", "canal", "filtro", "n", "n_ativos",
                  "ativos_T", "prec_20", "prec_50", "prec_100",
                  "auroc", "ef_1%", "bedroc_20.5"]
        with open(SAIDA, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=campos, delimiter="\t",
                               extrasaction="ignore")
            w.writeheader(); w.writerows(linhas)

    print(f"\n{len(linhas)} linhas em {SAIDA} "
          f"({(time.time()-t0)/60:.1f} min)", flush=True)


if __name__ == "__main__":
    main()
