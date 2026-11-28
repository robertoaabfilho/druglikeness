"""chembl_multitask_model (ChEMBL 36, ONNX, oficial do EMBL-EBI) no mesmo
arame de LIT-PCBA usado para o PP4CKD.

Comparacao limpa: mesma versao do ChEMBL (36), mesmo dataset, mesmas metricas,
mesma subamostra de inativos. ECFP4 raio 2, 1024 bits, como o README do modelo
especifica.
"""
import csv
import random
import sys
import zlib
import time
from pathlib import Path

import caminhos_dados as cd_

sys.path.insert(0, str(Path(__file__).parent))

LITPCBA = cd_.de_ambiente("LITPCBA_FULL")
MT = Path(__file__).parent / "mt/trained_models/chembl_36_model"
ONNX = MT / "chembl_36_multitask.onnx"
ALVOS_JSON = MT / "targets_36_all.json"
SAIDA = Path(__file__).parent / "resultado_multitask.tsv"

SUBAMOSTRA = 60_000
LOTE = 1024
SEMENTE = 20261127          # a mesma do benchmark do PP4CKD
FP_SIZE, RADIUS = 1024, 2

MAPA = {
    "ADRB2": "CHEMBL210", "ALDH1": "CHEMBL3577", "ESR1_ago": "CHEMBL206",
    "ESR1_ant": "CHEMBL206", "FEN1": "CHEMBL5027", "GBA": "CHEMBL2179",
    "IDH1": "CHEMBL2007625", "KAT2A": "CHEMBL5501", "MAPK1": "CHEMBL4040",
    "MTORC1": "CHEMBL2842", "OPRK1": "CHEMBL237", "PKM2": "CHEMBL1075189",
    "PPARG": "CHEMBL235", "TP53": "CHEMBL4096", "VDR": "CHEMBL1977",
}

import json
import multiprocessing as mp

import numpy as np
import onnxruntime
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
RDLogger.DisableLog("rdApp.*")

import metricas_vs as mv


def ler_smi(caminho):
    out = []
    with open(caminho) as fh:
        for linha in fh:
            p = linha.split()
            if p:
                out.append(p[0])
    return out


def _um_alvo(args):
    """Roda um alvo num processo proprio. Devolve a linha de metricas."""
    alvo, cid = args
    alvos_mt = json.load(open(ALVOS_JSON))
    if cid not in set(alvos_mt):
        return None
    opc = onnxruntime.SessionOptions()
    opc.intra_op_num_threads = 1          # 4 processos, 1 thread cada
    opc.inter_op_num_threads = 1
    ses = onnxruntime.InferenceSession(
        str(ONNX), sess_options=opc, providers=["CPUExecutionProvider"])
    nome_in = ses.get_inputs()[0].name

    from rdkit import Chem as _C, RDLogger as _R
    from rdkit.Chem import rdMolDescriptors as _rd
    _R.DisableLog("rdApp.*")

    def escore(smiles):
        mol = _C.MolFromSmiles(smiles)
        if mol is None:
            return None
        f = np.array(_rd.GetMorganFingerprintAsBitVect(
            mol, RADIUS, nBits=FP_SIZE), dtype=np.float32)
        return float(ses.run([cid], {nome_in: f})[0][0])

    d = LITPCBA / alvo
    ativos = ler_smi(d / "actives.smi")
    inativos = ler_smi(d / "inactives.smi")
    # ATENCAO: esta subamostra NAO e a mesma lista do benchmark do PP4CKD.
    # La a semente era global e o estado do gerador era consumido em
    # sequencia, alvo a alvo; aqui cada processo sorteia do zero. Comparar
    # contra a coluna "completo" do PP4CKD, que e o numero sem subamostra --
    # e que bate com a subamostra dele na terceira decimal, mostrando que a
    # escolha da subamostra nao move o resultado.
    # zlib.crc32 e estavel entre execucoes; hash() de str nao e (PYTHONHASHSEED).
    rnd = random.Random(SEMENTE + zlib.crc32(alvo.encode()))
    idx = (rnd.sample(range(len(inativos)), SUBAMOSTRA)
           if len(inativos) > SUBAMOSTRA else list(range(len(inativos))))

    t = time.time()
    esc, rot = [], []
    for s_ in ativos:
        v = escore(s_)
        if v is not None:
            esc.append(v); rot.append(1)
    for i in idx:
        v = escore(inativos[i])
        if v is not None:
            esc.append(v); rot.append(0)
    r = mv.resumo(esc, rot)
    r.update(alvo=alvo, chembl=cid, metodo="chembl_mt36",
             conjunto=f"sub{SUBAMOSTRA}", segundos=round(time.time() - t, 1))
    return r


def main():
    t0 = time.time()
    alvos_mt = set(json.load(open(ALVOS_JSON)))
    tarefas = [(a, MAPA[a]) for a in sorted(MAPA)]
    fora = [a for a, c in tarefas if c not in alvos_mt]
    if fora:
        print(f"fora do painel de 786, pulados: {', '.join(fora)}", flush=True)
    tarefas = [(a, c) for a, c in tarefas if c in alvos_mt]
    print(f"{len(tarefas)} alvos em 4 processos", flush=True)

    with mp.Pool(4) as pool:
        linhas = [r for r in pool.imap_unordered(_um_alvo, tarefas)
                  if r is not None]
    linhas.sort(key=lambda r: r["alvo"])
    for r in linhas:
        print(f"{r['alvo']:<10} {r['chembl']:<14} AUROC={r['auroc']:.3f} "
              f"EF1%={r['ef_1%']:6.2f} BEDROC={r['bedroc_20.5']:.3f} "
              f"({r['segundos']}s)", flush=True)
    campos = ["alvo", "chembl", "metodo", "conjunto", "n", "n_ativos",
              "auroc", "ef_0.5%", "ef_1%", "ef_5%", "bedroc_20.5",
              "bedroc_0.85", "segundos"]
    with open(SAIDA, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=campos, delimiter="\t")
        w.writeheader(); w.writerows(linhas)
    print(f"\n{len(linhas)} linhas em {SAIDA} "
          f"({(time.time()-t0)/60:.1f} min)", flush=True)


if __name__ == "__main__":
    main()
