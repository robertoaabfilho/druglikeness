#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
validar_alvos_medidos.py

Mede o `prever_alvos.py` contra alvos MEDIDOS experimentalmente, em vez de
contra âncoras de literatura ou decoys aleatórios.

A pergunta: quando um composto tem bioatividade medida contra um alvo, o
modelo coloca esse alvo no topo do ranking?

ENTRADA
-------
Um CSV de evidência de bioatividade com, no mínimo, estas colunas (é o
formato que o pipeline de prospecção de plantas escreve em
`evidencias_chembl.csv`; os dois repositórios continuam separados -- aqui só
conhecemos os nomes das colunas):

    InChIKey          identidade do composto, para juntar com os SMILES
    alvo              nome preferencial do alvo, como o ChEMBL o escreve
    organismo_alvo    organismo do alvo (obrigatório: há ortólogos)
    pchembl           potência em escala log
    relacao           '=', '<', '<=', '~', '>' ...
    tipo_ensaio       'B' (ligação), 'F' (funcional) ...
    descricao_ensaio  usada para identificar triagem HTS

E um CSV de compostos com `InChIKey`, `smiles` e `composto_nome`.

POR QUE O ALVO É CASADO POR (nome, organismo)
---------------------------------------------
O arquivo de evidência não traz ID de alvo do ChEMBL, só o nome preferencial.
Nome sozinho é ambíguo -- 'Sodium-dependent noradrenaline transporter' existe
em humano e camundongo --, então a chave de junção é o PAR (nome, organismo).
Os dois lados vêm do ChEMBL, então a grafia coincide.

DOIS CUIDADOS QUE MUDAM O RESULTADO
-----------------------------------
1. TRIAGEM HTS NÃO É ATRIBUIÇÃO DE ALVO. Linhas de `PUBCHEM_BIOASSAY`/`qHTS`
   são varreduras de promiscuidade, muitas com veredito 'Inconclusive'.
   Medidas nelas, o recall@10 cai para 0,11; em ensaio dirigido, 0,27. Elas
   entram como estrato próprio, nunca somadas ao resto sem aviso.

2. DENSIDADE DE RÓTULOS. O recall@10 é limitado por construção quando o
   composto tem mais de 10 alvos medidos -- e aqui a mediana é 2, mas a cauda
   chega a 848 (painéis de triagem). Por isso reportamos também `hit@k`
   (algum alvo verdadeiro apareceu no top-k?) e a POSIÇÃO do alvo no ranking,
   que não depende da densidade.

LIMITE DE VAZAMENTO -- LEIA ANTES DE CITAR OS NÚMEROS
-----------------------------------------------------
O PP4CKD usou como rótulo positivo de treino todo par (composto, alvo) com
`pchembl >= 5.0` (relação =, <, <=, ~ ou nula) ou inibição >= 50%, em alvo
com >= 5 moléculas ativas. Logo, pares com `pchembl < 5` nunca foram
ensinados e seriam o conjunto limpo.

Só que a divisão de folds (`data/cv_splits_10fold.npz`) NÃO foi publicada --
não está no repositório nem no tarball do Zenodo --, e depois de remover as
linhas de HTS sobram 21 pares com pchembl < 5, amostra pequena demais.

Conclusão honesta: NÃO há como produzir uma estimativa livre de vazamento com
os artefatos publicados. Todo número aqui é LIMITE SUPERIOR. O que se pode
fazer, e é o que este script faz, é mostrar como o desempenho varia com o
treino do alvo e a qualidade do ensaio -- e isso é informativo mesmo com
vazamento, porque a variação é monotônica.

USO
---
    python mods/validar_alvos_medidos.py \
        --evidencias ~/prospector/.../evidencias_chembl.csv \
        --compostos  ~/prospector/.../atividades_compostos_parcial.csv \
        --rotulos    ~/modelos_baixados/pp4ckd/target_labels.tsv \
        --checkpoint ~/modelos_baixados/pp4ckd/models/ecfp4_fold_4.pt
"""

import argparse
import csv
import sys
from math import comb
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from prever_alvos import carregar_modelo, carregar_rotulos, fingerprint_de  # noqa: E402

try:
    import torch
except ImportError:
    print("[ERRO] Precisa do PyTorch: pip install torch", file=sys.stderr)
    sys.exit(1)

# Relações que o PP4CKD aceitou como atividade medida (não um limite inferior).
RELACOES_VALIDAS = {"=", "<", "<=", "~", ""}
# Limiar que o PP4CKD usou para transformar atividade em rótulo positivo.
PCHEMBL_TREINO = 5.0


# Heuristica unica de triagem HTS, compartilhada com a coleta de evidencias.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from mods.evidencias.medido import e_triagem  # noqa: E402,F401


def ler_evidencias(caminho: str, rotulos: list[dict]) -> tuple[dict, dict]:
    """-> ({inchikey: {idx_alvo: (pchembl, hts, tipo)}}, contadores)"""
    por_nome_org: dict[tuple[str, str], list[int]] = {}
    for i, r in enumerate(rotulos):
        chave = (r["nome"].casefold(), r["organismo"].casefold())
        por_nome_org.setdefault(chave, []).append(i)

    pares: dict[str, dict[int, tuple]] = {}
    c = {"linhas": 0, "sem_alvo": 0, "sem_pchembl": 0, "relacao_ruim": 0, "hts": 0}
    with open(Path(caminho).expanduser(), encoding="utf-8-sig", newline="") as f:
        for linha in csv.DictReader(f):
            c["linhas"] += 1
            nome = (linha.get("alvo") or "").strip()
            org = (linha.get("organismo_alvo") or "").strip()
            if not nome or nome == "Unchecked":
                c["sem_alvo"] += 1
                continue
            achados = por_nome_org.get((nome.casefold(), org.casefold()))
            if not achados:
                c["sem_alvo"] += 1
                continue
            try:
                pch = float(linha.get("pchembl") or "nan")
            except ValueError:
                pch = float("nan")
            if pch != pch:
                c["sem_pchembl"] += 1
                continue
            if (linha.get("relacao") or "").strip() not in RELACOES_VALIDAS:
                c["relacao_ruim"] += 1
                continue
            hts = e_triagem(linha.get("descricao_ensaio") or "")
            if hts:
                c["hts"] += 1
            tipo = (linha.get("tipo_ensaio") or "").strip().upper()
            ik = (linha.get("InChIKey") or "").strip()
            if not ik:
                continue
            d = pares.setdefault(ik, {})
            for i in achados:
                ant = d.get(i)
                # Guarda a evidência mais forte do par; no empate, prefere
                # ensaio dirigido a triagem.
                if ant is None or pch > ant[0] or (pch == ant[0] and ant[1] and not hts):
                    d[i] = (pch, hts, tipo)
    return pares, c


def ler_smiles(caminho: str, desejados: set[str], teto: int) -> dict[str, str]:
    smiles: dict[str, str] = {}
    with open(Path(caminho).expanduser(), encoding="utf-8-sig", newline="") as f:
        for linha in csv.DictReader(f):
            ik = (linha.get("InChIKey") or "").strip()
            if ik not in desejados or ik in smiles:
                continue
            s = (linha.get("smiles") or "").strip()
            if s and s.lower() not in ("nan", "none", "-", "null"):
                smiles[ik] = s
            if len(smiles) >= teto:
                break
    return smiles


def posicoes_no_ranking(iks, fps, rede, dispositivo, n_alvos, lote=512) -> dict:
    """{inchikey: vetor de posições, 0 = mais provável}"""
    ordem = {}
    with torch.no_grad():
        for i in range(0, len(fps), lote):
            x = torch.from_numpy(np.asarray(fps[i:i + lote], dtype="float32")).to(dispositivo)
            p = torch.sigmoid(rede(x))
            rk = torch.argsort(p, dim=1, descending=True).cpu().numpy()
            for j, ik in enumerate(iks[i:i + lote]):
                pos = np.empty(n_alvos, dtype=np.int32)
                pos[rk[j]] = np.arange(n_alvos)
                ordem[ik] = pos
    return ordem


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[3],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--evidencias", required=True)
    ap.add_argument("--compostos", required=True)
    ap.add_argument("--rotulos", required=True)
    ap.add_argument("--checkpoint", required=True)
    ap.add_argument("--fp", default="ECFP4")
    ap.add_argument("--max-compostos", type=int, default=6000)
    args = ap.parse_args()

    rotulos = carregar_rotulos(args.rotulos)
    n_alvos = len(rotulos)
    pares, c = ler_evidencias(args.evidencias, rotulos)
    print(f"[INFO] {c['linhas']:,} linhas de evidência lidas.")
    print(f"[INFO]   sem alvo resolvível no modelo : {c['sem_alvo']:,}")
    print(f"[INFO]   sem pchembl                   : {c['sem_pchembl']:,}")
    print(f"[INFO]   relação é limite inferior     : {c['relacao_ruim']:,}")
    print(f"[INFO]   triagem HTS/PubChem           : {c['hts']:,} (estrato próprio)")
    print(f"[INFO] {len(pares):,} compostos com alvo medido resolvido.")

    smiles = ler_smiles(args.compostos, set(pares), args.max_compostos)
    rede, cfg, dispositivo, ckpt = carregar_modelo(None, args.fp, checkpoint=args.checkpoint)

    iks, fps = [], []
    for ik, s in smiles.items():
        fp = fingerprint_de(s, cfg)
        if fp is not None:
            iks.append(ik)
            fps.append(fp)
    print(f"[INFO] {len(iks):,} compostos previstos com '{ckpt}' em '{dispositivo}'.\n")
    if not iks:
        print("[ERRO] Nenhum composto utilizável.", file=sys.stderr)
        sys.exit(1)

    ordem = posicoes_no_ranking(iks, fps, rede, dispositivo, n_alvos)

    densidade = [len(pares[ik]) for ik in iks]
    print(f"[DENSIDADE] alvos medidos por composto: média {np.mean(densidade):.1f}, "
          f"mediana {np.median(densidade):.0f}, máx {max(densidade)}")
    print("            (recall@10 é limitado por construção quando passa de 10)\n")

    def medir(nome, aceita):
        rec = {k: [0, 0] for k in (1, 5, 10)}
        hit = {k: [0, 0] for k in (1, 5, 10)}
        for ik in iks:
            verd = [i for i, (pch, hts, tipo) in pares[ik].items() if aceita(i, pch, hts, tipo)]
            if not verd:
                continue
            pos = ordem[ik]
            for k in rec:
                n = sum(1 for i in verd if pos[i] < k)
                rec[k][0] += n
                rec[k][1] += len(verd)
                hit[k][0] += (n > 0)
                hit[k][1] += 1
        if not rec[10][1]:
            print(f"{nome:<46} {'sem dados suficientes':>35}")
            return
        print(f"{nome:<46} {rec[10][1]:>7,} {rec[1][0]/rec[1][1]:>7.3f} "
              f"{rec[10][0]/rec[10][1]:>8.3f} {hit[10][0]/hit[10][1]:>8.3f}")

    print(f"{'estrato':<46} {'pares':>7} {'rec@1':>7} {'rec@10':>8} {'hit@10':>8}")
    print("-" * 80)
    medir("todos os pares com pchembl", lambda i, p, h, t: True)
    print("-- qualidade do ensaio --")
    medir("triagem HTS/PubChem (ruído)", lambda i, p, h, t: h)
    medir("ensaio dirigido (não-HTS)", lambda i, p, h, t: not h)
    medir("não-HTS, ensaio de ligação (tipo B)", lambda i, p, h, t: not h and t == "B")
    print("-- treino do alvo (não-HTS, pchembl>=6) --")
    for lo, hi, lab in [(0, 50, "< 50 moléculas ativas"), (50, 500, "50-500"),
                        (500, 5000, "500-5.000"), (5000, 10**9, ">= 5.000")]:
        medir(f"alvo com {lab}",
              lambda i, p, h, t, lo=lo, hi=hi: (not h and p >= 6
                                                and lo <= rotulos[i]["n_ativos"] < hi))
    print("-- organismo (não-HTS, pchembl>=6) --")
    medir("alvo humano", lambda i, p, h, t: not h and p >= 6
          and rotulos[i]["organismo"] == "Homo sapiens")
    medir("alvo não humano", lambda i, p, h, t: not h and p >= 6
          and rotulos[i]["organismo"] != "Homo sapiens")
    print("-- vazamento (ver docstring: amostra insuficiente) --")
    medir("não-HTS, pchembl<5 (nunca treinado)",
          lambda i, p, h, t: not h and p < PCHEMBL_TREINO)
    print("-" * 80)

    # Posição no ranking: a medida que não depende da densidade de rótulos.
    posicoes, melhores = [], []
    for ik in iks:
        ps = [int(ordem[ik][i]) + 1 for i, (p, h, t) in pares[ik].items()
              if not h and p >= 6]
        if ps:
            posicoes.extend(ps)
            melhores.append(min(ps))
    if posicoes:
        posicoes = np.array(posicoes)
        melhores = np.array(melhores)
        print(f"\nPosição do alvo medido entre os {n_alvos:,} alvos "
              f"(ensaio dirigido, pchembl>=6; {len(posicoes):,} pares):")
        for q in (10, 25, 50, 75, 90):
            print(f"  percentil {q:>2}: posição {np.percentile(posicoes, q):>7.0f}"
                  f"   (melhor alvo do composto: {np.percentile(melhores, q):>6.0f})")
        print(f"\n  mediana: posição {np.median(posicoes):.0f} de {n_alvos:,} "
              f"= top {100*np.median(posicoes)/n_alvos:.1f}%")
        print(f"  melhor alvo de cada composto: posição {np.median(melhores):.0f}")

        print("\nTraduzindo para o uso real (você fornece uma lista curta de alvos):")
        print("fração dos pares em que o alvo verdadeiro ficaria no top-3 DA LISTA")
        p_vence = (n_alvos - posicoes) / (n_alvos - 1)
        for tam in (10, 20, 50, 100):
            p3 = np.mean([sum(comb(tam - 1, k) * ((1 - pv) ** k) * (pv ** (tam - 1 - k))
                              for k in range(3)) for pv in p_vence])
            print(f"  lista de {tam:>3} alvos -> {100*p3:>5.1f}%")
        print("\n(aproximação: supõe que os outros alvos da lista se comportam como")
        print(" alvos sorteados do modelo. Uma lista toda de alvos plausíveis para a")
        print(" mesma doença é mais difícil que isso.)")


if __name__ == "__main__":
    main()
