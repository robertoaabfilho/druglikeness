#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
comparar_metodos_ranking.py

Compara empiricamente o método de ranking usado neste projeto contra
âncoras consolidadas da literatura de quimioinformática/MCDM:

  1. 'perfil'          -- similaridade ponderada à média das referências
                           (método próprio, gerar_ranking em ranking_multibase.py)
  2. 'topsis_classico'  -- TOPSIS PADRÃO: ideal/anti-ideal vêm do próprio
                           conjunto de candidatos (farmaco.csv), com direção
                           benefício/custo por critério definida pela
                           literatura de Lipinski/Veber. Implementado aqui
                           do zero, de forma independente (não depende de
                           nenhum "TOPSIS híbrido" -- esse foi removido do
                           ranking_multibase.py por não ser mais defensável
                           como TOPSIS reconhecível; ver histórico do projeto).
  3. 'qed'              -- Quantitative Estimate of Druglikeness
                           (Bickerton et al., Nature Chemistry 2012),
                           via RDKit (rdkit.Chem.QED). CHECAGEM DE
                           COERÊNCIA, não um comparador direto -- mede um
                           construto diferente (drug-likeness geral).
  4. 'tanimoto_max'      -- maior similaridade de Tanimoto (fingerprint de
                           Morgan/ECFP4, raio 2, 2048 bits) entre o composto
                           e QUALQUER droga de referência. Também uma
                           CHECAGEM DE COERÊNCIA (proximidade estrutural),
                           não um comparador direto.

Para cada par de métodos, calcula a correlação de Spearman entre as
pontuações. Ao interpretar os resultados: correlação alta com
'topsis_classico' é evidência direta de robustez (mesmo tipo de problema,
resolvido de duas formas). Correlação com 'qed'/'tanimoto_max' deve ser
lida como coerência, não como validação -- baixa correlação não invalida
o método por si só, já que medem construtos diferentes.

USO
---
    python comparar_metodos_ranking.py \
        --farmaco farmaco.csv \
        --referencia ranking.referencia.csv \
        --pesos ranking.pesos.json \
        --output comparacao_metodos.csv

REQUISITOS
----------
    pip install rdkit numpy scipy
    ranking_multibase.py precisa estar na mesma pasta (reaproveita
    gerar_ranking e a leitura do farmaco.csv de lá).
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr
from rdkit import Chem
from rdkit.Chem import QED, AllChem, DataStructs

try:
    from ranking_multibase import (
        ler_farmaco_csv,
        gerar_ranking,
        _para_float,
        _para_bool,
    )
except ImportError:
    print(
        "[ERRO] Não foi possível importar 'ranking_multibase.py'. "
        "Coloque este script na mesma pasta.",
        file=sys.stderr,
    )
    sys.exit(1)


# =========================================================================
# TOPSIS clássico (ideal vem do PRÓPRIO farmaco.csv) -- escopo restrito aos
# descritores de Lipinski/Veber, cuja direção benefício/custo é consenso
# estabelecido na literatura, evitando qualquer heurística de classificação
# ad-hoc para as colunas ADMET dinâmicas.
# =========================================================================

# 'custo': quanto MENOR, mais "drug-like" (Lipinski/Veber). 'beneficio':
# quanto MAIOR (aqui, codificado como 1.0), melhor.
DIRECAO_LIPINSKI_VEBER = {
    "PesoMolecular": "custo",
    "DoadoresHB": "custo",
    "AceptoresHB": "custo",
    "TPSA": "custo",
    "LigacoesRotacionaveis": "custo",
    "RefratividadeMolar": "custo",
    "AtomosPesados": "custo",
    "Passa_Ro5_parcial": "beneficio",
    "Passa_Veber": "beneficio",
}


def _valor_topsis_classico(linha: dict, coluna: str) -> float | None:
    if coluna in ("Passa_Ro5_parcial", "Passa_Veber"):
        v = _para_bool(linha.get(coluna))
        return None if v is None else (1.0 if v else 0.0)
    return _para_float(linha.get(coluna))


def calcular_topsis_classico(linhas_farmaco: list[dict], pesos: dict[str, dict]) -> list[dict]:
    criterios = [c for c in DIRECAO_LIPINSKI_VEBER if c in pesos]
    if not criterios:
        print("[AVISO] Nenhum descritor de Lipinski/Veber encontrado nos pesos -- TOPSIS clássico usará peso 1.0 para todos.")
        criterios = list(DIRECAO_LIPINSKI_VEBER)

    matriz = []
    for linha in linhas_farmaco:
        vetor = [_valor_topsis_classico(linha, c) for c in criterios]
        matriz.append(vetor)

    matriz = np.array(matriz, dtype=float)
    # Imputa ausentes com a média da própria coluna (não há "alvo externo" no TOPSIS clássico).
    for j in range(matriz.shape[1]):
        coluna_vals = matriz[:, j]
        validos = coluna_vals[~np.isnan(coluna_vals)]
        media = validos.mean() if len(validos) else 0.0
        matriz[np.isnan(matriz[:, j]), j] = media

    normas = np.sqrt((matriz ** 2).sum(axis=0))
    normas[normas == 0] = 1.0
    r = matriz / normas

    pesos_vetor = np.array([pesos.get(c, {}).get("peso", 1.0) for c in criterios])
    v = r * pesos_vetor

    ideal_mais = np.zeros(len(criterios))
    ideal_menos = np.zeros(len(criterios))
    for j, c in enumerate(criterios):
        if DIRECAO_LIPINSKI_VEBER[c] == "beneficio":
            ideal_mais[j], ideal_menos[j] = v[:, j].max(), v[:, j].min()
        else:  # custo
            ideal_mais[j], ideal_menos[j] = v[:, j].min(), v[:, j].max()

    d_mais = np.sqrt(((v - ideal_mais) ** 2).sum(axis=1))
    d_menos = np.sqrt(((v - ideal_menos) ** 2).sum(axis=1))
    soma = d_mais + d_menos
    proximidade = np.divide(d_menos, soma, out=np.zeros_like(d_menos), where=soma != 0)

    resultado = []
    for i, linha in enumerate(linhas_farmaco):
        resultado.append({
            "Nome": linha.get("Nome", ""),
            "SMILES": linha.get("SMILES", ""),
            "Pontuacao_Final": round(float(proximidade[i]) * 100, 2),
        })
    return resultado


# =========================================================================
# QED (Bickerton et al. 2012) -- via RDKit
# =========================================================================

def calcular_qed(linhas_farmaco: list[dict]) -> dict[str, float]:
    resultado = {}
    for linha in linhas_farmaco:
        smiles = linha.get("SMILES", "")
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        try:
            resultado[linha.get("Nome", "")] = round(QED.qed(mol) * 100, 2)  # escala 0-100, como os outros métodos
        except Exception:
            continue
    return resultado


# =========================================================================
# Tanimoto máximo contra as drogas de referência (fingerprint de Morgan/ECFP4)
# =========================================================================

def _fingerprint(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)


def calcular_tanimoto_max(linhas_farmaco: list[dict], nomes_smiles_referencia: list[tuple[str, str]]) -> dict[str, float]:
    fps_referencia = []
    for _nome, smiles in nomes_smiles_referencia:
        fp = _fingerprint(smiles)
        if fp is not None:
            fps_referencia.append(fp)

    if not fps_referencia:
        print("[AVISO] Nenhum fingerprint válido nas drogas de referência -- Tanimoto não pôde ser calculado.")
        return {}

    resultado = {}
    for linha in linhas_farmaco:
        fp = _fingerprint(linha.get("SMILES", ""))
        if fp is None:
            continue
        similaridades = DataStructs.BulkTanimotoSimilarity(fp, fps_referencia)
        resultado[linha.get("Nome", "")] = round(max(similaridades) * 100, 2)  # escala 0-100
    return resultado


# =========================================================================
# Leitura auxiliar
# =========================================================================

def ler_referencia_nomes_smiles(caminho: str) -> list[tuple[str, str]]:
    delimitador = ";"
    with open(caminho, "r", encoding="utf-8-sig", newline="") as f:
        primeira_linha = f.readline()
        f.seek(0)
        if "," in primeira_linha and ";" not in primeira_linha:
            delimitador = ","
        reader = csv.DictReader(f, delimiter=delimitador)
        return [(l["Nome"], l["SMILES"]) for l in reader if l.get("Nome") and l.get("SMILES")]


# =========================================================================
# Montagem da tabela comparativa + correlação de Spearman
# =========================================================================

def montar_tabela(resultados_por_metodo: dict[str, dict[str, float]]) -> tuple[list[str], dict[str, dict[str, float]]]:
    """resultados_por_metodo = {metodo: {nome_composto: pontuacao}}"""
    todos_nomes = sorted(set().union(*[set(d.keys()) for d in resultados_por_metodo.values()]))
    tabela = {nome: {} for nome in todos_nomes}
    for metodo, dados in resultados_por_metodo.items():
        for nome in todos_nomes:
            tabela[nome][metodo] = dados.get(nome)
    return todos_nomes, tabela


def salvar_tabela_csv(nomes: list[str], tabela: dict[str, dict], metodos: list[str], caminho: str) -> None:
    with open(caminho, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter=";")
        writer.writerow(["Nome"] + metodos)
        for nome in nomes:
            writer.writerow([nome] + [tabela[nome].get(m, "") for m in metodos])
    print(f"[INFO] Tabela comparativa salva em: {caminho}")


def calcular_matriz_spearman(nomes: list[str], tabela: dict[str, dict], metodos: list[str]) -> dict[tuple[str, str], tuple[float, float]]:
    correlacoes = {}
    for i, m1 in enumerate(metodos):
        for m2 in metodos[i + 1:]:
            pares = [(tabela[n][m1], tabela[n][m2]) for n in nomes if tabela[n].get(m1) is not None and tabela[n].get(m2) is not None]
            if len(pares) < 3:
                correlacoes[(m1, m2)] = (float("nan"), float("nan"))
                continue
            x, y = zip(*pares)
            rho, p_valor = spearmanr(x, y)
            correlacoes[(m1, m2)] = (rho, p_valor)
    return correlacoes


def imprimir_matriz_spearman(correlacoes: dict, metodos: list[str]) -> None:
    print("\n[INFO] Correlação de Spearman entre os métodos (rho, p-valor):")
    print("       ('topsis_classico' é o comparador direto; 'qed'/'tanimoto_max' são checagens de coerência, não comparadores diretos.)")
    for (m1, m2), (rho, p) in correlacoes.items():
        rho_str = f"{rho:.3f}" if rho == rho else "N/D"  # NaN check
        p_str = f"{p:.4f}" if p == p else "N/D"
        print(f"         {m1:16s} x {m2:16s}   rho={rho_str}   p={p_str}")


def main():
    parser = argparse.ArgumentParser(description="Compara os métodos de ranking (perfil, TOPSIS clássico, QED, Tanimoto).")
    parser.add_argument("--farmaco", default="farmaco.csv", help="CSV de compostos candidatos.")
    parser.add_argument("--referencia", default="ranking.referencia.csv", help="CSV de drogas de referência (saída do ranking_multibase.py).")
    parser.add_argument("--pesos", default="ranking.pesos.json", help="JSON de pesos por característica (saída do ranking_multibase.py).")
    parser.add_argument("--output", default="comparacao_metodos.csv", help="CSV comparativo de saída.")
    args = parser.parse_args()

    for caminho, nome in [(args.farmaco, "farmaco"), (args.referencia, "referência"), (args.pesos, "pesos")]:
        if not Path(caminho).is_file():
            print(f"[ERRO] Arquivo de {nome} não encontrado: {caminho}", file=sys.stderr)
            sys.exit(1)

    _colunas, linhas_farmaco = ler_farmaco_csv(args.farmaco)
    nomes_smiles_referencia = ler_referencia_nomes_smiles(args.referencia)
    with open(args.pesos, "r", encoding="utf-8") as f:
        pesos = json.load(f)

    print(f"[INFO] {len(linhas_farmaco)} composto(s) candidato(s), {len(nomes_smiles_referencia)} droga(s) de referência, {len(pesos)} peso(s) carregado(s).")

    # ---- 4 métodos ----
    print("[INFO] Calculando método 'perfil'...")
    ranking_perfil = gerar_ranking(linhas_farmaco, pesos)
    perfil = {r["Nome"]: r["Pontuacao_Final"] for r in ranking_perfil}

    print("[INFO] Calculando método 'topsis_classico'...")
    ranking_topsis_c = calcular_topsis_classico(linhas_farmaco, pesos)
    topsis_classico = {r["Nome"]: r["Pontuacao_Final"] for r in ranking_topsis_c}

    print("[INFO] Calculando 'qed'...")
    qed = calcular_qed(linhas_farmaco)

    print("[INFO] Calculando 'tanimoto_max'...")
    tanimoto_max = calcular_tanimoto_max(linhas_farmaco, nomes_smiles_referencia)

    resultados_por_metodo = {
        "perfil": perfil,
        "topsis_classico": topsis_classico,
        "qed": qed,
        "tanimoto_max": tanimoto_max,
    }
    metodos = list(resultados_por_metodo.keys())

    nomes, tabela = montar_tabela(resultados_por_metodo)
    salvar_tabela_csv(nomes, tabela, metodos, args.output)

    correlacoes = calcular_matriz_spearman(nomes, tabela, metodos)
    imprimir_matriz_spearman(correlacoes, metodos)

    caminho_corr = str(Path(args.output).with_suffix("")) + ".spearman.csv"
    with open(caminho_corr, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter=";")
        writer.writerow(["Metodo_A", "Metodo_B", "Rho_Spearman", "P_Valor"])
        for (m1, m2), (rho, p) in correlacoes.items():
            writer.writerow([m1, m2, rho, p])
    print(f"[INFO] Matriz de correlação salva em: {caminho_corr}")


if __name__ == "__main__":
    main()
