#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
main.py

Orquestrador do pipeline de drogabilidade. Lê 'config.txt' e roda as
etapas em sequência, cada uma como um subprocesso do script correspondente
em mods/. Todos os arquivos gerados (CSVs, JSON de pesos, etc.) vão para
a pasta results/ (criada automaticamente se não existir).

ESTRUTURA DE PASTAS ESPERADA
------------------------------
    ./main.py
    ./config.txt          <- você edita
    ./input.txt            <- você edita (lista de SMILES)
    ./mods/                <- todos os scripts do pipeline
        configuracao.py
        farmaco_completo.py
        ranking_multibase.py
        comparar_metodos_ranking.py
        validar_enriquecimento.py
        ranking_embeddings.py
        calibrar_deeppk.py
    ./results/              <- criada automaticamente; todas as saídas vão aqui
        farmaco.csv
        ranking.csv
        ranking.referencia.csv
        ranking.pesos.json
        ...

ETAPAS (controladas pelo config.txt, seção "Etapas do pipeline")
-------------------------------------------------------------------
  1. mods/farmaco_completo.py   -- input.txt -> results/farmaco.csv     (rodar_farmaco)
  2. mods/ranking_multibase.py  -- sintoma(s) -> results/ranking.csv    (rodar_ranking)
     Pergunta interativamente se você quer remover alguma droga de
     referência antes de prosseguir (a menos que 'pular_revisao_referencias
     = true' no config.txt).
  3. mods/comparar_metodos_ranking.py -- validação vs. TOPSIS/QED/Tanimoto (rodar_comparacao_metodos)
  4. mods/validar_enriquecimento.py   -- validação leave-one-out + decoys (rodar_validacao_enriquecimento)
  5. mods/ranking_embeddings.py       -- ranking alternativo via ChemBERTa (rodar_embeddings)
     Verifica automaticamente se 'transformers' e 'torch' estão instalados
     antes de rodar essa etapa, e instala via pip se faltarem (são pesados
     -- alguns GB de download na primeira vez -- por isso essa etapa fica
     desativada por padrão no config.txt).

Cada etapa só roda se a anterior (da qual ela depende) tiver terminado
com sucesso. Se uma etapa falhar, o pipeline para ali e reporta
claramente qual etapa falhou.

USO
---
    python main.py                      # usa config.txt na pasta atual
    python main.py --config outro.txt   # usa outro arquivo de config

Se 'config.txt' não existir, main.py cria um modelo comentado e pede para
o usuário preenchê-lo antes de rodar de novo.
"""

import argparse
import importlib
import subprocess
import sys
from pathlib import Path

MODS_DIR = Path(__file__).parent / "mods"
sys.path.insert(0, str(MODS_DIR))

from configuracao import carregar_config, get_str, get_int, get_float, get_bool, get_list  # noqa: E402

MODELO_CONFIG = """\
# config.txt -- configuração do pipeline de drogabilidade
# Linhas começando com '#' são comentários. Formato: chave = valor

# --- Sintoma(s)/condição a pesquisar nas bases de referência (em inglês) ---
# Pode ser mais de um, separado por vírgula.
sintomas = neuropathic pain, peripheral neuropathy

# --- Pasta onde todos os resultados são salvos ---
results_dir = results

# --- farmaco_completo.py (input.txt -> farmaco.csv) ---
farmaco_input = input.txt
farmaco_output = farmaco.csv
pred_type = admet
batch_size = 5

# --- ranking_multibase.py (sintoma -> ranking.csv) ---
ranking_output = ranking.csv
top_por_fonte = 10
no_admet = false
pular_revisao_referencias = false

# --- validar_enriquecimento.py ---
n_decoys = 50

# --- ranking_embeddings.py (opcional -- requer transformers+torch, pesado) ---
modelo_embeddings = seyonec/ChemBERTa-zinc-base-v1
similaridade_embeddings = centroide
peso_embedding = 0.5
peso_admet_embeddings = 0.5

# --- Quais etapas rodar (main.py) ---
rodar_farmaco = true
rodar_ranking = true
rodar_comparacao_metodos = false
rodar_validacao_enriquecimento = false
rodar_embeddings = false
"""


def criar_config_modelo(caminho: str) -> None:
    with open(caminho, "w", encoding="utf-8") as f:
        f.write(MODELO_CONFIG)
    print(f"[INFO] '{caminho}' não existia -- criei um modelo. Edite-o e rode de novo.")


def verificar_e_instalar_dependencias(pacotes: dict[str, str]) -> bool:
    """
    'pacotes' é {nome_do_modulo_para_import: nome_do_pacote_pip}, ex:
    {'torch': 'torch', 'transformers': 'transformers'}. Verifica se cada
    módulo já está instalado (importável); se não, instala via pip
    automaticamente. Retorna True se, ao final, todos os módulos estão
    disponíveis (já estavam, ou foram instalados com sucesso).
    """
    faltando = []
    for modulo, pacote_pip in pacotes.items():
        try:
            importlib.import_module(modulo)
        except ImportError:
            faltando.append((modulo, pacote_pip))

    if not faltando:
        print(f"[INFO] Dependência(s) já instalada(s): {', '.join(pacotes.values())}.")
        return True

    print(f"[INFO] Dependência(s) faltando: {', '.join(p for _, p in faltando)}. Instalando via pip (pode demorar alguns minutos)...")
    for modulo, pacote_pip in faltando:
        print(f"[INFO]   Instalando '{pacote_pip}'...")
        resultado = subprocess.run([sys.executable, "-m", "pip", "install", pacote_pip])
        if resultado.returncode != 0:
            print(f"[ERRO] Falha ao instalar '{pacote_pip}'. Instale manualmente: pip install {pacote_pip}")
            return False

    print("[INFO] Dependência(s) instalada(s) com sucesso.")
    return True


def verificar_dependencias_opcionais(pacotes: dict[str, str]) -> None:
    """
    Como verificar_e_instalar_dependencias, mas para dependências cuja
    ausência não deveria travar o pipeline -- o próprio script sabe
    degradar graciosamente sem elas (ex: psycopg2 para a fonte
    DrugCentral). Tenta instalar; se falhar, só avisa e segue.
    """
    for modulo, pacote_pip in pacotes.items():
        try:
            importlib.import_module(modulo)
        except ImportError:
            print(f"[INFO] Dependência opcional '{pacote_pip}' não encontrada. Tentando instalar (não é obrigatória)...")
            resultado = subprocess.run([sys.executable, "-m", "pip", "install", pacote_pip])
            if resultado.returncode != 0:
                print(f"[AVISO] Não foi possível instalar '{pacote_pip}' automaticamente. A funcionalidade associada será pulada (o script já trata essa ausência).")


# Dependências obrigatórias e opcionais de cada etapa, mapeadas
# {nome_do_módulo_para_import: nome_do_pacote_pip}.
DEPENDENCIAS_OBRIGATORIAS = {
    "farmaco": {"rdkit": "rdkit", "requests": "requests"},
    "ranking": {"rdkit": "rdkit", "requests": "requests"},
    "comparacao": {"rdkit": "rdkit", "numpy": "numpy", "scipy": "scipy"},
    "validacao": {"numpy": "numpy", "scipy": "scipy", "requests": "requests"},
    "embeddings": {"numpy": "numpy", "torch": "torch", "transformers": "transformers"},
}
DEPENDENCIAS_OPCIONAIS = {
    "ranking": {"psycopg2": "psycopg2-binary"},  # só para a fonte DrugCentral; degrada sozinha sem ela
}


def rodar_etapa(nome_etapa: str, comando: list[str]) -> bool:
    """Roda uma etapa como subprocesso, herdando stdin/stdout/stderr (para
    que prompts interativos, como a revisão de referências do
    ranking_multibase.py, continuem funcionando normalmente).
    Retorna True se terminou com sucesso (código de saída 0)."""
    print("\n" + "=" * 70)
    print(f"ETAPA: {nome_etapa}")
    print("=" * 70)
    print(f"[INFO] Comando: {' '.join(comando)}\n")

    try:
        resultado = subprocess.run(comando)
    except FileNotFoundError:
        print(f"[ERRO] Não foi possível executar '{comando[1]}' -- o arquivo existe em mods/?")
        return False

    if resultado.returncode != 0:
        print(f"\n[ERRO] Etapa '{nome_etapa}' terminou com código {resultado.returncode}. Pipeline interrompido.")
        return False

    print(f"\n[INFO] Etapa '{nome_etapa}' concluída com sucesso.")
    return True


def main():
    parser = argparse.ArgumentParser(description="Orquestrador do pipeline de drogabilidade (lê config.txt e roda as etapas em sequência).")
    parser.add_argument("--config", default="config.txt", help="Caminho do config.txt (padrão: config.txt).")
    args = parser.parse_args()

    if not Path(args.config).is_file():
        criar_config_modelo(args.config)
        sys.exit(1)

    if not MODS_DIR.is_dir():
        print(f"[ERRO] Pasta de módulos não encontrada: {MODS_DIR}", file=sys.stderr)
        sys.exit(1)

    config = carregar_config(args.config)
    python_exe = sys.executable

    results_dir = Path(get_str(config, "results_dir", "results"))
    results_dir.mkdir(parents=True, exist_ok=True)
    print(f"[INFO] Resultados serão salvos em: {results_dir}/")

    sintomas = get_list(config, "sintomas")
    farmaco_input = get_str(config, "farmaco_input", "input.txt")  # entrada do usuário, fica na raiz
    farmaco_output = str(results_dir / get_str(config, "farmaco_output", "farmaco.csv"))
    pred_type = get_str(config, "pred_type", "admet")
    batch_size = get_int(config, "batch_size", 5)
    checkpoint = str(results_dir / "deeppk_progress.json")

    ranking_output = str(results_dir / get_str(config, "ranking_output", "ranking.csv"))
    top_por_fonte = get_int(config, "top_por_fonte", 10)
    no_admet = get_bool(config, "no_admet", False)
    pular_revisao = get_bool(config, "pular_revisao_referencias", False)

    caminho_referencia = str(Path(ranking_output).with_suffix("")) + ".referencia.csv"
    caminho_pesos = str(Path(ranking_output).with_suffix("")) + ".pesos.json"

    rodar_farmaco = get_bool(config, "rodar_farmaco", True)
    rodar_ranking = get_bool(config, "rodar_ranking", True)
    rodar_comparacao = get_bool(config, "rodar_comparacao_metodos", False)
    rodar_validacao = get_bool(config, "rodar_validacao_enriquecimento", False)
    rodar_embeddings = get_bool(config, "rodar_embeddings", False)

    print("[INFO] Pipeline de drogabilidade -- configuração carregada de:", args.config)
    print(f"[INFO] Sintoma(s): {sintomas}")

    # ---- Etapa 1: mods/farmaco_completo.py ----
    if rodar_farmaco:
        if not verificar_e_instalar_dependencias(DEPENDENCIAS_OBRIGATORIAS["farmaco"]):
            print("[ERRO] Dependências obrigatórias da etapa farmaco_completo.py não puderam ser instaladas.", file=sys.stderr)
            sys.exit(1)
        comando = [
            python_exe, str(MODS_DIR / "farmaco_completo.py"),
            "--input", farmaco_input,
            "--output", farmaco_output,
            "--pred-type", pred_type,
            "--batch-size", str(batch_size),
            "--checkpoint", checkpoint,
        ]
        if not rodar_etapa("farmaco_completo.py (input.txt -> farmaco.csv)", comando):
            sys.exit(1)
    else:
        print("\n[INFO] 'rodar_farmaco = false' no config.txt -- etapa pulada.")
        if not Path(farmaco_output).is_file():
            print(f"[ERRO] '{farmaco_output}' não existe e a etapa foi pulada. Rode com 'rodar_farmaco = true' ao menos uma vez.", file=sys.stderr)
            sys.exit(1)

    # ---- Etapa 2: mods/ranking_multibase.py ----
    if rodar_ranking:
        if not sintomas:
            print("[ERRO] Nenhum sintoma configurado ('sintomas' vazio no config.txt).", file=sys.stderr)
            sys.exit(1)
        if not verificar_e_instalar_dependencias(DEPENDENCIAS_OBRIGATORIAS["ranking"]):
            print("[ERRO] Dependências obrigatórias da etapa ranking_multibase.py não puderam ser instaladas.", file=sys.stderr)
            sys.exit(1)
        verificar_dependencias_opcionais(DEPENDENCIAS_OPCIONAIS["ranking"])
        comando = [
            python_exe, str(MODS_DIR / "ranking_multibase.py"),
            "--sintoma", *sintomas,
            "--farmaco", farmaco_output,
            "--output", ranking_output,
            "--top-por-fonte", str(top_por_fonte),
        ]
        if no_admet:
            comando.append("--no-admet")
        if pular_revisao:
            comando.append("--pular-revisao-referencias")
        if not rodar_etapa("ranking_multibase.py (sintoma -> ranking.csv)", comando):
            sys.exit(1)
    else:
        print("\n[INFO] 'rodar_ranking = false' no config.txt -- etapa pulada.")

    # ---- Etapa 3 (opcional): mods/comparar_metodos_ranking.py ----
    if rodar_comparacao:
        if not verificar_e_instalar_dependencias(DEPENDENCIAS_OBRIGATORIAS["comparacao"]):
            print("[AVISO] Dependências da etapa comparar_metodos_ranking.py não puderam ser instaladas. Etapa pulada.")
        else:
            comando = [
                python_exe, str(MODS_DIR / "comparar_metodos_ranking.py"),
                "--farmaco", farmaco_output,
                "--referencia", caminho_referencia,
                "--pesos", caminho_pesos,
                "--output", str(results_dir / "comparacao_metodos.csv"),
            ]
            if not rodar_etapa("comparar_metodos_ranking.py (validação vs. TOPSIS/QED/Tanimoto)", comando):
                print("[AVISO] Etapa de comparação de métodos falhou, mas o pipeline principal já tinha concluído -- não é um erro fatal.")

    # ---- Etapa 4 (opcional): mods/validar_enriquecimento.py ----
    if rodar_validacao:
        if not verificar_e_instalar_dependencias(DEPENDENCIAS_OBRIGATORIAS["validacao"]):
            print("[AVISO] Dependências da etapa validar_enriquecimento.py não puderam ser instaladas. Etapa pulada.")
        else:
            n_decoys = get_int(config, "n_decoys", 50)
            comando = [
                python_exe, str(MODS_DIR / "validar_enriquecimento.py"),
                "--referencia", caminho_referencia,
                "--n-decoys", str(n_decoys),
                "--output", str(results_dir / "validacao_enriquecimento.csv"),
            ]
            if not rodar_etapa("validar_enriquecimento.py (validação leave-one-out + decoys)", comando):
                print("[AVISO] Etapa de validação por enriquecimento falhou, mas o pipeline principal já tinha concluído -- não é um erro fatal.")

    # ---- Etapa 5 (opcional): mods/ranking_embeddings.py ----
    if rodar_embeddings:
        print("\n[INFO] 'rodar_embeddings = true' -- verificando dependências pesadas (transformers/torch)...")
        if not verificar_e_instalar_dependencias(DEPENDENCIAS_OBRIGATORIAS["embeddings"]):
            print("[AVISO] Não foi possível instalar transformers/torch automaticamente. Etapa de embeddings pulada -- "
                  "instale manualmente com 'pip install transformers torch' e rode de novo.")
        else:
            modelo_embeddings = get_str(config, "modelo_embeddings", "seyonec/ChemBERTa-zinc-base-v1")
            similaridade_embeddings = get_str(config, "similaridade_embeddings", "centroide")
            peso_embedding = get_float(config, "peso_embedding", 0.5)
            peso_admet_embeddings = get_float(config, "peso_admet_embeddings", 0.5)
            comando = [
                python_exe, str(MODS_DIR / "ranking_embeddings.py"),
                "--farmaco", farmaco_output,
                "--referencia", caminho_referencia,
                "--pesos", caminho_pesos,
                "--output", str(results_dir / "ranking_embeddings.csv"),
                "--modelo", modelo_embeddings,
                "--similaridade", similaridade_embeddings,
                "--peso-embedding", str(peso_embedding),
                "--peso-admet", str(peso_admet_embeddings),
            ]
            if not rodar_etapa("ranking_embeddings.py (ranking alternativo via ChemBERTa)", comando):
                print("[AVISO] Etapa de embeddings falhou, mas o pipeline principal já tinha concluído -- não é um erro fatal.")

    print("\n" + "=" * 70)
    print("PIPELINE CONCLUÍDO")
    print("=" * 70)
    print(f"  Ranking final:        {ranking_output}")
    print(f"  Drogas de referência: {caminho_referencia}")
    print(f"  Pesos calculados:     {caminho_pesos}")


if __name__ == "__main__":
    main()
