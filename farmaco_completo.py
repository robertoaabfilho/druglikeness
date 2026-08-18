#!/usr/bin/env python3
"""
farmaco_completo.py

Une o cálculo de descritores de drogabilidade com a consulta à API ADMET do Deep-PK em um único fluxo:

    1. Lê 'input.txt' (um SMILES por linha, com ou sem nome — ver formatos
       abaixo).
    2. Calcula localmente, via RDKit, os descritores de drogabilidade
       (PM, doadores/aceptores de HB, TPSA, ligações rotacionáveis,
       refratividade molar, átomos pesados) e as avaliações da Regra dos 5
       parcial (sem LogP, que vem do Deep-PK) e da Regra de Veber.
    3. Consulta a API do Deep-PK para os mesmos SMILES (predições ADMET
       completas: Geral/Absorção/Distribuição/Metabolismo/Excreção/
       Toxicidade), em lotes pequenos, com checkpoint e retomada.
    4. Junta os dois conjuntos de resultados pela posição original de cada
       molécula no arquivo de entrada e salva tudo em 'farmaco.csv'.

Formatos aceitos por linha do input.txt (pode misturar entre linhas):
    CC(=O)Oc1ccccc1C(=O)O
    Aspirina,CC(=O)Oc1ccccc1C(=O)O
    Aspirina;CC(=O)Oc1ccccc1C(=O)O
    Aspirina<TAB>CC(=O)Oc1ccccc1C(=O)O

Linhas vazias ou iniciadas com '#' são ignoradas. Se nenhum nome for
fornecido, a molécula é numerada automaticamente (Molecula_1, Molecula_2...).

Uso:
    python farmaco_completo.py --input input.txt --output farmaco.csv --pred-type admet --batch-size 5

Requisitos:
    pip install rdkit requests
"""

import argparse
import csv
import json
import os
import re
import sys
import tempfile
import time
from pathlib import Path

import requests
from rdkit import Chem
from rdkit.Chem import Descriptors

# =========================================================================
# PARTE 1 — druglikeness_rdkit.py (descritores locais via RDKit)
# =========================================================================

def calcular_descritores(smiles: str, nome: str = "") -> dict:
    """
    Calcula os descritores de drogabilidade para uma molécula a partir do SMILES.
    Retorna um dicionário com os valores, ou None nos campos se o SMILES for inválido.
    """
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return {
            "Nome": nome,
            "SMILES": smiles,
            "Valido": False,
            "PesoMolecular": None,
            "DoadoresHB": None,
            "AceptoresHB": None,
            "TPSA": None,
            "LigacoesRotacionaveis": None,
            "RefratividadeMolar": None,
            "AtomosPesados": None,
        }

    return {
        "Nome": nome,
        "SMILES": smiles,
        "Valido": True,
        "PesoMolecular": round(Descriptors.MolWt(mol), 2),
        "DoadoresHB": Descriptors.NumHDonors(mol),
        "AceptoresHB": Descriptors.NumHAcceptors(mol),
        "TPSA": round(Descriptors.TPSA(mol), 2),
        "LigacoesRotacionaveis": Descriptors.NumRotatableBonds(mol),
        "RefratividadeMolar": round(Descriptors.MolMR(mol), 2),
        "AtomosPesados": Descriptors.HeavyAtomCount(mol),
    }


def avaliar_regra_dos_5(descritores: dict) -> dict:
    """
    Avalia a Regra dos 5 de Lipinski (sem o critério de LogP, pois este
    já vem do Deep-PK). Considera até 1 violação como aceitável (padrão Lipinski).
    """
    if not descritores["Valido"]:
        return {"Violacoes_Ro5_parcial": None, "Passa_Ro5_parcial": None}

    violacoes = 0
    if descritores["PesoMolecular"] > 500:
        violacoes += 1
    if descritores["DoadoresHB"] > 5:
        violacoes += 1
    if descritores["AceptoresHB"] > 10:
        violacoes += 1

    return {
        "Violacoes_Ro5_parcial": violacoes,
        "Passa_Ro5_parcial": violacoes <= 1,
    }


def avaliar_regra_de_veber(descritores: dict) -> dict:
    """
    Regra de Veber: boa biodisponibilidade oral quando
    TPSA <= 140 e ligações rotacionáveis <= 10.
    """
    if not descritores["Valido"]:
        return {"Passa_Veber": None}

    passa = descritores["TPSA"] <= 140 and descritores["LigacoesRotacionaveis"] <= 10
    return {"Passa_Veber": passa}


def ler_smiles_do_arquivo(caminho: str) -> list[tuple[str, str]]:
    """
    Lê um arquivo de texto com um SMILES por linha e retorna uma lista
    de tuplas (nome, smiles).

    Formatos aceitos por linha:
        SMILES
        Nome,SMILES
        Nome;SMILES
        Nome<TAB>SMILES

    Linhas vazias ou iniciadas com '#' são ignoradas.
    """
    moleculas = []
    contador = 0

    with open(caminho, "r", encoding="utf-8") as f:
        for linha in f:
            linha = linha.strip()

            if not linha or linha.startswith("#"):
                continue

            # Tenta separar por vírgula, ponto-e-vírgula ou tab, nessa ordem
            partes = None
            for separador in (",", ";", "\t"):
                if separador in linha:
                    partes = linha.split(separador, 1)
                    break

            if partes and len(partes) == 2:
                nome, smiles = partes[0].strip(), partes[1].strip()
            else:
                contador += 1
                nome, smiles = f"Molecula_{contador}", linha

            moleculas.append((nome, smiles))

    return moleculas


def processar_lista_druglikeness(moleculas: list[tuple[str, str]]) -> list[dict]:
    """
    Processa uma lista de tuplas (nome, smiles) e retorna os descritores
    de drogabilidade combinados com as avaliações de regras, na MESMA
    ORDEM da lista de entrada (posição = indice_original usado depois
    para juntar com os resultados do Deep-PK).
    """
    resultados = []
    for nome, smiles in moleculas:
        desc = calcular_descritores(smiles, nome)
        desc.update(avaliar_regra_dos_5(desc))
        desc.update(avaliar_regra_de_veber(desc))
        resultados.append(desc)
    return resultados


CAMPOS_DRUGLIKENESS = [
    "Nome",
    "SMILES",
    "Valido",
    "PesoMolecular",
    "DoadoresHB",
    "AceptoresHB",
    "TPSA",
    "LigacoesRotacionaveis",
    "RefratividadeMolar",
    "AtomosPesados",
    "Violacoes_Ro5_parcial",
    "Passa_Ro5_parcial",
    "Passa_Veber",
]


# =========================================================================
# PARTE 2 — deeppk_query.py (consulta à API Deep-PK)
# =========================================================================

API_URL = "https://biosig.lab.uq.edu.au/deeppk/api/predict"
POLL_INTERVAL_SECONDS = 15      # intervalo entre consultas de status de um job
MAX_WAIT_SECONDS = 30 * 60      # tempo máximo de espera por um job (30 min)
CSV_DELIMITADOR = ";"           # separador usado no CSV de saída final


def barra_progresso(atual: int, total: int, largura: int = 30) -> str:
    """Barra de progresso simples em texto, ex: '[#####-----] 50% (5/10)'."""
    proporcao = 1.0 if total == 0 else min(atual / total, 1.0)
    preenchido = int(largura * proporcao)
    barra = "#" * preenchido + "-" * (largura - preenchido)
    return f"[{barra}] {int(proporcao * 100)}% ({atual}/{total})"


def barra_espera(decorrido_s: int, maximo_s: int, largura: int = 20) -> str:
    """Barra de espera para o tempo de processamento de um job na API."""
    proporcao = min(decorrido_s / maximo_s, 1.0) if maximo_s else 1.0
    preenchido = int(largura * proporcao)
    barra = "#" * preenchido + "-" * (largura - preenchido)
    return f"Aguardando Deep-PK [{barra}] {decorrido_s}s"

# A API do Deep-PK devolve cada parâmetro em 3 chaves com o formato:
#   "[Categoria/Nome do Parâmetro] Predictions"
#   "[Categoria/Nome do Parâmetro] Probability"
#   "[Categoria/Nome do Parâmetro] Interpretation"
REGEX_CHAVE = re.compile(r"^\[(?P<categoria>[^/\]]+)/(?P<propriedade>.+)\]\s+(?P<tipo>Predictions|Probability|Interpretation)$")

# Mapeia o nome de categoria usado pela API para o prefixo de coluna no CSV.
# A ordem deste dicionário define a ordem das colunas no arquivo final.
CATEGORIA_PREFIXOS = {
    "General Properties": "Geral",
    "Absorption": "Absorcao",
    "Distribution": "Distribuicao",
    "Metabolism": "Metabolismo",
    "Excretion": "Excrecao",
    "Toxicity": "Toxicidade",
}
COLUNA_OUTRAS = "Outras_Propriedades"  # rede de segurança para parâmetros não reconhecidos


def carregar_checkpoint(caminho_checkpoint: str) -> int:
    """Retorna o índice (0-based) do próximo SMILES a processar."""
    if not os.path.exists(caminho_checkpoint):
        return 0
    with open(caminho_checkpoint, "r", encoding="utf-8") as f:
        estado = json.load(f)
    return estado.get("proximo_indice", 0)


def salvar_checkpoint(caminho_checkpoint: str, proximo_indice: int) -> None:
    with open(caminho_checkpoint, "w", encoding="utf-8") as f:
        json.dump({"proximo_indice": proximo_indice}, f)


def slugify(texto: str) -> str:
    """Converte um nome de parâmetro (ex: 'Log(D) at pH=7.4') em um nome de coluna seguro."""
    texto = re.sub(r"[^\w]+", "_", texto, flags=re.UNICODE)
    return texto.strip("_")


def limpar_texto(texto) -> str:
    """Remove apenas as tags/entidades HTML específicas que a API às vezes
    inclui ('<br/>', '&nbsp;'), sem mexer em sinais de comparação como
    '<' e '>' que aparecem legitimamente em textos como 'log vp < 4'."""
    texto = str(texto)
    texto = re.sub(r"<br\s*/?>", " ", texto, flags=re.IGNORECASE)
    texto = texto.replace("&nbsp;", " ")
    return re.sub(r"\s+", " ", texto).strip()


def descobrir_colunas_existentes(caminho_output: str) -> list[str] | None:
    """
    Se o CSV intermediário do Deep-PK já existir (retomada), lê o cabeçalho
    e recupera a lista de colunas de parâmetros já usada, para manter a
    mesma estrutura.
    """
    if not os.path.exists(caminho_output) or os.path.getsize(caminho_output) == 0:
        return None
    with open(caminho_output, "r", newline="", encoding="utf-8") as f:
        header = next(csv.reader(f, delimiter=CSV_DELIMITADOR), None)
    if not header:
        return None
    fixas = {"indice_original", "SMILES", COLUNA_OUTRAS}
    return [c for c in header if c not in fixas]


def extrair_colunas(resultados: dict) -> list[str]:
    """
    Descobre dinamicamente, a partir da resposta da API, uma coluna para
    cada parâmetro presente, agrupadas por categoria (na ordem definida em
    CATEGORIA_PREFIXOS) e, dentro de cada categoria, na ordem em que
    aparecem na resposta.
    """
    encontrados: dict[str, dict[str, str]] = {cat: {} for cat in CATEGORIA_PREFIXOS}

    for mol in resultados.values():
        if not isinstance(mol, dict):
            continue
        for chave in mol:
            m = REGEX_CHAVE.match(chave)
            if not m or m.group("tipo") != "Predictions":
                continue
            categoria = m.group("categoria").strip()
            propriedade = m.group("propriedade").strip()
            if categoria in encontrados:
                encontrados[categoria].setdefault(slugify(propriedade), propriedade)

    colunas = []
    for categoria, prefixo in CATEGORIA_PREFIXOS.items():
        for slug in encontrados[categoria]:
            colunas.append(f"{prefixo}_{slug}")
    return colunas


def construir_colunas_deeppk(colunas_parametros: list[str]) -> list[str]:
    return ["indice_original", "SMILES"] + colunas_parametros + [COLUNA_OUTRAS]


MSG_MAX = 150  # tamanho máximo de trecho de resposta bruta incluído nas mensagens de erro


def _resumo(obj, limite: int = MSG_MAX) -> str:
    """Formata um trecho curto de um objeto/resposta para usar em mensagens de erro."""
    texto = str(obj).replace("\n", " ")
    return texto if len(texto) <= limite else texto[:limite] + "..."


def submeter_job(smiles_lote: list[str], pred_type: str, email: str | None) -> str:
    """
    Envia um pequeno lote de SMILES (arquivo temporário) para a API do Deep-PK
    e retorna o job_id.
    """
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".txt", delete=False, encoding="utf-8"
    ) as tmp:
        tmp.write("\n".join(smiles_lote))
        caminho_tmp = tmp.name

    try:
        with open(caminho_tmp, "rb") as f:
            files = {"smiles_file": f}
            data = {"pred_type": pred_type}
            if email:
                data["email"] = email
            resposta = requests.post(API_URL, files=files, data=data, timeout=120)
    finally:
        os.remove(caminho_tmp)

    resposta.raise_for_status()

    try:
        payload = _decodificar_json(resposta.text)
    except (json.JSONDecodeError, TypeError):
        raise RuntimeError(f"resposta inválida ao submeter job: {_resumo(resposta.text)}")

    job_id = payload.get("job_id") if isinstance(payload, dict) else None
    if not job_id:
        raise RuntimeError(f"resposta inesperada ao submeter job: {_resumo(payload)}")

    return job_id


def _decodificar_json(texto: str):
    """
    Decodifica a resposta da API. O Deep-PK às vezes retorna o resultado
    como uma STRING JSON aninhada (JSON dentro de JSON) em vez de um objeto
    direto, então aqui fazemos json.loads() uma segunda vez quando necessário.
    """
    dado = json.loads(texto)
    if isinstance(dado, str):
        dado = json.loads(dado)
    return dado


def aguardar_resultado(job_id: str, mostrar_espera: bool = True) -> dict:
    """
    Consulta periodicamente o status do job até que o processamento termine,
    retornando o dicionário de resultados. Enquanto o job está rodando,
    mostra uma barra de espera de uma linha só (sobrescrita a cada consulta).
    """
    inicio = time.time()

    try:
        while True:
            resposta = requests.get(API_URL, data={"job_id": job_id}, timeout=60)
            resposta.raise_for_status()

            try:
                payload = _decodificar_json(resposta.text)
            except (json.JSONDecodeError, TypeError):
                raise RuntimeError(f"resposta inválida para o job {job_id}: {_resumo(resposta.text)}")

            if isinstance(payload, dict) and payload.get("status") == "running":
                decorrido = int(time.time() - inicio)
                if mostrar_espera:
                    print(f"\r[INFO]   {barra_espera(decorrido, MAX_WAIT_SECONDS)}", end="", flush=True)
                if decorrido > MAX_WAIT_SECONDS:
                    raise TimeoutError(f"tempo máximo de espera ({MAX_WAIT_SECONDS}s) excedido")
                time.sleep(POLL_INTERVAL_SECONDS)
                continue

            if not isinstance(payload, dict):
                raise RuntimeError(f"resposta inesperada (esperava dict, veio {type(payload).__name__}): {_resumo(payload)}")

            # A API sinaliza falha do job de formas variadas (status explícito
            # de erro, chave "error"/"message", ou payload sem nenhuma
            # molécula processada). Antes isso era tratado como resultado
            # válido e vazio — agora falha alto, sem gerar um CSV "completo"
            # e silenciosamente sem dados.
            status = str(payload.get("status", "")).lower()
            if status in ("error", "failed", "failure"):
                raise RuntimeError(f"status de erro ({payload.get('status')!r}): {_resumo(payload)}")
            if ("error" in payload or "message" in payload) and not any(k.isdigit() for k in payload):
                raise RuntimeError(f"erro da API: {_resumo(payload)}")
            if not any(k.isdigit() for k in payload):
                raise RuntimeError(f"nenhuma molécula retornada: {_resumo(payload)}")

            return payload
    finally:
        if mostrar_espera:
            print("\r" + " " * 60 + "\r", end="", flush=True)  # limpa a linha da barra de espera


def resumir_molecula(dados_molecula: dict, colunas_parametros: list[str]) -> dict:
    """
    Monta a linha de uma molécula com uma coluna para cada parâmetro
    (formato "predição (interpretação)" quando há interpretação disponível,
    senão só a predição). Parâmetros fora do conjunto de colunas conhecido
    caem em 'Outras_Propriedades'.
    """
    valores: dict[str, str] = {c: "" for c in colunas_parametros}
    outras: list[str] = []

    for chave, valor in dados_molecula.items():
        m = REGEX_CHAVE.match(chave)
        if not m or m.group("tipo") != "Predictions":
            continue  # Probability/Interpretation são lidas via get() abaixo, ancoradas em Predictions

        categoria = m.group("categoria").strip()
        propriedade = m.group("propriedade").strip()

        chave_interp = f"[{categoria}/{propriedade}] Interpretation"
        interpretacao = dados_molecula.get(chave_interp, "-")
        tem_interpretacao = interpretacao is not None and str(interpretacao).strip() not in ("-", "None", "")
        if tem_interpretacao:
            interpretacao = limpar_texto(interpretacao)

        texto = f"{valor} ({interpretacao})" if tem_interpretacao else f"{valor}"

        prefixo = CATEGORIA_PREFIXOS.get(categoria)
        nome_coluna = f"{prefixo}_{slugify(propriedade)}" if prefixo else None

        if nome_coluna and nome_coluna in valores:
            valores[nome_coluna] = texto
        else:
            outras.append(f"{categoria}/{propriedade}: {texto}")

    valores[COLUNA_OUTRAS] = " | ".join(outras)
    return valores


def resultados_para_linhas(resultados: dict, lote: list[tuple[int, str]], colunas_parametros: list[str]) -> list[dict]:
    """
    Converte o dicionário de resultados de um lote em uma lista de linhas.

    'lote' é a lista de (indice_original, smiles) que foi de fato enviada
    à API nesse job — usada para traduzir o índice local devolvido pela
    API (0, 1, 2...) de volta para a posição original no input.txt (o que
    é essencial para juntar depois com os descritores de drogabilidade,
    já que moléculas com SMILES inválido são puladas antes de chegar aqui).
    """
    linhas = []
    indices_moleculas = sorted((k for k in resultados if k.isdigit()), key=lambda x: int(x))
    for indice_local in indices_moleculas:
        dados_molecula = resultados[indice_local]
        if not isinstance(dados_molecula, dict):
            continue

        pos = int(indice_local)
        if pos >= len(lote):
            continue  # resposta da API com mais entradas que o lote enviado; ignora sobra
        indice_original_real = lote[pos][0]

        linha = {
            "indice_original": indice_original_real,
            "SMILES": dados_molecula.get("SMILES", ""),
        }
        linha.update(resumir_molecula(dados_molecula, colunas_parametros))
        linhas.append(linha)
    return linhas


def anexar_csv(linhas: list[dict], caminho_output: str, colunas_finais: list[str]) -> None:
    """
    Anexa linhas ao CSV intermediário do Deep-PK. Escreve o cabeçalho apenas
    se o arquivo ainda não existir (ou estiver vazio) — permite retomar sem
    duplicar header.
    """
    if not linhas:
        return

    arquivo_existe = os.path.exists(caminho_output) and os.path.getsize(caminho_output) > 0

    modo = "a" if arquivo_existe else "w"
    with open(caminho_output, modo, newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=colunas_finais, delimiter=CSV_DELIMITADOR)
        if not arquivo_existe:
            writer.writeheader()
        writer.writerows(linhas)


def _processar_lote_com_bissecao(
    lote: list[tuple[int, str]],
    pred_type: str,
    email: str | None,
    estado: dict,
    falhas: list[tuple[int, str, str]],
    caminho_raw: str,
    descritores_druglikeness: list[dict],
    caminho_saida: str,
    profundidade: int = 0,
) -> list[dict]:
    """
    Tenta consultar o Deep-PK para 'lote' inteiro de uma vez. Se o job
    falhar (erro da API, timeout, resposta sem moléculas etc.), divide o
    lote ao meio e tenta cada metade recursivamente — isso isola
    automaticamente qual(is) SMILES específico(s) está(ão) quebrando o job,
    sem descartar o restante do lote que seria perfeitamente válido.

    Quando um SMILES falha mesmo sozinho (lote de tamanho 1), ele é
    registrado em 'falhas' e fica sem dados ADMET no resultado final — o
    restante do arquivo continua normalmente.

    A cada sucesso (lote inteiro OU qualquer sublote resultante de uma
    divisão), o resultado é imediatamente anexado ao CSV intermediário do
    Deep-PK ('caminho_raw') e 'caminho_saida' (farmaco.csv) é regravado do
    zero já com esse progresso — assim o arquivo final vai crescendo aos
    poucos em vez de só aparecer no fim de tudo.

    A API do Deep-PK rejeita arquivos de upload em lote com apenas 1 SMILES
    (erro "No valid molecules were provided..."), então quando o lote cai
    para tamanho 1 aqui, a molécula é enviada duplicada (2x) só para
    satisfazer esse mínimo do endpoint — apenas o primeiro resultado é
    usado. Isso permite diferenciar um SMILES genuinamente rejeitado pelo
    Deep-PK de um falso positivo causado por esse detalhe de formato.

    'estado' é um dict mutável com a chave "colunas" (list[str] | None),
    usado para descobrir as colunas de parâmetros ADMET na primeira
    resposta bem-sucedida e reaproveitá-las nas chamadas seguintes.
    """
    smiles_lote = [smiles for _, smiles in lote]
    indent = "  " * (profundidade + 1)

    # Contorno do mínimo de 2 SMILES exigido pelo endpoint de lote da API.
    smiles_enviados = smiles_lote * 2 if len(smiles_lote) == 1 else smiles_lote

    try:
        job_id = submeter_job(smiles_enviados, pred_type, email)
        resultados = aguardar_resultado(job_id)

        if not any(k.isdigit() for k in resultados):
            raise RuntimeError(f"sem moléculas na resposta: {_resumo(resultados)}")

        if estado["colunas"] is None:
            estado["colunas"] = extrair_colunas(resultados)
            print(f"[INFO]{indent} {len(estado['colunas'])} colunas de par\u00e2metros ADMET detectadas.")
            if not estado["colunas"]:
                print(
                    f"[AVISO]{indent} A API respondeu mas nenhum parâmetro ADMET reconhecido foi "
                    f"encontrado nas chaves esperadas. Os dados brutos irão para '{COLUNA_OUTRAS}'."
                )

        linhas = resultados_para_linhas(resultados, lote, estado["colunas"])
        if not linhas:
            raise RuntimeError("nenhuma linha associada ao lote enviado")
        if profundidade > 0:
            print(f"[INFO]{indent} Sublote de {len(lote)} SMILES OK após divisão.")

        colunas_finais = construir_colunas_deeppk(estado["colunas"] or [])
        anexar_csv(linhas, caminho_raw, colunas_finais)

        linhas_deeppk_por_indice = ler_csv_deeppk(caminho_raw)
        combinar_e_salvar(
            descritores_druglikeness=descritores_druglikeness,
            linhas_deeppk_por_indice=linhas_deeppk_por_indice,
            colunas_parametros_deeppk=estado["colunas"] or [],
            caminho_saida=caminho_saida,
        )
        print(f"[INFO]{indent} '{caminho_saida}' atualizado ({len(linhas)} molécula(s) deste sublote).")

        return linhas

    except Exception as erro:
        erro_curto = _resumo(erro, 100)
        if len(lote) == 1:
            idx, smiles = lote[0]
            print(
                f"[AVISO]{indent} índice {idx} falhou sozinho no Deep-PK ({erro_curto}). "
                f"Sem dados ADMET. SMILES: {_resumo(smiles, 60)}"
            )
            falhas.append((idx, smiles, str(erro)))
            return []

        meio = len(lote) // 2
        print(f"[AVISO]{indent} lote de {len(lote)} falhou ({erro_curto}) — dividindo...")
        linhas_esq = _processar_lote_com_bissecao(
            lote[:meio], pred_type, email, estado, falhas,
            caminho_raw, descritores_druglikeness, caminho_saida, profundidade + 1,
        )
        linhas_dir = _processar_lote_com_bissecao(
            lote[meio:], pred_type, email, estado, falhas,
            caminho_raw, descritores_druglikeness, caminho_saida, profundidade + 1,
        )
        return linhas_esq + linhas_dir


def consultar_deeppk(
    pares_validos: list[tuple[int, str]],
    caminho_raw: str,
    caminho_checkpoint: str,
    pred_type: str,
    batch_size: int,
    email: str | None,
    descritores_druglikeness: list[dict],
    caminho_saida: str,
) -> tuple[list[str], str]:
    """
    Consulta a API do Deep-PK para os pares (indice_original, smiles) já
    filtrados (apenas SMILES válidos — ver main()), em lotes, salvando
    incrementalmente em 'caminho_raw' (CSV intermediário) com checkpoint em
    'caminho_checkpoint'. A cada lote/sublote bem-sucedido, 'caminho_saida'
    (farmaco.csv) também é regravado com o progresso combinado até aquele
    ponto. Retorna (colunas_parametros, caminho_raw).

    Se um lote falhar como um todo, as moléculas dentro dele são isoladas
    por bissecção (ver _processar_lote_com_bissecao) — só quem realmente
    quebra o Deep-PK sozinho fica sem dados ADMET; o resto do arquivo é
    processado normalmente.
    """
    total = len(pares_validos)

    if total == 0:
        print("[AVISO] Nenhum SMILES válido para consultar no Deep-PK.")
        return descobrir_colunas_existentes(caminho_raw) or [], caminho_raw

    proximo_indice = carregar_checkpoint(caminho_checkpoint)
    if proximo_indice > 0:
        print(f"[INFO] Checkpoint do Deep-PK encontrado. Retomando a partir do SMILES válido nº {proximo_indice + 1} de {total}.")

    colunas_parametros = descobrir_colunas_existentes(caminho_raw)

    if proximo_indice >= total:
        print("[INFO] Todos os SMILES já foram consultados no Deep-PK anteriormente. Nada a fazer.")
        if colunas_parametros is None:
            colunas_parametros = []
        return colunas_parametros, caminho_raw

    estado = {"colunas": colunas_parametros}
    falhas: list[tuple[int, str, str]] = []

    indice_atual = proximo_indice
    while indice_atual < total:
        lote = pares_validos[indice_atual: indice_atual + batch_size]
        n_lote = len(lote)

        _processar_lote_com_bissecao(
            lote, pred_type, email, estado, falhas,
            caminho_raw, descritores_druglikeness, caminho_saida,
        )

        indice_atual += n_lote
        salvar_checkpoint(caminho_checkpoint, indice_atual)

        print(f"[INFO] Deep-PK {barra_progresso(indice_atual, total)}")

    if falhas:
        print(f"[AVISO] {len(falhas)} SMILES sem dados ADMET (falharam sozinhos no Deep-PK):")
        for idx, smiles, erro in falhas:
            print(f"  - índice {idx}: {_resumo(smiles, 60)} — {_resumo(erro, 100)}")

    print(f"[INFO] Consulta ao Deep-PK concluída para os {total} SMILES válidos.")
    if os.path.exists(caminho_checkpoint):
        os.remove(caminho_checkpoint)

    return estado["colunas"] or [], caminho_raw


def ler_csv_deeppk(caminho_raw: str) -> dict[int, dict]:
    """Lê o CSV intermediário do Deep-PK e retorna um dict indice_original -> linha."""
    linhas_por_indice: dict[int, dict] = {}
    if not os.path.exists(caminho_raw):
        return linhas_por_indice
    with open(caminho_raw, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=CSV_DELIMITADOR)
        for linha in reader:
            try:
                idx = int(linha["indice_original"])
            except (KeyError, ValueError, TypeError):
                continue
            linhas_por_indice[idx] = linha
    return linhas_por_indice


# =========================================================================
# PARTE 3 — combinação e salvamento final em farmaco.csv
# =========================================================================

def combinar_e_salvar(
    descritores_druglikeness: list[dict],
    linhas_deeppk_por_indice: dict[int, dict],
    colunas_parametros_deeppk: list[str],
    caminho_saida: str,
) -> None:
    colunas_deeppk_sem_fixas = colunas_parametros_deeppk + [COLUNA_OUTRAS]
    colunas_finais = CAMPOS_DRUGLIKENESS + colunas_deeppk_sem_fixas

    linhas_finais = []
    for idx, desc in enumerate(descritores_druglikeness):
        linha = {campo: desc.get(campo) for campo in CAMPOS_DRUGLIKENESS}

        linha_deeppk = linhas_deeppk_por_indice.get(idx)
        if linha_deeppk:
            for coluna in colunas_deeppk_sem_fixas:
                linha[coluna] = linha_deeppk.get(coluna, "")
        else:
            for coluna in colunas_deeppk_sem_fixas:
                linha[coluna] = ""

        linhas_finais.append(linha)

    with open(caminho_saida, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=colunas_finais, delimiter=CSV_DELIMITADOR)
        writer.writeheader()
        writer.writerows(linhas_finais)

    print(f"[INFO] Arquivo final salvo em: {caminho_saida}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Calcula descritores de drogabilidade (RDKit) e consulta a API "
            "ADMET do Deep-PK para os SMILES de 'input.txt', salvando tudo "
            "combinado em 'farmaco.csv'."
        )
    )
    parser.add_argument("--input", default="input.txt", help="Arquivo de entrada com os SMILES (um por linha).")
    parser.add_argument("--output", default="farmaco.csv", help="Arquivo CSV final de saída.")
    parser.add_argument(
        "--raw-deeppk",
        default=None,
        help="CSV intermediário onde o progresso do Deep-PK é salvo (padrão: <output>.deeppk_raw.csv).",
    )
    parser.add_argument(
        "--checkpoint",
        default="deeppk_progress.json",
        help="Arquivo de checkpoint da consulta ao Deep-PK, para retomada.",
    )
    parser.add_argument(
        "--pred-type",
        default="admet",
        choices=["absorption", "distribution", "metabolism", "excretion", "admet"],
        help="Tipo de predição do Deep-PK a ser realizada (padrão: admet, traz todas as categorias).",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=5,
        help="Quantos SMILES enviar por job ao Deep-PK / salvar por vez (padrão: 5).",
    )
    parser.add_argument("--email", default=None, help="E-mail opcional para notificação ao final do job do Deep-PK.")
    args = parser.parse_args()

    if not Path(args.input).is_file():
        print(f"Arquivo não encontrado: {args.input}")
        sys.exit(1)

    caminho_raw = args.raw_deeppk or (str(Path(args.output).with_suffix("")) + ".deeppk_raw.csv")

    try:
        # ---- Etapa 1: druglikeness (RDKit, local) ----
        moleculas = ler_smiles_do_arquivo(args.input)
        if not moleculas:
            print(f"Nenhum SMILES encontrado em: {args.input}")
            sys.exit(1)

        print(f"[INFO] {len(moleculas)} SMILES lidos de '{args.input}'.")

        descritores_druglikeness = processar_lista_druglikeness(moleculas)
        n_invalidos = sum(1 for r in descritores_druglikeness if not r["Valido"])
        print(f"[INFO] Descritores de drogabilidade calculados ({len(descritores_druglikeness)} moléculas).")
        if n_invalidos:
            print(f"[AVISO] {n_invalidos} SMILES inválido(s) para RDKit encontrado(s).")

        # Salva o farmaco.csv já com os descritores locais, mesmo antes de
        # começar a consultar o Deep-PK — assim o arquivo existe desde já e
        # vai sendo atualizado lote a lote a partir daqui.
        combinar_e_salvar(
            descritores_druglikeness=descritores_druglikeness,
            linhas_deeppk_por_indice=ler_csv_deeppk(caminho_raw),
            colunas_parametros_deeppk=descobrir_colunas_existentes(caminho_raw) or [],
            caminho_saida=args.output,
        )

        # ---- Etapa 2: consulta ADMET ao Deep-PK ----
        # Só moléculas com SMILES válido (segundo o RDKit) são enviadas à
        # API — um único SMILES inválido dentro de um lote costuma derrubar
        # o job inteiro no lado do Deep-PK, o que antes zerava silenciosamente
        # todas as colunas ADMET do arquivo inteiro.
        pares_validos = [
            (idx, smiles)
            for idx, (desc, (_, smiles)) in enumerate(zip(descritores_druglikeness, moleculas))
            if desc["Valido"]
        ]
        n_pulados = len(moleculas) - len(pares_validos)
        if n_pulados:
            print(f"[AVISO] {n_pulados} molécula(s) com SMILES inválido não serão enviadas ao Deep-PK.")

        # A partir daqui, 'args.output' (farmaco.csv) é regravado
        # automaticamente a cada lote/sublote bem-sucedido dentro de
        # consultar_deeppk — não é preciso esperar o fim de tudo para ver
        # o progresso no arquivo.
        colunas_parametros_deeppk, caminho_raw = consultar_deeppk(
            pares_validos=pares_validos,
            caminho_raw=caminho_raw,
            caminho_checkpoint=args.checkpoint,
            pred_type=args.pred_type,
            batch_size=args.batch_size,
            email=args.email,
            descritores_druglikeness=descritores_druglikeness,
            caminho_saida=args.output,
        )

        linhas_deeppk_por_indice = ler_csv_deeppk(caminho_raw)

        if pares_validos and not colunas_parametros_deeppk:
            print(
                "[AVISO] Nenhuma coluna ADMET foi obtida do Deep-PK. O 'farmaco.csv' final "
                "só terá os descritores de drogabilidade locais."
            )

        # ---- Etapa 3: combinar e salvar em farmaco.csv (garantia final,
        # cobre também o caso de 0 SMILES válidos ou retomada já concluída) ----
        combinar_e_salvar(
            descritores_druglikeness=descritores_druglikeness,
            linhas_deeppk_por_indice=linhas_deeppk_por_indice,
            colunas_parametros_deeppk=colunas_parametros_deeppk,
            caminho_saida=args.output,
        )

        print(f"[INFO] Concluído! Resultados combinados de {len(moleculas)} moléculas salvos em '{args.output}'.")

    except KeyboardInterrupt:
        print("\n[AVISO] Interrompido pelo usuário. Progresso do Deep-PK salvo — rode o script novamente para retomar.")
        sys.exit(1)
    except Exception as erro:
        print(f"[ERRO] {erro}", file=sys.stderr)
        print("[AVISO] Progresso do Deep-PK (se houver) foi salvo — rode o script novamente para retomar.")
        sys.exit(1)


if __name__ == "__main__":
    main()
