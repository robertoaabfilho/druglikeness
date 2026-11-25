# -*- coding: utf-8 -*-
"""
admet.py

ADMET pelo ADMET-AI, LOCAL. Substitui a consulta a API do Deep-PK.

Por que trocar
--------------
O caminho antigo (`mods/farmaco_completo.py`) manda lotes de 5 SMILES a uma
API externa: ~21.800 requisicoes sequenciais para a base inteira, com tres
defeitos no mesmo trecho -- falha de rede rotulada como "rejeitado pelo
Deep-PK" (linhas 583-585), checkpoint avancado de todo jeito (607) e apagado
no fim (613). O ADMET-AI resolve os tres de uma vez: e MIT, os dez modelos
Chemprop vem DENTRO do pacote e nada toca a rede.

num_workers: a medicao que decide o padrao
------------------------------------------
O `ADMETModel` escolhe `num_workers=8` quando ha GPU. Medido nesta maquina,
1.600 moleculas:

    num_workers=8 (padrao do ADMET-AI) ... 384 s
    num_workers=0 ......................... 8,6 s

45x. Sao dez modelos, cada um pagando o nascimento de 8 processos de
carregamento para um lote que a GPU resolve em milissegundos. Este modulo passa
`num_workers=0` por padrao; `--workers` muda.

Degradacao
----------
Sem o pacote `admet_ai` a etapa RODA: `disponivel()` e falso, as colunas de
ADMET saem vazias e o criterio de toxicidade sai `indeterminado`. Ausencia
nunca vira aprovacao nem reprovacao.
"""

import importlib.util
import logging
import warnings

LOTE = 2048

# Endpoints de toxicidade que o criterio `toxicidade` olha. Nomes do ADMET-AI,
# nao nossos: a coluna tem de ser citavel no artigo que cita a ferramenta.
ENDPOINTS_TOXICIDADE = ("hERG", "AMES", "DILI", "ClinTox", "Carcinogens_Lagunin")
SUFIXO_PERCENTIL = "_drugbank_approved_percentile"

_modelo = None


def disponivel() -> bool:
    return importlib.util.find_spec("admet_ai") is not None


def versao() -> str:
    from importlib.metadata import PackageNotFoundError, version
    try:
        return version("admet_ai")
    except PackageNotFoundError:
        return "desconhecida"


def _silenciar():
    """O Lightning imprime um bloco de progresso e dicas por modelo e por lote
    -- dez vezes por lote. O ruido esconderia o progresso da propria etapa.
    So os registradores dele; nada global."""
    for nome in ("lightning", "lightning.pytorch", "pytorch_lightning",
                 "lightning.pytorch.utilities.rank_zero"):
        logging.getLogger(nome).setLevel(logging.ERROR)
    warnings.filterwarnings("ignore", module="lightning.*")
    warnings.filterwarnings("ignore", module="chemprop.*")


def carregar(num_workers: int = 0):
    """Modelo carregado uma vez por processo (0,1 s de carga, 6 s de import)."""
    global _modelo
    if _modelo is None:
        _silenciar()
        from admet_ai import ADMETModel
        _modelo = ADMETModel(num_workers=num_workers)
    return _modelo


def colunas(num_workers: int = 0) -> list[str]:
    """Os nomes das 104 colunas do ADMET-AI, descobertos da instalacao e nao
    fixados no codigo: fixa-los faria o repositorio mentir sobre a versao.
    Custa uma previsao de uma molecula trivial."""
    modelo = carregar(num_workers)
    return list(modelo.predict("C").keys())


def prever(smiles: list[str], *, num_workers: int = 0,
           lote: int = LOTE) -> dict[str, dict]:
    """{smiles: {coluna: valor}} para os SMILES que o ADMET-AI aceitou.

    SMILES ausente do resultado e um SMILES que a ferramenta descartou (ela
    avisa "N invalid molecules removed" e os remove da tabela). Quem chama
    trata a ausencia como `indeterminado`; aqui nao se inventa valor."""
    modelo = carregar(num_workers)
    saida: dict[str, dict] = {}
    unicos = list(dict.fromkeys(s for s in smiles if s))
    for inicio in range(0, len(unicos), lote):
        bloco = unicos[inicio:inicio + lote]
        tabela = modelo.predict(bloco)
        if isinstance(tabela, dict):   # um SMILES so: a API devolve dict
            saida[bloco[0]] = tabela
            continue
        for chave, linha in tabela.to_dict(orient="index").items():
            saida[str(chave)] = linha
    return saida


def veredito_de_toxicidade(linha: dict | None, percentil: float):
    """(veredito, valor) do criterio `toxicidade`.

    Olha o PERCENTIL contra os aprovados do DrugBank, nao a probabilidade crua:
    0,67 de DILI nao quer dizer nada sozinho -- a aspirina tem 0,67 e esta
    aprovada. Percentil 90 quer dizer "pior que 90% dos farmacos aprovados",
    que e uma afirmacao comparavel entre endpoints de escalas diferentes.

    Sem linha (ADMET-AI ausente ou SMILES descartado) o veredito e
    `indeterminado`: nao se aprova nem se reprova o que nao foi medido."""
    if not linha:
        return "indeterminado", "sem previsao de ADMET"
    piores = []
    for endpoint in ENDPOINTS_TOXICIDADE:
        valor = linha.get(endpoint + SUFIXO_PERCENTIL)
        if valor is None or valor != valor:
            continue
        if valor >= percentil:
            piores.append(f"{endpoint}={valor:.0f}")
    if not piores:
        return "passa", f"nenhum endpoint acima do percentil {percentil:g}"
    return "marca", ";".join(piores)
