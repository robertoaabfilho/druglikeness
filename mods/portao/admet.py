# -*- coding: utf-8 -*-
"""
admet.py

ADMET pelo ADMET-AI, LOCAL: nada toca a rede.

Por que local
-------------
Calcular ADMET por API externa custaria ~21.800 requisicoes sequenciais para
a base inteira, em lotes de 5 SMILES -- e, pior, torna indistinguiveis falha
de rede e molecula recusada pelo servico, o que faz um erro de conexao virar
uma propriedade ausente sem aviso. O ADMET-AI resolve isso: e MIT, os dez
modelos Chemprop vem DENTRO do pacote, e a execucao e reproduzivel porque nao
depende do estado de um servidor de terceiro.

num_workers: a medicao que decide o padrao
------------------------------------------
O `ADMETModel` escolhe `num_workers=8` quando ha GPU. Medido nesta maquina,
1.600 moleculas:

    num_workers=8 (padrao do ADMET-AI) ... 384 s
    num_workers=0 ......................... 8,6 s

45x. Sao dez modelos, cada um pagando o nascimento de 8 processos de
carregamento para um lote que a GPU resolve em milissegundos. Este modulo passa
`num_workers=0` por padrao; `--workers` muda.

CPU por padrao, e a medicao que decide
--------------------------------------
A inferencia na GPU NAO e deterministica: duas execucoes da mesma molecula no
mesmo processo diferem em 36 de 312 colunas, na ultima casa (ordem de reducao
nos kernels). A diferenca e de ~1e-8 relativo nos valores, mas o PERCENTIL a
amplifica para ~0,035, porque um percentil e uma posicao entre os ~2.850
aprovados do DrugBank e um fio de diferenca inverte uma comparacao. Na CPU: 0
de 312, bit-exato.

Isso quebraria o modo de repeticao desta etapa, que existe para reproduzir a
saida. Medido o preco (200 moleculas de `input.exemplo.tsv`):

    GPU ... 20,5 ms/mol -> base de 108.935 em ~37 min
    CPU ... 36,4 ms/mol -> base de 108.935 em ~66 min

1,8x, nao 10x. Por 29 minutos numa execucao que precisa ser citavel, a CPU e o
padrao; `--gpu` troca reprodutibilidade por tempo, e o manifesto registra qual
dos dois rodou -- pelo dispositivo REAL, nao pelo pedido.

Degradacao
----------
Sem o pacote `admet_ai` a etapa RODA: `disponivel()` e falso, as colunas de
ADMET saem vazias e o criterio de toxicidade sai `indeterminado`. Ausencia
nunca vira aprovacao nem reprovacao.
"""

import importlib.util
import logging
import os
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


def carregar(num_workers: int = 0, usar_gpu: bool = False):
    """Modelo carregado uma vez por processo (0,1 s de carga, 6 s de import).

    `ADMETModel` nao aceita dispositivo, e o Lightning escolhe "auto" -- entao
    a unica forma de exigir CPU e esconder a GPU ANTES de o torch ser
    importado. Se o torch ja estiver inicializado com CUDA neste processo, isto
    nao tem efeito: por isso quem grava o manifesto usa `dispositivo()`, que
    responde o que de fato esta valendo."""
    global _modelo
    if _modelo is None:
        if not usar_gpu:
            os.environ["CUDA_VISIBLE_DEVICES"] = ""
        _silenciar()
        from admet_ai import ADMETModel
        _modelo = ADMETModel(num_workers=num_workers)
    return _modelo


def dispositivo() -> str:
    """"gpu" ou "cpu": o que a inferencia VAI usar, nao o que foi pedido.

    Um manifesto que registre o pedido mentiria quando o pedido nao pegou --
    e um manifesto que mente sobre o dispositivo promete uma repeticao
    bit-exata que a GPU nao entrega."""
    try:
        import torch
    except ImportError:
        return "indeterminado"
    return "gpu" if torch.cuda.is_available() else "cpu"


def colunas(num_workers: int = 0, usar_gpu: bool = False) -> list[str]:
    """Os nomes das 104 colunas do ADMET-AI, descobertos da instalacao e nao
    fixados no codigo: fixa-los faria o repositorio mentir sobre a versao.
    Custa uma previsao de uma molecula trivial."""
    modelo = carregar(num_workers, usar_gpu)
    return list(modelo.predict("C").keys())


def prever(smiles: list[str], *, num_workers: int = 0, usar_gpu: bool = False,
           lote: int = LOTE) -> dict[str, dict]:
    """{smiles: {coluna: valor}} para os SMILES que o ADMET-AI aceitou.

    SMILES ausente do resultado e um SMILES que a ferramenta descartou (ela
    avisa "N invalid molecules removed" e os remove da tabela). Quem chama
    trata a ausencia como `indeterminado`; aqui nao se inventa valor."""
    modelo = carregar(num_workers, usar_gpu)
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
