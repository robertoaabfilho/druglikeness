# -*- coding: utf-8 -*-
"""Fixtures compartilhadas. Tarefas seguintes estendem este arquivo."""
import pytest

N_ALVOS_FALSO = 6


def pytest_sessionstart(session):
    """Sem torch a suite NAO e parcial, e inutilizavel: o painel, o CLI e a
    evidencia prevista dependem dele, e o resultado seria dezenas de falhas
    espurias (ou um INTERNALERROR no import) ou, pior, um 'passed' que omite
    os testes de regressao, como o do guarda de largura do checkpoint. Falhar
    cedo, com o remedio dito, e o unico resultado que nao engana."""
    import importlib.util
    if importlib.util.find_spec("torch") is None:
        pytest.exit(
            "torch NAO esta instalado: a suite completa (incluindo os testes da "
            "evidencia prevista e o guarda de largura do checkpoint) nao pode "
            "rodar. Instale com: pip install -r requirements-opcionais.txt",
            returncode=2)


@pytest.fixture(autouse=True)
def sem_rede(monkeypatch):
    """Nenhum teste fala com a rede: conectar um socket de internet falha o
    teste. Quem esquecer de injetar `obter_json` descobre aqui, nao na maquina
    de outra pessoa. (Sockets UNIX e subprocessos ficam de fora.)"""
    import socket
    conectar = socket.socket.connect

    def bloqueado(self, endereco, *a, **k):
        if self.family in (socket.AF_INET, socket.AF_INET6):
            raise RuntimeError(f"teste tentou usar a rede: {endereco}")
        return conectar(self, endereco, *a, **k)
    monkeypatch.setattr(socket.socket, "connect", bloqueado)


@pytest.fixture
def modelo_falso(tmp_path):
    """state_dict aleatorio com a arquitetura do PP4CKD (4096 -> 1000 -> 500 ->
    N), so que com poucos alvos. Exercita o mesmo caminho de carga
    (torch.load com weights_only=True) sem os 785 MB dos pesos reais; os
    valores das probabilidades nao significam nada, so seu formato e faixa."""
    torch = pytest.importorskip("torch")
    from mods.prever_alvos import construir_rede, FINGERPRINTS
    torch.manual_seed(0)
    rede = construir_rede(FINGERPRINTS["ECFP4"]["bits"], N_ALVOS_FALSO)
    caminho = tmp_path / "ecfp4_falso.pt"
    torch.save(rede.state_dict(), caminho)
    return str(caminho)


@pytest.fixture
def rotulos_falsos():
    """Rotulos no formato de carregar_rotulos, indexados como o modelo falso."""
    return [
        {"chembl_id": f"CHEMBL{100 + i}", "nome": f"alvo falso {i}",
         "organismo": "Homo sapiens" if i % 2 == 0 else "Rattus norvegicus",
         "tipo": "SINGLE PROTEIN", "n_ativos": 4000 if i == 0 else 8 + i}
        for i in range(N_ALVOS_FALSO)
    ]


BASE_CHEMBL = "https://www.ebi.ac.uk/chembl/api/data/"
BASE_REACTOME = "https://reactome.org/ContentService/"
BASE_OPENFDA = "https://api.fda.gov/"
URL_OPEN_TARGETS = "https://api.platform.opentargets.org/api/v4/graphql"

# Acima disto o sistema de arquivos recusa o nome (255 bytes no ext4), e a
# consulta dirigida de indicacao do openFDA passa de 250 so de termos. O corte
# e alto de proposito: nenhum nome ja gravado chega perto, entao ligar isto nao
# renomeia fixture nenhuma.
LIMITE_NOME = 150


def arquivo_da_resposta(url: str, corpo: dict | None = None) -> str:
    """Nome do arquivo gravado para uma requisicao: a URL sem a base, com tudo
    que nao for [A-Za-z0-9._=-] trocado por '_'. Deterministico e legivel.
    Open Targets e POST num unico endereco: o nome sai da operacao GraphQL e
    de suas variaveis.

    Nome comprido demais e cortado e assinado com 12 digitos do SHA-256 do
    nome inteiro: fica legivel no inicio e continua injetivo."""
    import hashlib
    import re
    if url == URL_OPEN_TARGETS:
        assert corpo, "Open Targets exige o corpo GraphQL"
        operacao = re.search(r"query\s+(\w+)", corpo["query"]).group(1)
        variaveis = "_".join(f"{k}={v}" for k, v in corpo.get("variables", {}).items())
        nome = re.sub(r"[^A-Za-z0-9._=-]", "_", f"opentargets_{operacao}_{variaveis}")
    else:
        base = next((b for b in (BASE_CHEMBL, BASE_REACTOME, BASE_OPENFDA)
                     if url.startswith(b)), None)
        assert base, url
        nome = re.sub(r"[^A-Za-z0-9._=-]", "_", url[len(base):])
        if base == BASE_REACTOME:
            nome = "reactome_" + nome
        elif base == BASE_OPENFDA:
            nome = "openfda_" + nome
    nome = nome if nome.endswith(".json") else nome + ".json"
    if len(nome) > LIMITE_NOME:
        tronco = nome[:-len(".json")]
        assinatura = hashlib.sha256(tronco.encode("utf-8")).hexdigest()[:12]
        nome = f"{tronco[:LIMITE_NOME - len('.json') - 13]}_{assinatura}.json"
    return nome


@pytest.fixture
def respostas():
    """`obter_json` falso: le a resposta gravada em tests/dados/respostas/.
    URL sem arquivo gravado FALHA o teste; nunca cai para a rede. Aceita o
    corpo GraphQL opcional que o Open Targets usa."""
    import json
    from pathlib import Path
    pasta = Path(__file__).parent / "dados" / "respostas"
    vistas: dict[str, str] = {}

    def obter(url: str, corpo: dict | None = None) -> dict:
        nome = arquivo_da_resposta(url, corpo)
        # Injetividade: duas URLs distintas no mesmo arquivo falham alto.
        assert vistas.setdefault(nome, url) == url, (
            f"colisao: {url} e {vistas[nome]} -> {nome}")
        caminho = pasta / nome
        assert caminho.exists(), f"resposta nao gravada: {url} -> {caminho.name}"
        return json.loads(caminho.read_text(encoding="utf-8"))
    return obter
