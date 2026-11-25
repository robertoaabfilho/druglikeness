# -*- coding: utf-8 -*-
"""Fixtures compartilhadas. Tarefas seguintes estendem este arquivo."""
import pytest

N_ALVOS_FALSO = 6


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
URL_OPEN_TARGETS = "https://api.platform.opentargets.org/api/v4/graphql"


def arquivo_da_resposta(url: str, corpo: dict | None = None) -> str:
    """Nome do arquivo gravado para uma requisicao: a URL sem a base, com tudo
    que nao for [A-Za-z0-9._=-] trocado por '_'. Deterministico e legivel.
    Open Targets e POST num unico endereco: o nome sai da operacao GraphQL e
    de suas variaveis."""
    import re
    if url == URL_OPEN_TARGETS:
        assert corpo, "Open Targets exige o corpo GraphQL"
        operacao = re.search(r"query\s+(\w+)", corpo["query"]).group(1)
        variaveis = "_".join(f"{k}={v}" for k, v in corpo.get("variables", {}).items())
        nome = re.sub(r"[^A-Za-z0-9._=-]", "_", f"opentargets_{operacao}_{variaveis}")
    else:
        base = next((b for b in (BASE_CHEMBL, BASE_REACTOME) if url.startswith(b)), None)
        assert base, url
        nome = re.sub(r"[^A-Za-z0-9._=-]", "_", url[len(base):])
        if base == BASE_REACTOME:
            nome = "reactome_" + nome
    return nome if nome.endswith(".json") else nome + ".json"


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
