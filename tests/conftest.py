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


def arquivo_da_resposta(url: str) -> str:
    """Nome do arquivo gravado para uma URL: a URL sem a base, com tudo que
    nao for [A-Za-z0-9._=-] trocado por '_'. Determinístico e legivel."""
    import re
    assert url.startswith(BASE_CHEMBL), url
    nome = re.sub(r"[^A-Za-z0-9._=-]", "_", url[len(BASE_CHEMBL):])
    return nome if nome.endswith(".json") else nome + ".json"


@pytest.fixture
def respostas():
    """`obter_json` falso: le a resposta gravada em tests/dados/respostas/.
    URL sem arquivo gravado FALHA o teste; nunca cai para a rede."""
    import json
    from pathlib import Path
    pasta = Path(__file__).parent / "dados" / "respostas"
    vistas: dict[str, str] = {}

    def obter(url: str) -> dict:
        nome = arquivo_da_resposta(url)
        # Injetividade: duas URLs distintas no mesmo arquivo falham alto.
        assert vistas.setdefault(nome, url) == url, (
            f"colisao: {url} e {vistas[nome]} -> {nome}")
        caminho = pasta / nome
        assert caminho.exists(), f"resposta nao gravada: {url} -> {caminho.name}"
        return json.loads(caminho.read_text(encoding="utf-8"))
    return obter
