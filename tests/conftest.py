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
