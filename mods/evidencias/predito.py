# -*- coding: utf-8 -*-
"""
predito.py

Evidencia PREVISTA: probabilidade de o ligante atuar sobre o alvo, dada pelos
pesos pre-treinados do PP4CKD. Funciona sem rede; so precisa do .pt e do
target_labels.tsv.

Regras
------
* `pchembl` fica sempre vazio: e a unica coluna comparavel entre linhas e
  significa potencia MEDIDA. Previsao nao entra nela.
* Alvo ausente do modelo nao gera linha. Probabilidade zero afirmaria algo que
  o modelo nunca disse.
* `n_moleculas_treino` (n_ativos do alvo) acompanha toda linha: alvo nao humano
  e muito mais esparso, e 0.9 sobre 8 exemplos nao vale 0.9 sobre 4000.
* Entrega por `yield`, em lotes de ligantes: a carga real e ~108 mil ligantes
  vezes dezenas de alvos.
"""

from typing import Iterator

from mods.evidencias.esquema import formatar_atividade
from mods.prever_alvos import carregar_modelo, fingerprint_de

LOTE = 256


def _indices_dos_alvos(alvos: list[dict], rotulos: list[dict]) -> list[int]:
    """Indices do modelo para os alvos pedidos; ausentes sao omitidos."""
    por_id = {r["chembl_id"].upper(): i for i, r in enumerate(rotulos)
              if r["chembl_id"]}
    indices: list[int] = []
    for alvo in alvos:
        chembl_id = (alvo.get("alvo_chembl_id") or "").strip().upper()
        if chembl_id:
            i = por_id.get(chembl_id)
        else:
            # Sem ID: nome exato e organismo, porque o mesmo nome existe em
            # varias especies.
            nome = (alvo.get("alvo") or "").strip().casefold()
            org = (alvo.get("organismo") or "").strip().casefold()
            i = next((k for k, r in enumerate(rotulos)
                      if nome and r["nome"].casefold() == nome
                      and org and org in r["organismo"].casefold()), None)
        if i is not None and i not in indices:
            indices.append(i)
    return indices


def linhas_preditas(ligantes: list[dict], alvos: list[dict], *, checkpoint: str,
                    rotulos: list[dict], corte: float = 0.1,
                    contagens: dict | None = None) -> Iterator[dict]:
    """Uma linha por par ligante/alvo com probabilidade >= `corte`.

    `contagens`, se dado, recebe `smiles_ilegiveis`: ligantes cujo fingerprint
    nao pode ser calculado e que por isso nunca foram pontuados."""
    if contagens is not None:
        contagens.setdefault("smiles_ilegiveis", 0)
    import numpy as np
    import torch

    indices = _indices_dos_alvos(alvos, rotulos)
    if not indices:
        return
    rede, cfg, dispositivo, nome_ckpt = carregar_modelo(
        None, "ECFP4", checkpoint=checkpoint)
    n_saida = rede[-1].out_features
    indices = [i for i in indices if i < n_saida]
    procedencia = f"predito (PP4CKD ECFP4 {nome_ckpt})"

    for inicio in range(0, len(ligantes), LOTE):
        validos, vetores = [], []
        for lig in ligantes[inicio:inicio + LOTE]:
            fp = fingerprint_de(lig.get("smiles", ""), cfg)
            if fp is None:
                if contagens is not None:
                    contagens["smiles_ilegiveis"] += 1
                continue
            validos.append(lig)
            vetores.append(fp)
        if not validos:
            continue
        with torch.no_grad():
            x = torch.from_numpy(np.asarray(vetores, dtype="float32")).to(dispositivo)
            probs = torch.sigmoid(rede(x)[:, indices]).cpu().tolist()
        for lig, linha_p in zip(validos, probs):
            for i, p in zip(indices, linha_p):
                if p < corte:
                    continue
                r = rotulos[i]
                yield {
                    "ligante": lig.get("nome", ""),
                    "ligante_inchikey_parental": lig.get("inchikey_parental", ""),
                    "ligante_chembl_id": lig.get("chembl_id", ""),
                    "alvo_chembl_id": r["chembl_id"],
                    "alvo": r["nome"],
                    "organismo": r["organismo"],
                    "tipo_alvo": r["tipo"],
                    "tipo_evidencia": "predito",
                    "fonte": "PP4CKD",
                    "atividade_relatada": formatar_atividade(
                        "p", "=", round(p, 2), "(ECFP4)"),
                    # Numero de verdade (nao o texto de exibicao): a fase de
                    # ranking ordena por ele; 6 casas preservam a ordem.
                    "medida": "probabilidade",
                    "relacao": "=",
                    "valor": round(p, 6),
                    "unidade": "",
                    "pchembl": "",
                    "n_moleculas_treino": r["n_ativos"],
                    "procedencia": procedencia,
                }
