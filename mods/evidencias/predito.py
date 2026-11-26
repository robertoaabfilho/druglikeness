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
* `origem_alvo` acompanha a linha, igual a evidencia medida: so assim se sabe
  se o alvo veio do pesquisador ou do enriquecimento.
* Checkpoint e rotulos precisam descrever o MESMO espaco de alvos. Largura de
  saida diferente do numero de rotulos e RECUSADA, nomeando as duas: a coluna
  i de um nao e o alvo i do outro, e o README manda trocar o checkpoint
  (`ecfp4_fold_4.pt`), o que torna o par divergente o caso esperado.
* Alvo pedido por nome casa por nome EXATO e organismo EXATO. Nome que casa
  com mais de um rotulo nao escolhe nenhum: e relatado, como
  `chembl.resolver_alvo_por_nome` faz.
* Entrega por `yield`, em lotes de ligantes: a carga real e ~108 mil ligantes
  vezes dezenas de alvos.
"""

from typing import Iterator

from mods.evidencias import esquema
from mods.evidencias.esquema import formatar_atividade
from mods.prever_alvos import carregar_modelo, fingerprint_de

LOTE = 256


def casar_alvos(alvos: list[dict], rotulos: list[dict],
                problemas: list[str] | None = None,
                avisos: list[str] | None = None) -> list[tuple[int, dict]]:
    """[(indice no modelo, alvo pedido)], na ordem dos alvos; alvo ausente do
    modelo e omitido e indice repetido aparece uma vez so.

    O ID do ChEMBL ja chega normalizado da leitura (`entrada.ler_alvos`), que e
    o unico ponto onde se normaliza. Sem ID, casa por nome EXATO e organismo
    EXATO: casar organismo por trecho traria o ortologo errado em silencio, que
    e justamente o anti-objetivo do desenho. Nome+organismo que casa com mais
    de um rotulo e AMBIGUO: nenhum e escolhido e o caso vai para `problemas`."""
    por_id: dict[str, int] = {}
    por_nome: dict[tuple[str, str], list[int]] = {}
    for i, r in enumerate(rotulos):
        if r["chembl_id"]:
            por_id.setdefault(r["chembl_id"], i)
        por_nome.setdefault((r["nome"].strip().casefold(),
                             r["organismo"].strip().casefold()), []).append(i)
    pares: list[tuple[int, dict]] = []
    usados: set[int] = set()
    for alvo in alvos:
        chembl_id = (alvo.get("alvo_chembl_id") or "").strip()
        if chembl_id:
            i = por_id.get(chembl_id)
        else:
            nome = (alvo.get("alvo") or "").strip().casefold()
            org = (alvo.get("organismo") or "").strip().casefold()
            achados = por_nome.get((nome, org), []) if nome and org else []
            if len(achados) > 1:
                if problemas is not None:
                    problemas.append(
                        f"alvo {alvo.get('alvo')!r} ({alvo.get('organismo')}) casa "
                        f"com {len(achados)} alvos do modelo "
                        f"({', '.join(rotulos[k]['chembl_id'] for k in achados)}): "
                        "ambiguo, NENHUMA linha prevista; informe o alvo_chembl_id")
                continue
            i = achados[0] if achados else None
        if i is not None and i in usados:
            # Pedido repetido (por ID e por nome): usado uma vez, e dito como
            # aviso, nao como falha -- o mesmo tratamento da evidencia medida.
            if avisos is not None:
                avisos.append(
                    f"alvo {chembl_id or alvo.get('alvo')}: e o mesmo alvo de um "
                    "pedido anterior; usado uma vez so")
        elif i is not None:
            usados.add(i)
            pares.append((i, alvo))
    return pares



def linhas_preditas(ligantes: list[dict], alvos: list[dict], *, checkpoint: str,
                    rotulos: list[dict], corte: float = 0.1,
                    contagens: dict | None = None) -> Iterator[dict]:
    """Uma linha por par ligante/alvo com probabilidade >= `corte`.

    `contagens`, se dado, recebe `smiles_ilegiveis` (ligantes cujo fingerprint
    nao pode ser calculado e que por isso nunca foram pontuados) e `problemas`
    (alvo ambiguo, nenhum alvo no modelo).

    ValueError se a largura de saida do checkpoint nao for o numero de rotulos:
    o par nao descreve o mesmo espaco de alvos."""
    if contagens is not None:
        contagens.setdefault("smiles_ilegiveis", 0)
        contagens.setdefault("problemas", [])
        contagens.setdefault("avisos", [])
    problemas = None if contagens is None else contagens["problemas"]
    import numpy as np
    import torch

    pares = casar_alvos(alvos, rotulos, problemas,
                        None if contagens is None else contagens["avisos"])
    if not pares:
        if problemas is not None:
            problemas.append(
                f"nenhum dos {len(alvos)} alvos pedidos existe no espaco do "
                "modelo: NENHUMA linha prevista (zero linhas nao significa "
                "'sem atividade')")
        return
    indices = [i for i, _ in pares]
    origens = {i: (a.get("origem_alvo") or "") for i, a in pares}
    rede, cfg, dispositivo, nome_ckpt = carregar_modelo(
        None, "ECFP4", checkpoint=checkpoint)
    n_saida = rede[-1].out_features
    if n_saida != len(rotulos):
        # Larguras diferentes nao sao "alguns alvos a menos": a coluna i de um
        # nao e o alvo i do outro, e pontuar assim atribuiria a probabilidade
        # de um alvo ao nome de outro.
        raise ValueError(
            f"checkpoint e rotulos descrevem espacos de alvos diferentes: "
            f"{nome_ckpt} tem {n_saida} saidas e o arquivo de rotulos tem "
            f"{len(rotulos)} alvos. Use o target_labels.tsv do mesmo PP4CKD do "
            "checkpoint (o tarball publicado esta truncado e o README manda "
            "trocar o checkpoint; o par tem de ser trocado junto)")
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
                    "origem_alvo": origens[i],
                    "tipo_evidencia": esquema.TIPO_PREDITO,
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
