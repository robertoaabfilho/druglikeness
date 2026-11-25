# -*- coding: utf-8 -*-
"""
resumo.py

Contagens que o painel mostra assim que um arquivo e enviado, antes de rodar
qualquer coisa. Ficam fora de `painel.py` para serem testaveis sem Streamlit.

`resumo_alvos` usa a mesma regra de casamento da evidencia prevista
(`predito.indices_dos_alvos`): o que o painel chama "fora do modelo" e
exatamente o que a previsao vai omitir.
"""

from mods.evidencias.entrada import Problema
from mods.evidencias.predito import indices_dos_alvos
from mods.importar_prospector import molecula

# Abaixo disto a previsao e fragil (README, validacao: recall@10 0,127 com <50
# ativos contra 0,477 com >=5000).
MINIMO_TREINO = 50


def resumo_ligantes(ligantes: list[dict], problemas: list[Problema]) -> dict:
    """`validos`: ligantes cujo SMILES o RDKit le. `invalidos`: SMILES ilegiveis
    na lista mais as linhas que a leitura ja descartou por erro de SMILES
    (ausente ou ilegivel); as duas origens nao se sobrepoem."""
    ilegiveis = sum(molecula(l.get("smiles") or "") is None for l in ligantes)
    descartados = sum(p.campo == "smiles" and p.gravidade == "erro" for p in problemas)
    return {"total": len(ligantes), "validos": len(ligantes) - ilegiveis,
            "invalidos": ilegiveis + descartados,
            "erros": sum(p.gravidade == "erro" for p in problemas),
            "avisos": sum(p.gravidade == "aviso" for p in problemas)}


def resumo_alvos(alvos: list[dict], rotulos: list[dict],
                 minimo_treino: int = MINIMO_TREINO,
                 problemas: list[Problema] = ()) -> dict:
    """Problemas de cada alvo pedido, para mostrar antes de rodar.
    `sem_organismo`: dado por nome sem organismo (erro na leitura). `ler_alvos`
    ja descarta essas linhas, entao `problemas` (os da leitura) tambem conta.
    `fora_do_modelo`: tem identificacao completa mas nao existe no modelo.
    `pouco_treino`: existe no modelo com menos de `minimo_treino` ativos."""
    sem_organismo, fora, fracos = [], [], []
    for a in alvos:
        if not (a.get("alvo_chembl_id") or "").strip() and not (
                a.get("organismo") or "").strip():
            sem_organismo.append(a)
            continue
        indices = indices_dos_alvos([a], rotulos)
        if not indices:
            fora.append(a)
        elif rotulos[indices[0]]["n_ativos"] < minimo_treino:
            fracos.append({**a, "alvo_chembl_id": rotulos[indices[0]]["chembl_id"],
                           "n_ativos": rotulos[indices[0]]["n_ativos"]})
    ja_descartados = sum(p.campo == "organismo" and p.gravidade == "erro"
                         for p in problemas)
    return {"total": len(alvos) + ja_descartados,
            "sem_organismo": len(sem_organismo) + ja_descartados,
            "fora_do_modelo": len(fora), "pouco_treino": len(fracos),
            "no_modelo": len(alvos) - len(sem_organismo) - len(fora),
            "lista_sem_organismo": sem_organismo, "lista_fora_do_modelo": fora,
            "lista_pouco_treino": fracos, "minimo_treino": minimo_treino}


def pares_previstos(n_ligantes: int, n_alvos_no_modelo: int) -> int:
    """Pares ligante x alvo que a evidencia prevista vai pontuar."""
    return n_ligantes * n_alvos_no_modelo


def montar_comando(modo: str, *, python: str, ligantes: str, alvos: str,
                   sintomas: str, saida: str, rotulos: str,
                   checkpoint: str | None = None,
                   evidencias_chembl: str | None = None,
                   alvos_enriquecidos: str | None = None,
                   efos: dict[str, str] | None = None, cache: str | None = None,
                   corte: float = 0.1, corte_associacao: float = 0.4,
                   teto_via: int = 150, limite: int | None = None) -> list[str]:
    """Linha de comando da CLI que o painel dispara. O painel nao tem logica
    propria de coleta: monta a mesma chamada que o pesquisador faria a mao.
    `enriquecer` nao usa pesos nem ensaios (so Open Targets e Reactome);
    `coletar` nao enriquece."""
    if modo not in ("enriquecer", "coletar"):
        raise ValueError(f"modo desconhecido: {modo!r}")
    cmd = [python, "-m", "mods.coletar_evidencias", "--ligantes", ligantes,
           "--alvos", alvos, "--sintomas", sintomas, "--saida", saida,
           "--rotulos", rotulos]
    if cache:
        cmd += ["--cache", cache]
    if limite:
        cmd += ["--limite", str(limite)]
    if modo == "enriquecer":
        if not alvos_enriquecidos:
            raise ValueError("enriquecer exige alvos_enriquecidos")
        cmd += ["--alvos-enriquecidos", alvos_enriquecidos,
                "--corte-associacao", str(corte_associacao),
                "--teto-via", str(teto_via)]
        for sintoma, efo in (efos or {}).items():
            cmd += ["--efo", f"{sintoma}={efo}"]
    else:
        cmd += ["--corte", str(corte)]
        if checkpoint:
            cmd += ["--checkpoint", checkpoint]
        if evidencias_chembl:
            cmd += ["--evidencias-chembl", evidencias_chembl]
    return cmd
