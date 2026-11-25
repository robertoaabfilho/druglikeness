# -*- coding: utf-8 -*-
"""
resumo.py

Contagens que o painel mostra assim que um arquivo e enviado, antes de rodar
qualquer coisa. Ficam fora de `painel.py` para serem testaveis sem Streamlit.

`resumo_alvos` usa a mesma regra de casamento da evidencia prevista
(`predito.casar_alvos`): o que o painel chama "fora do modelo" e exatamente o
que a previsao vai omitir, e o que ele chama "ambiguo" e o que ela se recusa a
escolher.
"""

import hashlib
import json
from pathlib import Path

from mods.evidencias.entrada import Problema
from mods.evidencias.esquema import desescapar
from mods.evidencias.open_targets import COLUNAS_ENRIQUECIDOS
from mods.evidencias.predito import casar_alvos
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
    `ambiguo`: o nome+organismo casa com mais de um alvo do modelo, e por isso
    nenhum e escolhido (nao e "fora do modelo": e indecidivel sem o ID).
    `pouco_treino`: existe no modelo com menos de `minimo_treino` ativos."""
    sem_organismo, fora, fracos, ambiguos = [], [], [], []
    for a in alvos:
        if not (a.get("alvo_chembl_id") or "").strip() and not (
                a.get("organismo") or "").strip():
            sem_organismo.append(a)
            continue
        duvidas: list[str] = []
        indices = [i for i, _ in casar_alvos([a], rotulos, duvidas)]
        if duvidas:
            ambiguos.append({**a, "motivo": duvidas[0]})
        elif not indices:
            fora.append(a)
        elif rotulos[indices[0]]["n_ativos"] < minimo_treino:
            fracos.append({**a, "alvo_chembl_id": rotulos[indices[0]]["chembl_id"],
                           "n_ativos": rotulos[indices[0]]["n_ativos"]})
    ja_descartados = sum(p.campo == "organismo" and p.gravidade == "erro"
                         for p in problemas)
    return {"total": len(alvos) + ja_descartados,
            "sem_organismo": len(sem_organismo) + ja_descartados,
            "fora_do_modelo": len(fora), "pouco_treino": len(fracos),
            "ambiguo": len(ambiguos),
            "no_modelo": len(alvos) - len(sem_organismo) - len(fora) - len(ambiguos),
            "lista_sem_organismo": sem_organismo, "lista_fora_do_modelo": fora,
            "lista_ambiguo": ambiguos,
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


def chave_enriquecimento(chave_entrada: str, efos: dict[str, str],
                         corte_associacao: float, teto_via: int,
                         limite: int | None) -> str:
    """Identifica um enriquecimento pelo que o determina: os tres arquivos
    enviados (`chave_entrada`, MD5 deles), os EFO escolhidos e os parametros.
    Nomeia a pasta da execucao: lista de outra escolha nunca e confundida com
    a atual."""
    texto = json.dumps([chave_entrada, sorted(efos.items()), corte_associacao,
                        teto_via, limite or 0])
    return hashlib.md5(texto.encode("utf-8")).hexdigest()[:12]


def chave_coleta(chave_entrada: str, assinatura_dos_alvos: str, corte: float,
                 limite: int | None, pesos: str, medido: str) -> str:
    """Identifica uma coleta pelo que a determina: os tres arquivos enviados
    (`chave_entrada`), a lista de alvos realmente usada (a enviada ou a
    revisada, por conteudo), o corte, o limite e as fontes de evidencia.

    Nomeia a pasta da execucao, como `chave_enriquecimento` faz com a lista
    enriquecida. Sem isso o painel lia `coleta/evidencias.tsv` num caminho
    fixo e mostrava a tabela de uma execucao anterior como se fosse o
    resultado dos arquivos que estao na tela."""
    texto = json.dumps([chave_entrada, assinatura_dos_alvos, corte, limite or 0,
                        pesos, medido])
    return hashlib.md5(texto.encode("utf-8")).hexdigest()[:12]


def assinatura_de_arquivo(caminho) -> str:
    """Identidade barata de um arquivo grande (caminho, tamanho e mtime), para
    entrar numa chave de pasta sem ler 785 MB de pesos a cada rerun. O MD5 do
    conteudo entra no manifesto, que e quem precisa da garantia forte."""
    if not caminho:
        return ""
    p = Path(caminho).expanduser()
    try:
        st = p.stat()
    except OSError:
        return f"{p}:ausente"
    return f"{p}:{st.st_size}:{int(st.st_mtime)}"


def chave_do_editor(pasta: str, geracao: int) -> str:
    """Chave do `st.data_editor` da lista de alvos. A geracao existe para o
    botao de descartar: subir a geracao cria um widget NOVO, cujo estado nasce
    vazio -- em vez de depender de remover o estado do widget antigo do
    session_state, que o Streamlit reenvia no rerun."""
    return f"editor_alvos_{pasta}_g{geracao}"


def edicao_pendente(estado_do_editor) -> bool:
    """O `st.data_editor` guarda suas edicoes em session_state como
    {edited_rows, added_rows, deleted_rows}. True se ha alguma: o painel nao
    pode redesenhar a pagina por cima do trabalho manual do pesquisador."""
    if not isinstance(estado_do_editor, dict):
        return False
    return any(estado_do_editor.get(k) for k in
               ("edited_rows", "added_rows", "deleted_rows"))


def ler_tabela_de_alvos(caminho) -> tuple[list[str], list[list[str]]]:
    """(colunas, linhas) de um alvos_enriquecidos.tsv. Arquivo vazio ou so com
    cabecalho nao e erro: um enriquecimento sem resultado e plausivel."""
    linhas = Path(caminho).read_text(encoding="utf-8").splitlines()
    if not linhas:
        return list(COLUNAS_ENRIQUECIDOS), []
    colunas = linhas[0].lstrip("#").split("\t")
    return colunas, [[desescapar(c) for c in l.split("\t")] for l in linhas[1:] if l]
