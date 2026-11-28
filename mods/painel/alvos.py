# -*- coding: utf-8 -*-
"""
alvos.py - revisar a lista de alvos enriquecidos. Ferramenta, nao portao.

A etapa 1 publica `alvos_enriquecidos.tsv`: os alvos que o pesquisador
declarou mais os que o Open Targets e o Reactome acrescentaram. Esta lista e
o que a fase 3 considera alvo afirmativo, e e por isso que ela e revisavel: o
enriquecimento e automatico, e um alvo que entrou por via metabolica pode nao
ser o que o pesquisador queria dizer.

**Editar nao reprova ligante nenhum.** Tirar um alvo estreita o conjunto de
alvos afirmativos. Um ligante que perde um alvo pode ter outros; se nao tiver,
ele cai para a faixa 3 e aparece em `nao_ranqueados` COM O MOTIVO -- nunca
descartado em silencio. Nada aqui toca em `portao.tsv` nem em veredito.

**A gravacao diz o que mudou**, como o upload das entradas: quais IDs sairam,
quais entraram, e o MD5 de antes e de depois. Uma lista que o pesquisador
acredita estar valendo e nao esta e pior que uma mensagem.
"""
import hashlib
from pathlib import Path

from mods import caminhos
from mods.evidencias import atomico
from mods.evidencias.open_targets import COLUNAS_ENRIQUECIDOS
from mods.portao.esquema import EscritorTSV, ler_tsv

CHAVE = "alvo_chembl_id"


def _md5(caminho: Path) -> str | None:
    if not caminho.is_file():
        return None
    return hashlib.md5(caminho.read_bytes()).hexdigest()


def caminho_enriquecidos(cfg: dict) -> Path | None:
    """Onde a lista sai, pela chave `etapa1.alvos_enriquecidos`.

    O caminho vem da configuracao e nao de um valor fixo porque o
    enriquecimento e OPCIONAL: sem a chave, a etapa 1 nao escreve lista
    nenhuma, e apontar para um caminho fixo faria a tela prometer um arquivo
    que nunca vai existir.

    Caminho relativo na configuracao e resolvido contra a RAIZ do projeto, nao
    contra o diretorio de onde o painel foi chamado: e assim que o `main.py`
    o entrega as CLIs (`cwd=RAIZ`), e um painel que dependesse do `cwd`
    mostraria a lista de um lugar e a etapa escreveria em outro."""
    valor = (cfg.get("etapa1") or {}).get("alvos_enriquecidos")
    if not valor:
        return None
    caminho = Path(valor).expanduser()
    return caminho if caminho.is_absolute() else caminhos.RAIZ / caminho


def ler_enriquecidos(caminho):
    """A lista como DataFrame; `None` se o arquivo nao existe."""
    import pandas as pd
    caminho = Path(caminho)
    if not caminho.is_file():
        return None
    linhas = list(ler_tsv(caminho))
    colunas = list(linhas[0]) if linhas else list(COLUNAS_ENRIQUECIDOS)
    return pd.DataFrame(linhas, columns=colunas, dtype=str)


def _texto_dos_ids(df) -> list[str]:
    """A coluna de ID como texto, com celula ausente virando string vazia.

    `str(None)` e `"None"` e `str(nan)` e `"nan"`: passar a coluna por `str`
    direto fazia a linha em branco que o editor cria parecer um alvo chamado
    "None". A ausencia tem de ser reconhecida ANTES da conversao."""
    if df is None or CHAVE not in getattr(df, "columns", []):
        return []
    coluna = df[CHAVE]
    presentes = coluna.notna() if hasattr(coluna, "notna") else None
    if presentes is None:
        return [str(v).strip() for v in coluna]
    return [str(v).strip() if ok else ""
            for v, ok in zip(coluna.tolist(), presentes.tolist())]


def _ids(df) -> list[str]:
    return [i for i in _texto_dos_ids(df) if i]


def diferenca(antes, depois) -> dict:
    """Quais alvos sairam e quais entraram, por `alvo_chembl_id`.

    Pelo ID e nao pelo nome: o nome do alvo varia entre fontes, e o ID e a
    chave que liga a lista a evidencia medida e ao rotulo do modelo."""
    de, para = set(_ids(antes)), set(_ids(depois))
    return {"removidos": sorted(de - para), "acrescentados": sorted(para - de)}


def _colunas_na_ordem(df) -> list[str]:
    """As do esquema da etapa 1 primeiro, as acrescentadas pelo editor depois.

    A etapa 2 le esta tabela pelo nome das colunas, nao pela posicao, mas a
    ordem e o que faz o arquivo continuar legivel para quem o abre a mao."""
    presentes = [c for c in COLUNAS_ENRIQUECIDOS if c in df.columns]
    return presentes + [c for c in df.columns if c not in presentes]


def gravar_enriquecidos(df, caminho) -> dict:
    """Grava a lista revisada e relata o que mudou.

    Recusa, antes de tocar no arquivo:

    * **lista vazia** -- jogaria a fase 3 inteira para a faixa 3, e a tela
      diria "rodou e nada qualificou" sobre um resultado causado por uma linha
      apagada sem querer;
    * **linha sem `alvo_chembl_id`** -- sem a chave o alvo nao se liga a
      evidencia medida nem ao rotulo do modelo, e em silencio ele pareceria um
      alvo a mais na contagem;
    * **ID repetido** -- duplicaria a evidencia daquele alvo no ranking.
    """
    caminho = Path(caminho)
    if df is None or not len(df):
        raise ValueError(
            "a lista de alvos nao pode ficar vazia: sem alvo afirmativo, "
            "nenhum ligante ranqueia e tudo cai para a faixa 3")
    if CHAVE not in df.columns:
        raise ValueError(f"a lista precisa da coluna {CHAVE!r}")
    ids = _texto_dos_ids(df)
    if any(not i for i in ids):
        vazias = [n + 1 for n, i in enumerate(ids) if not i]
        raise ValueError(
            f"ha linha sem {CHAVE} (linha(s) {', '.join(map(str, vazias))}): "
            "o alvo nao se ligaria a evidencia medida nem ao rotulo do "
            "modelo. Apague a linha ou preencha o ID.")
    repetidos = sorted({i for i in ids if ids.count(i) > 1})
    if repetidos:
        raise ValueError(
            f"{CHAVE} repetido: {', '.join(repetidos)} -- duplicaria a "
            "evidencia daquele alvo no ranking")

    anterior = ler_enriquecidos(caminho)
    md5_anterior = _md5(caminho)
    mudou = diferenca(anterior, df)
    colunas = _colunas_na_ordem(df)
    caminho.parent.mkdir(parents=True, exist_ok=True)
    tmp = atomico.temporario(caminho)
    try:
        with EscritorTSV(tmp, colunas) as escritor:
            for linha in df.to_dict("records"):
                escritor.escrever({c: linha.get(c) for c in colunas})
        atomico.publicar(tmp, caminho)
    finally:
        tmp.unlink(missing_ok=True)  # so o temporario, nunca a lista do usuario
    return {"removidos": mudou["removidos"],
            "acrescentados": mudou["acrescentados"],
            "total": len(ids), "md5_anterior": md5_anterior,
            "md5": _md5(caminho), "caminho": caminho}
