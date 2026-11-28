# -*- coding: utf-8 -*-
"""
T7 -- revisar a lista de alvos enriquecidos. FERRAMENTA, nao portao.

A etapa 1 publica `alvos_enriquecidos.tsv`: os alvos que o pesquisador
declarou mais os que o Open Targets e o Reactome acrescentaram. Editar essa
lista estreita ou amplia o que a fase 3 considera alvo afirmativo -- e
**nao reprova ligante nenhum**. Um ligante que perde um alvo pode ter outros;
se nao tiver, ele cai para a faixa 3 e aparece em `nao_ranqueados` com o
motivo, nunca descartado em silencio.

Como o upload, a gravacao **diz o que mudou**, e pelo mesmo motivo: uma lista
que o pesquisador acredita estar valendo e nao esta e pior que uma mensagem.
"""
import pandas as pd
import pytest

from mods.evidencias.open_targets import COLUNAS_ENRIQUECIDOS
from mods.painel import alvos as palvos


def _lista(ids=("CHEMBL233", "CHEMBL1951")) -> pd.DataFrame:
    linhas = []
    for i, cid in enumerate(ids):
        linha = {c: "" for c in COLUNAS_ENRIQUECIDOS}
        linha.update({"alvo_chembl_id": cid, "alvo": f"Alvo {i}",
                      "organismo": "Homo sapiens", "origem": "usuario",
                      "no_modelo": "sim"})
        linhas.append(linha)
    return pd.DataFrame(linhas)


def _gravado(tmp_path, df=None):
    caminho = tmp_path / "alvos_enriquecidos.tsv"
    palvos.gravar_enriquecidos(df if df is not None else _lista(), caminho)
    return caminho


# --------------------------------------------------------------------------
# caminho e leitura
# --------------------------------------------------------------------------
def test_caminho_vem_da_configuracao():
    cfg = {"etapa1": {"alvos_enriquecidos": "tsv/etapa1/alvos_enriquecidos.tsv"}}
    assert palvos.caminho_enriquecidos(cfg).name == "alvos_enriquecidos.tsv"


def test_sem_a_chave_na_configuracao_nao_ha_lista():
    """O enriquecimento e opcional: sem a chave, a etapa 1 nao escreve lista
    nenhuma, e o bloco nao aparece em vez de apontar para um caminho fixo."""
    assert palvos.caminho_enriquecidos({"etapa1": {}}) is None


def test_ler_enriquecidos_de_arquivo_ausente_devolve_none(tmp_path):
    assert palvos.ler_enriquecidos(tmp_path / "nao_existe.tsv") is None


def test_ler_de_volta_o_que_gravou(tmp_path):
    caminho = _gravado(tmp_path)
    lido = palvos.ler_enriquecidos(caminho)
    assert lido["alvo_chembl_id"].tolist() == ["CHEMBL233", "CHEMBL1951"]


# --------------------------------------------------------------------------
# gravar
# --------------------------------------------------------------------------
def test_gravar_relata_o_que_saiu(tmp_path):
    caminho = _gravado(tmp_path)
    relato = palvos.gravar_enriquecidos(_lista(("CHEMBL233",)), caminho)
    assert relato["removidos"] == ["CHEMBL1951"]
    assert relato["acrescentados"] == []
    assert relato["total"] == 1


def test_gravar_relata_o_que_entrou(tmp_path):
    caminho = _gravado(tmp_path, _lista(("CHEMBL233",)))
    relato = palvos.gravar_enriquecidos(_lista(("CHEMBL233", "CHEMBL4")),
                                        caminho)
    assert relato["acrescentados"] == ["CHEMBL4"]


def test_gravar_relata_o_md5_anterior(tmp_path):
    caminho = _gravado(tmp_path)
    relato = palvos.gravar_enriquecidos(_lista(("CHEMBL233",)), caminho)
    assert relato["md5_anterior"] and relato["md5_anterior"] != relato["md5"]


def test_gravar_preserva_a_ordem_das_colunas(tmp_path):
    """A etapa 2 le esta tabela pelo nome das colunas; a ordem e a do esquema
    da etapa 1, e colunas extras do editor vao para o fim."""
    caminho = _gravado(tmp_path)
    cabecalho = caminho.read_text(encoding="utf-8").splitlines()[0]
    assert cabecalho.startswith("#" + "\t".join(COLUNAS_ENRIQUECIDOS))


def test_gravar_mantem_coluna_extra_no_fim(tmp_path):
    df = _lista()
    df["minha_anotacao"] = ["x", "y"]
    caminho = _gravado(tmp_path, df)
    cabecalho = caminho.read_text(encoding="utf-8").splitlines()[0]
    assert cabecalho.endswith("minha_anotacao")


def test_gravar_recusa_lista_vazia(tmp_path):
    """Lista sem nenhum alvo jogaria a fase 3 inteira para a faixa 3: nada
    ranquearia, e a tela diria "rodou e nada qualificou" sobre um resultado
    que foi causado por uma linha apagada sem querer."""
    caminho = _gravado(tmp_path)
    with pytest.raises(ValueError):
        palvos.gravar_enriquecidos(_lista(()), caminho)


def test_gravar_recusa_linha_sem_identidade(tmp_path):
    """Sem `alvo_chembl_id` o alvo nao se liga a evidencia medida nem ao
    rotulo do modelo: a linha nao serve para nada, e gravada em silencio ela
    pareceria um alvo a mais na contagem."""
    caminho = _gravado(tmp_path)
    df = _lista()
    df.loc[0, "alvo_chembl_id"] = ""
    with pytest.raises(ValueError):
        palvos.gravar_enriquecidos(df, caminho)


def test_gravar_recusa_id_repetido(tmp_path):
    caminho = _gravado(tmp_path)
    with pytest.raises(ValueError, match="CHEMBL233"):
        palvos.gravar_enriquecidos(_lista(("CHEMBL233", "CHEMBL233")), caminho)


def test_gravar_nao_deixa_a_lista_pela_metade(tmp_path, monkeypatch):
    caminho = _gravado(tmp_path)
    antes = caminho.read_bytes()

    def explode(*a, **k):
        raise OSError("disco cheio")
    monkeypatch.setattr(palvos.atomico, "publicar", explode)
    with pytest.raises(OSError):
        palvos.gravar_enriquecidos(_lista(("CHEMBL233",)), caminho)
    assert caminho.read_bytes() == antes


def test_gravar_em_lugar_novo_nao_exige_anterior(tmp_path):
    relato = palvos.gravar_enriquecidos(_lista(), tmp_path / "nova.tsv")
    assert relato["md5_anterior"] is None and relato["total"] == 2


# --------------------------------------------------------------------------
# diferenca
# --------------------------------------------------------------------------
def test_diferenca_identifica_por_chembl_id():
    d = palvos.diferenca(_lista(("CHEMBL233", "CHEMBL1951")),
                         _lista(("CHEMBL233", "CHEMBL4")))
    assert d["removidos"] == ["CHEMBL1951"] and d["acrescentados"] == ["CHEMBL4"]


def test_diferenca_sem_lista_anterior():
    d = palvos.diferenca(None, _lista(("CHEMBL233",)))
    assert d["acrescentados"] == ["CHEMBL233"] and d["removidos"] == []


def test_diferenca_de_listas_iguais_e_vazia():
    d = palvos.diferenca(_lista(), _lista())
    assert d["removidos"] == [] and d["acrescentados"] == []


def test_editar_a_lista_nao_reprova_ligante():
    """A regra que faz disto ferramenta e nao portao: nada nesta API toca em
    ligante, em `portao.tsv` ou em veredito. Se algum dia tocar, este teste e
    o lugar onde o nome aparece."""
    import inspect
    fonte = inspect.getsource(palvos)
    for proibido in ("portao", "aprovado", "reprova", "ligante"):
        assert f"def {proibido}" not in fonte


# --------------------------------------------------------------------------
# a casca: o bloco fica em Operar, porque e escrita
# --------------------------------------------------------------------------
def test_o_editor_aparece_quando_a_lista_existe(tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest
    from mods import caminhos, config
    monkeypatch.setattr(caminhos, "TSV", tmp_path)
    monkeypatch.setattr(caminhos, "EXECUCAO", tmp_path / ".execucao")
    lista = tmp_path / "alvos_enriquecidos.tsv"
    palvos.gravar_enriquecidos(_lista(), lista)
    verdadeiro = config.ler

    def com_a_lista(caminho=None):
        cfg = verdadeiro(caminho)
        cfg["etapa1"]["alvos_enriquecidos"] = str(lista)
        return cfg
    monkeypatch.setattr(config, "ler", com_a_lista)
    at = AppTest.from_file(str(caminhos.RAIZ / "painel.py"), default_timeout=30)
    at.run()
    assert not at.exception
    assert at.get("data_editor"), "o editor da lista de alvos nao apareceu"


def test_sem_a_lista_o_painel_abre_igual(tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest
    from mods import caminhos
    monkeypatch.setattr(caminhos, "TSV", tmp_path)
    monkeypatch.setattr(caminhos, "EXECUCAO", tmp_path / ".execucao")
    at = AppTest.from_file(str(caminhos.RAIZ / "painel.py"), default_timeout=30)
    at.run()
    assert not at.exception


def test_caminho_relativo_resolve_contra_a_raiz_do_projeto():
    """Caminho relativo na configuracao vale a partir da raiz, como o main.py
    entrega as CLIs (`cwd=RAIZ`). Resolver contra o `cwd` faria o painel
    mostrar a lista de um lugar e a etapa escrever em outro."""
    from mods import caminhos
    cfg = {"etapa1": {"alvos_enriquecidos": "tsv/etapa1/alvos_enriquecidos.tsv"}}
    assert palvos.caminho_enriquecidos(cfg) == \
        caminhos.RAIZ / "tsv/etapa1/alvos_enriquecidos.tsv"


def test_caminho_absoluto_e_respeitado(tmp_path):
    cfg = {"etapa1": {"alvos_enriquecidos": str(tmp_path / "x.tsv")}}
    assert palvos.caminho_enriquecidos(cfg) == tmp_path / "x.tsv"


def test_gravar_recusa_linha_em_branco_do_editor(tmp_path):
    """I1: `st.data_editor(num_rows="dynamic")` devolve a linha nova com
    `None`. `str(None)` e `"None"`, que nao e vazio, e a guarda passava: o
    arquivo saia com um alvo de ID vazio e a tela dizia "2 alvos gravados"."""
    caminho = _gravado(tmp_path, _lista(("CHEMBL233",)))
    df = _lista(("CHEMBL233",))
    nova = {c: None for c in df.columns}
    df = pd.concat([df, pd.DataFrame([nova])], ignore_index=True)
    with pytest.raises(ValueError):
        palvos.gravar_enriquecidos(df, caminho)
    assert palvos.ler_enriquecidos(caminho)["alvo_chembl_id"].tolist() == \
        ["CHEMBL233"]


def test_gravar_recusa_linha_com_nan(tmp_path):
    """O mesmo caso quando o pandas converte a celula vazia para NaN."""
    import numpy as np
    caminho = _gravado(tmp_path, _lista(("CHEMBL233",)))
    df = _lista(("CHEMBL233",))
    nova = {c: np.nan for c in df.columns}
    df = pd.concat([df, pd.DataFrame([nova])], ignore_index=True)
    with pytest.raises(ValueError):
        palvos.gravar_enriquecidos(df, caminho)


def test_diferenca_nao_conta_none_como_alvo(tmp_path):
    """`acrescentados: ['None']` era o que a tela mostrava."""
    df = _lista(("CHEMBL233",))
    df = pd.concat([df, pd.DataFrame([{c: None for c in df.columns}])],
                   ignore_index=True)
    assert palvos.diferenca(_lista(("CHEMBL233",)), df)["acrescentados"] == []
