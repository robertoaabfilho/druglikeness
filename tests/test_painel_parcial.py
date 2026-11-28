# -*- coding: utf-8 -*-
"""
T2 -- ler o que a etapa ainda esta escrevendo.

A publicacao e atomica: a tabela so aparece no destino quando a etapa termina
inteira. Durante a execucao o que existe e o temporario que `atomico.temporario`
criou AO LADO do destino (`<nome>.<aleatorio>.tmp`). O painel copia esse
temporario e le a copia -- nunca o arquivo que a etapa tem aberto.

A ultima linha da copia e descartada SEMPRE. A etapa esta escrevendo; julgar
se a ultima linha "parece completa" e adivinhar, e `ler_tsv` levanta
`ValueError` quando o numero de campos nao bate. Uma linha a menos num
resultado parcial nao custa nada.
"""
import pytest

from mods import caminhos
from mods.evidencias import atomico
from mods.painel import leitura


@pytest.fixture
def etapa1(tmp_path, monkeypatch):
    """`tsv/` numa pasta temporaria, com a pasta da etapa 1 criada."""
    monkeypatch.setattr(caminhos, "TSV", tmp_path)
    (tmp_path / "etapa1").mkdir()
    return tmp_path / "etapa1"


def _temporario(pasta, texto):
    """Um temporario como a etapa o cria: ao lado do destino, pelo `atomico`."""
    tmp = atomico.temporario(pasta / "evidencias.tsv")
    tmp.write_text(texto, encoding="utf-8")
    return tmp


def test_temporario_da_etapa_devolve_none_sem_execucao(etapa1):
    assert leitura.temporario_da_etapa("etapa1") is None


def test_temporario_da_etapa_encontra_o_tmp(etapa1):
    tmp = _temporario(etapa1, "#a\tb\nA\t1\n")
    assert leitura.temporario_da_etapa("etapa1") == tmp


def test_temporario_da_etapa_ignora_tmp_de_outra_tabela(etapa1):
    """A etapa 1 publica tres arquivos. O parcial da tela e o da TABELA; o
    temporario do manifesto nao tem linhas para mostrar."""
    tmp_manifesto = atomico.temporario(etapa1 / "manifesto.json")
    tmp_manifesto.write_text("{}", encoding="utf-8")
    assert leitura.temporario_da_etapa("etapa1") is None


def test_copiar_parcial_descarta_a_ultima_linha(etapa1, tmp_path):
    """A ultima linha esta truncada no meio de um campo: e exatamente o que a
    etapa deixa no disco entre dois `write`."""
    _temporario(etapa1, "#ligante\tvalor\nA\t1\nB\t2\nC\t")
    copia = tmp_path / "copia.tsv"
    assert leitura.copiar_parcial("etapa1", copia) == 2
    assert leitura.ler_tabela(copia)["ligante"].tolist() == ["A", "B"]


def test_copiar_parcial_descarta_a_ultima_linha_mesmo_completa(etapa1,
                                                               tmp_path):
    """Nao ha como saber que a linha terminou: o `\\n` pode ter sido gravado e
    o resto da linha ainda estar a caminho no buffer do sistema. Descartar
    sempre e a unica regra que nao depende de adivinhar."""
    _temporario(etapa1, "#ligante\tvalor\nA\t1\nB\t2\n")
    copia = tmp_path / "copia.tsv"
    assert leitura.copiar_parcial("etapa1", copia) == 1


def test_copiar_parcial_de_arquivo_so_com_cabecalho_devolve_zero(etapa1,
                                                                 tmp_path):
    _temporario(etapa1, "#ligante\tvalor\n")
    copia = tmp_path / "copia.tsv"
    assert leitura.copiar_parcial("etapa1", copia) == 0
    assert leitura.ler_tabela(copia) is not None


def test_copiar_parcial_de_arquivo_vazio_devolve_zero(etapa1, tmp_path):
    """O temporario nasce vazio: `mkstemp` cria antes de a etapa escrever o
    cabecalho. Zero linhas, nao erro."""
    _temporario(etapa1, "")
    copia = tmp_path / "copia.tsv"
    assert leitura.copiar_parcial("etapa1", copia) == 0


def test_copiar_parcial_sem_temporario_devolve_none(etapa1, tmp_path):
    assert leitura.copiar_parcial("etapa1", tmp_path / "copia.tsv") is None


def test_copiar_parcial_nao_toca_no_temporario(etapa1, tmp_path):
    tmp = _temporario(etapa1, "#ligante\tvalor\nA\t1\nB\t2\n")
    antes = tmp.read_bytes()
    leitura.copiar_parcial("etapa1", tmp_path / "copia.tsv")
    assert tmp.read_bytes() == antes


def test_copiar_parcial_respeita_o_limite(etapa1, tmp_path):
    """Copiar e barato, carregar nao: o teto existe para a copia tambem, para
    que uma etapa 1 de 23 milhoes de pares nao se copie inteira a cada dois
    segundos de acompanhamento."""
    _temporario(etapa1, "#a\tb\n" + "".join(f"{i}\t{i}\n" for i in range(100)))
    copia = tmp_path / "copia.tsv"
    assert leitura.copiar_parcial("etapa1", copia, limite=10) == 10


def test_estado_da_etapa_e_parcial_quando_ha_temporario(etapa1):
    _temporario(etapa1, "#ligante\tvalor\nA\t1\n")
    assert leitura.estado_da_etapa("etapa1")["estado"] == "parcial"


def test_tabela_publicada_vence_o_temporario_de_outra_etapa(etapa1,
                                                            monkeypatch):
    """Tabela pronta E temporario ao lado significa que uma execucao nova esta
    reescrevendo a etapa: o estado e `parcial`, porque e a execucao em curso
    que o pesquisador esta acompanhando, e a tabela publicada ainda e a
    anterior."""
    (etapa1 / "evidencias.tsv").write_text("#ligante\nA\n", encoding="utf-8")
    _temporario(etapa1, "#ligante\nB\n")
    estado = leitura.estado_da_etapa("etapa1")
    assert estado["estado"] == "parcial"
    assert estado["tabela"] is not None  # a anterior continua citavel
