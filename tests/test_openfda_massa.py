# -*- coding: utf-8 -*-
"""
Plano C -- o acervo local do openFDA.

Por que existe
--------------
`alerta_bula` e um dos DOIS unicos criterios que reprovam no portao. Com
`--limite-bula 500` sobre 108.935 ligantes, ele sai `indeterminado` para
**99,5% da base**: o criterio existe e quase nao age. O teto era imposto pela
taxa do servico publico, nao pelo desenho -- e Roberto removeu esse motivo.

A ideia que faz isto caber em pouco codigo
------------------------------------------
O acervo responde as **mesmas URLs** que o `bula.py` monta. Assim o `bula.py`
nao muda uma linha: ele recebe outro `obter_json`, e toda a logica medida
(janela de 50 por data, consulta dirigida no acervo inteiro, desconto da
indicacao) continua valendo igual.

Nenhum teste aqui usa rede: o download recebe um baixador injetado, e o indice
e construido de arquivos locais.
"""
import gzip
import json
import sqlite3

import pytest

from mods import openfda_massa as massa


# --------------------------------------------------------------------------
# as particoes do download.json
# --------------------------------------------------------------------------
DOWNLOAD = {
    "results": {
        "drug": {
            "label": {
                "export_date": "2026-11-25",
                "total_records": 265009,
                "partitions": [
                    {"size_mb": "120.5", "records": "20000",
                     "display_name": "drug-label-0001-of-0014",
                     "file": "https://download.open.fda.gov/drug/label/"
                             "drug-label-0001-of-0014.json.zip"},
                    {"size_mb": "118.2", "records": "19000",
                     "display_name": "drug-label-0002-of-0014",
                     "file": "https://download.open.fda.gov/drug/label/"
                             "drug-label-0002-of-0014.json.zip"},
                ],
            }
        }
    }
}


def test_particoes_vem_do_download_json():
    partes = massa.particoes(DOWNLOAD)
    assert [p["nome"] for p in partes] == ["drug-label-0001-of-0014",
                                           "drug-label-0002-of-0014"]
    assert partes[0]["url"].endswith("drug-label-0001-of-0014.json.zip")
    assert partes[0]["mb"] == pytest.approx(120.5)


def test_particoes_traz_a_data_do_export():
    """E a versao do acervo: o openFDA nao tem numero de versao, e a data da
    carga e o que identifica o conteudo -- a mesma regra de `versao_openfda`."""
    assert massa.data_do_export(DOWNLOAD) == "2026-11-25"


def test_particoes_de_resposta_sem_label_e_erro():
    with pytest.raises(ValueError, match="drug/label"):
        massa.particoes({"results": {}})


# --------------------------------------------------------------------------
# baixar, com retomada
# --------------------------------------------------------------------------
def _zip_de(rotulos: list[dict]) -> bytes:
    """Uma particao como o openFDA a publica: zip com um .json dentro."""
    import io
    import zipfile
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("drug-label-0001.json",
                   json.dumps({"meta": {"last_updated": "2026-11-25"},
                               "results": rotulos}))
    return buffer.getvalue()


def _rotulo(generico="CIPROFLOXACIN", **secoes):
    base = {"openfda": {"generic_name": [generico],
                        "brand_name": [f"{generico.title()} XR"],
                        "substance_name": [generico]},
            "effective_time": "20260707"}
    base.update({k: [v] for k, v in secoes.items()})
    return base


def test_baixar_grava_cada_particao(tmp_path):
    dados = _zip_de([_rotulo()])
    chamadas = []

    def baixador(url, destino):
        chamadas.append(url)
        destino.write_bytes(dados)
    rel = massa.baixar(DOWNLOAD, tmp_path, baixador=baixador)
    assert len(chamadas) == 2
    assert rel["baixadas"] == 2 and rel["ja_tinha"] == 0
    assert (tmp_path / "drug-label-0001-of-0014.json.zip").exists()


def test_baixar_retoma_e_nao_rebaixa(tmp_path):
    dados = _zip_de([_rotulo()])

    def baixador(url, destino):
        destino.write_bytes(dados)
    massa.baixar(DOWNLOAD, tmp_path, baixador=baixador)
    chamadas = []

    def contando(url, destino):
        chamadas.append(url)
        destino.write_bytes(dados)
    rel = massa.baixar(DOWNLOAD, tmp_path, baixador=contando)
    assert chamadas == [] and rel["ja_tinha"] == 2


def test_baixar_refaz_particao_truncada(tmp_path):
    """Um zip que nao abre e download interrompido, nao acervo: retomar tem de
    refaze-lo, nao aceita-lo. Aceitar daria acervo incompleto apresentado
    como completo -- e `indeterminado` virando `sem alerta`."""
    (tmp_path / "drug-label-0001-of-0014.json.zip").write_bytes(b"truncado")
    dados = _zip_de([_rotulo()])
    chamadas = []

    def baixador(url, destino):
        chamadas.append(url)
        destino.write_bytes(dados)
    massa.baixar(DOWNLOAD, tmp_path, baixador=baixador)
    assert len(chamadas) == 2


def test_baixar_nao_deixa_particao_pela_metade(tmp_path):
    """Escrita atomica: falha no meio nao deixa um zip parcial que a retomada
    aceitaria depois."""
    def explode(url, destino):
        destino.write_bytes(b"meio")
        raise OSError("conexao caiu")
    with pytest.raises(OSError):
        massa.baixar(DOWNLOAD, tmp_path, baixador=explode)
    assert list(tmp_path.glob("*.json.zip")) == []


def test_baixar_registra_a_procedencia(tmp_path):
    def baixador(url, destino):
        destino.write_bytes(_zip_de([_rotulo()]))
    massa.baixar(DOWNLOAD, tmp_path, baixador=baixador)
    proc = json.loads((tmp_path / "procedencia.json").read_text(encoding="utf-8"))
    assert proc["export_date"] == "2026-11-25"
    assert proc["particoes"] == 2
    assert proc["total_records"] == 265009


# --------------------------------------------------------------------------
# indexar
# --------------------------------------------------------------------------
@pytest.fixture
def acervo(tmp_path):
    """Acervo de tres bulas, indexado: ciprofloxacino com neuropatia em secao
    forte, capsaicina indicada para a condicao, e um farmaco sem nada."""
    pasta = tmp_path / "bruto"
    pasta.mkdir()
    rotulos = [
        _rotulo("CIPROFLOXACIN",
                boxed_warning="Fluoroquinolones have been associated with "
                              "disabling peripheral neuropathy.",
                adverse_reactions="paresthesia was reported"),
        _rotulo("CAPSAICIN",
                indications_and_usage="for the management of neuropathic pain",
                warnings_and_cautions="application site paresthesia"),
        _rotulo("PARACETAMOL", warnings="hepatotoxicity"),
    ]
    (pasta / "p1.json.zip").write_bytes(_zip_de(rotulos))
    indice = tmp_path / "openfda.sqlite"
    rel = massa.indexar(pasta, indice, export_date="2026-11-25")
    return {"indice": indice, "relatorio": rel, "pasta": pasta}


def test_indexar_conta_as_bulas(acervo):
    assert acervo["relatorio"]["bulas"] == 3


def test_indexar_guarda_so_as_secoes_que_importam(acervo):
    """O acervo bruto tem 1,7 GB e a bula tem 34 campos; o indice guarda os que
    o criterio le. Guardar tudo multiplicaria o tamanho sem mudar resposta."""
    con = sqlite3.connect(acervo["indice"])
    colunas = {r[1] for r in con.execute("PRAGMA table_info(bulas)")}
    assert "boxed_warning" in colunas
    assert "indications_and_usage" in colunas
    assert "spl_unclassified_section" not in colunas


def test_indexar_registra_a_versao_do_acervo(acervo):
    con = sqlite3.connect(acervo["indice"])
    valor = con.execute(
        "SELECT valor FROM meta WHERE chave='export_date'").fetchone()[0]
    assert valor == "2026-11-25"


def test_indexar_recusa_pasta_sem_particao(tmp_path):
    with pytest.raises(ValueError, match="nenhuma particao"):
        massa.indexar(tmp_path, tmp_path / "x.sqlite", export_date="2026-01-01")


def test_indexar_le_json_simples_e_gz(tmp_path):
    """Quem baixou a mao pode ter .json ou .json.gz; aceitar os tres evita
    exigir que o pesquisador recompacte 1,7 GB."""
    pasta = tmp_path / "bruto"
    pasta.mkdir()
    (pasta / "a.json").write_text(json.dumps({"results": [_rotulo("A")]}),
                                  encoding="utf-8")
    with gzip.open(pasta / "b.json.gz", "wt", encoding="utf-8") as f:
        json.dump({"results": [_rotulo("B")]}, f)
    rel = massa.indexar(pasta, tmp_path / "i.sqlite", export_date="x")
    assert rel["bulas"] == 2


# --------------------------------------------------------------------------
# o acervo responde as URLs do bula.py
# --------------------------------------------------------------------------
def test_a_versao_sai_do_acervo(acervo):
    from mods.portao import bula
    local = massa.Acervo(acervo["indice"])
    assert bula.versao_openfda(obter_json=local.obter_json) == \
        "openfda-massa 2026-11-25"


def test_buscar_bulas_acha_pelo_nome_generico(acervo):
    from mods.portao import bula
    local = massa.Acervo(acervo["indice"])
    r = bula.buscar_bulas("CIPROFLOXACIN", obter_json=local.obter_json)
    assert r["campo_do_nome"] == "generic_name"
    assert r["total"] == 1 and len(r["bulas"]) == 1


def test_buscar_bulas_acha_pelo_nome_de_marca(acervo):
    """Indexado por generico E por marca: e o que o plano pedia, e e como o
    pesquisador escreve o nome na lista de ligantes."""
    from mods.portao import bula
    local = massa.Acervo(acervo["indice"])
    r = bula.buscar_bulas("Ciprofloxacin XR", obter_json=local.obter_json)
    assert r["bulas"] and r["campo_do_nome"] == "brand_name"


def test_nome_desconhecido_devolve_vazio_e_nao_erro(acervo):
    from mods.portao import bula
    local = massa.Acervo(acervo["indice"])
    r = bula.buscar_bulas("NAO EXISTE", obter_json=local.obter_json)
    assert r["bulas"] == [] and r["total"] == 0


def test_o_criterio_inteiro_roda_contra_o_acervo_local(acervo):
    """O teste que importa: o veredito medido do controle, sem rede.
    Ciprofloxacino tem neuropatia em `boxed_warning` e NAO e indicado para a
    condicao -- alerta FORTE."""
    from mods.portao import bula
    local = massa.Acervo(acervo["indice"])
    resultado = bula.buscar_bulas("CIPROFLOXACIN", obter_json=local.obter_json)
    dirigidos = {
        f: bula.achados_dirigidos("CIPROFLOXACIN", resultado["campo_do_nome"],
                                  ["neuropathy"], f,
                                  obter_json=local.obter_json)
        for f in ("forte", "fraca")}
    lista = bula.achados(resultado["bulas"], ["neuropathy"])
    veredito, valor = bula.veredito(resultado, lista, dirigidos)
    assert veredito == "reprova"


def test_a_indicacao_desconta_no_acervo_local(acervo):
    """Capsaicina bate `paresthesia` em secao forte E e indicada para dor
    neuropatica: a indicacao cancela o alerta. E a correcao 3 do criterio,
    medida no controle de 12, e tem de valer igual no acervo local."""
    from mods.portao import bula
    local = massa.Acervo(acervo["indice"])
    achados = bula.achados_de_indicacao(
        "CAPSAICIN", "generic_name", ["neuropathic pain"],
        obter_json=local.obter_json)
    assert achados and achados[0]["forca"] == "indicacao"


def test_a_consulta_dirigida_conta_o_acervo_inteiro(acervo):
    """`meta.results.total` da dirigida e o peso da evidencia: quantas bulas do
    principio ativo casam no ACERVO, nao na janela."""
    from mods.portao import bula
    local = massa.Acervo(acervo["indice"])
    r = bula.achados_dirigidos("CIPROFLOXACIN", "generic_name",
                               ["neuropathy"], "forte",
                               obter_json=local.obter_json)
    assert r["total"] == 1


def test_termo_ausente_da_dirigida_total_zero(acervo):
    from mods.portao import bula
    local = massa.Acervo(acervo["indice"])
    r = bula.achados_dirigidos("PARACETAMOL", "generic_name",
                               ["neuropathy"], "forte",
                               obter_json=local.obter_json)
    assert r["total"] == 0 and r["achados"] == []


def test_a_busca_por_termo_e_por_palavra_inteira(acervo):
    """A mesma regra do criterio na rede: sem isso "neuritis" casaria dentro de
    outras palavras e o falso positivo nao seria nem conferivel."""
    from mods.portao import bula
    local = massa.Acervo(acervo["indice"])
    r = bula.achados_dirigidos("CIPROFLOXACIN", "generic_name",
                               ["neuropath"], "forte",
                               obter_json=local.obter_json)
    assert r["total"] == 0  # 'neuropath' nao e palavra inteira em 'neuropathy'


def test_url_desconhecida_e_erro_nomeado(acervo):
    """O acervo responde as URLs que o `bula.py` monta. Qualquer outra forma e
    erro que NOMEIA a URL: devolver vazio em silencio faria uma consulta nova
    sair `sem alerta` sem ninguem notar."""
    local = massa.Acervo(acervo["indice"])
    with pytest.raises(ValueError, match="nao sei responder"):
        local.obter_json("https://api.fda.gov/drug/event.json?limit=1")


def test_acervo_de_arquivo_ausente_e_erro_com_remedio(tmp_path):
    with pytest.raises(FileNotFoundError, match="openfda_massa"):
        massa.Acervo(tmp_path / "nao_existe.sqlite")


def test_a_janela_sai_ordenada_por_data_decrescente(tmp_path):
    """`effective_time:desc` e a ordem medida: sem ela, as 50 bulas que a busca
    devolve sao arbitrarias e a unica que indica o farmaco fica de fora."""
    pasta = tmp_path / "bruto"
    pasta.mkdir()
    antiga = _rotulo("X", warnings="a")
    antiga["effective_time"] = "20100101"
    nova = _rotulo("X", warnings="b")
    nova["effective_time"] = "20260101"
    (pasta / "p.json.zip").write_bytes(_zip_de([antiga, nova]))
    indice = tmp_path / "i.sqlite"
    massa.indexar(pasta, indice, export_date="x")
    from mods.portao import bula
    local = massa.Acervo(indice)
    r = bula.buscar_bulas("X", obter_json=local.obter_json, limite=2)
    assert [b["effective_time"] for b in r["bulas"]] == ["20260101", "20100101"]


def test_o_limite_da_janela_e_respeitado(tmp_path):
    pasta = tmp_path / "bruto"
    pasta.mkdir()
    (pasta / "p.json.zip").write_bytes(
        _zip_de([_rotulo("Y", warnings=f"w{i}") for i in range(10)]))
    indice = tmp_path / "i.sqlite"
    massa.indexar(pasta, indice, export_date="x")
    from mods.portao import bula
    local = massa.Acervo(indice)
    r = bula.buscar_bulas("Y", obter_json=local.obter_json, limite=3)
    assert len(r["bulas"]) == 3 and r["total"] == 10
