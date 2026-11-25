# -*- coding: utf-8 -*-
import json

import pytest

import mods.coletar_evidencias as cli
from mods.coletar_evidencias import main
from mods.evidencias.cache import (
    Cache, FaltaNoCache, VersaoIndeterminada, abrir_ao_vivo,
    abrir_para_repeticao, consultar_versoes)
from mods.evidencias.chembl import ListaTruncada, _sem_truncar
from mods.evidencias.manifesto import ler_manifesto, verificar_manifesto

URL = "https://www.ebi.ac.uk/chembl/api/data/target/CHEMBL1.json"


def _buscador(resposta, chamadas):
    def buscar(url, corpo=None):
        chamadas.append((url, corpo))
        if isinstance(resposta, Exception):
            raise resposta
        return resposta
    return buscar


def test_repeticao_recusa_versao_divergente():
    anterior = {"versoes": {"open_targets": "24.09", "checkpoint_md5": "abc"}}
    atual = {"versoes": {"open_targets": "25.03", "checkpoint_md5": "abc"}}
    assert "open_targets" in " ".join(verificar_manifesto(anterior, atual))
    assert verificar_manifesto(anterior, anterior) == []
    atual2 = {"versoes": {"open_targets": "24.09", "checkpoint_md5": "zzz"}}
    assert "checkpoint_md5" in " ".join(verificar_manifesto(anterior, atual2))


def test_repeticao_recusa_parametro_que_muda_a_saida():
    anterior = {"parametros": {"corte": 0.1, "limite": None}}
    assert verificar_manifesto(anterior, {"parametros": {"corte": 0.1, "limite": None}}) == []
    assert "limite" in " ".join(verificar_manifesto(anterior, {"parametros": {"corte": 0.1, "limite": 5}}))


def test_queda_de_rede_nao_deixa_saida_parcial(tmp_path, rotulos_falsos, monkeypatch):
    # Falhar antes de escrever, nunca a meio arquivo.
    def cai(url, corpo=None):
        raise OSError("rede fora do ar")
    monkeypatch.setattr("mods.evidencias.cache._buscar_padrao", cai)
    rot = tmp_path / "labels.tsv"
    rot.write_text("target_idx\ttid\tchembl_id\tpref_name\ttarget_type\torganism\t"
                   "active_molecules\tis_ckd\n0\t1\tCHEMBL222\tNET\tSINGLE PROTEIN\t"
                   "Homo sapiens\t5\tTrue\n", encoding="utf-8")
    saida = tmp_path / "evidencias.tsv"
    codigo = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
                   "--alvos", "tests/dados/alvos_minimo.tsv",
                   "--sintomas", "tests/dados/sintomas.txt",
                   "--saida", str(saida), "--rotulos", str(rot),
                   "--alvos-enriquecidos", str(tmp_path / "enr.tsv")])
    assert codigo != 0
    assert not saida.exists()
    assert not (tmp_path / "manifesto.json").exists()


def test_cache_serve_a_segunda_chamada_sem_rede(tmp_path):
    chamadas = []
    cache = Cache(tmp_path, _buscador({"a": [1]}, chamadas))
    primeira = cache.obter_json("https://exemplo/teste")
    segunda = Cache(tmp_path, _buscador(OSError("sem rede"), chamadas)).obter_json(
        "https://exemplo/teste")
    assert primeira == segunda == {"a": [1]}
    assert len(chamadas) == 1


def test_chave_distingue_corpo_e_url(tmp_path):
    chamadas = []
    cache = Cache(tmp_path, _buscador({"x": 1}, chamadas))
    cache.obter_json("https://e/g", {"q": 1})
    cache.obter_json("https://e/g", {"q": 2})
    cache.obter_json("https://e/h", {"q": 1})
    cache.obter_json("https://e/g", {"q": 1})
    assert len(chamadas) == 3


@pytest.mark.parametrize("resposta", [
    {"errors": [{"message": "x"}]},
    {"error_message": "boom"},
    {"targets": [{}], "page_meta": {"total_count": 5}},   # pagina truncada
    None])
def test_resposta_ruim_nao_e_gravada(tmp_path, resposta):
    chamadas = []
    cache = Cache(tmp_path, _buscador(resposta, chamadas))
    assert cache.obter_json(URL) == resposta   # devolvida ao chamador...
    cache.obter_json(URL)
    assert len(chamadas) == 2                  # ...mas nunca servida do disco
    assert not list(tmp_path.glob("*.json"))


def test_excecao_propaga_e_nao_grava(tmp_path):
    cache = Cache(tmp_path, _buscador(OSError("timeout"), []))
    with pytest.raises(OSError):
        cache.obter_json(URL)
    assert not list(tmp_path.iterdir())


def test_pagina_truncada_segue_erro_de_item(tmp_path):
    resposta = {"targets": [{}], "page_meta": {"total_count": 5}}
    cache = Cache(tmp_path, _buscador(resposta, []))
    with pytest.raises(ListaTruncada):
        _sem_truncar(cache.obter_json(URL), "targets", "teste")


def test_somente_leitura_falta_e_erro_e_nao_busca(tmp_path):
    chamadas = []
    cache = Cache(tmp_path, _buscador({"a": 1}, chamadas), somente_leitura=True)
    with pytest.raises(FaltaNoCache):
        cache.obter_json(URL)
    assert chamadas == []


def test_entrada_corrompida_falha_alto(tmp_path):
    cache = Cache(tmp_path, _buscador({"a": 1}, []))
    cache.obter_json(URL)
    arquivo = next(tmp_path.glob("*.json"))
    arquivo.write_text("{ truncado", encoding="utf-8")
    with pytest.raises(RuntimeError, match="corrompida"):
        cache.obter_json(URL)


def _servicos(ot_mes=9, chembl="ChEMBL_35"):
    def buscar(url, corpo=None):
        if "opentargets" in url:
            return {"data": {"meta": {"dataVersion": {"year": 24, "month": ot_mes,
                                                      "iteration": None}}}}
        if url.endswith("status.json"):
            return {"chembl_db_version": chembl}
        if url.endswith("database/version"):
            return 93
        return {"x": 1}
    return buscar


def test_versoes_vem_da_api_e_separam_as_pastas(tmp_path):
    assert consultar_versoes(_servicos()) == {
        "open_targets": "24.09", "chembl": "ChEMBL_35", "reactome": "93"}
    c1, v1 = abrir_ao_vivo(tmp_path, _servicos(9))
    c2, v2 = abrir_ao_vivo(tmp_path, _servicos(12))
    assert c1.diretorio != c2.diretorio     # release nova nunca le o cache velho
    c1.obter_json("https://e/x")
    chamadas = []
    c2._buscar = _buscador({"novo": 1}, chamadas)
    assert c2.obter_json("https://e/x") == {"novo": 1} and len(chamadas) == 1
    r, versoes = abrir_para_repeticao(tmp_path, c1.diretorio.name)
    assert versoes == v1 and r.obter_json("https://e/x") == {"x": 1}


def test_versao_indeterminada_aborta(tmp_path):
    def quebrado(url, corpo=None):
        return {"inesperado": True}
    with pytest.raises(VersaoIndeterminada):
        consultar_versoes(quebrado)


def _rodar_enriquecido(tmp_path, nome, extra=()):
    rot = tmp_path / "labels.tsv"
    rot.write_text("target_idx\ttid\tchembl_id\tpref_name\ttarget_type\torganism\t"
                   "active_molecules\tis_ckd\n0\t1\tCHEMBL222\tNET\tSINGLE PROTEIN\t"
                   "Homo sapiens\t5\tTrue\n", encoding="utf-8")
    pasta = tmp_path / nome
    codigo = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
                   "--alvos", "tests/dados/alvos_minimo.tsv",
                   "--sintomas", "tests/dados/sintomas.txt",
                   "--saida", str(pasta / "e.tsv"), "--rotulos", str(rot),
                   "--alvos-enriquecidos", str(pasta / "enr.tsv"),
                   "--cache", str(tmp_path / "cache"), *extra])
    return codigo, pasta


def test_manifesto_grava_versoes_e_repeticao_funciona_sem_rede(tmp_path, monkeypatch):
    monkeypatch.setattr("mods.evidencias.cache._buscar_padrao", _servicos())
    monkeypatch.setattr(cli, "_enriquecer",
                        lambda *a: [{"alvo_chembl_id": "CHEMBL222",
                                                       "origem": "usuario"}])
    codigo, pasta = _rodar_enriquecido(tmp_path, "um")
    assert codigo == 0
    m = ler_manifesto(pasta / "manifesto.json")
    assert m["versoes"]["open_targets"] == "24.09" and m["versoes"]["chembl"] == "ChEMBL_35"
    assert m["parametros"]["cache_pasta"].startswith("open_targets-24.09")

    # Repeticao: a rede agora falha; mesmo assim roda, e nao toca o manifesto lido.
    def cai(url, corpo=None):
        raise OSError("fora do ar")
    monkeypatch.setattr("mods.evidencias.cache._buscar_padrao", cai)
    original = (pasta / "manifesto.json").read_text()
    codigo, pasta2 = _rodar_enriquecido(
        tmp_path, "dois", ["--manifesto", str(pasta / "manifesto.json")])
    assert codigo == 0
    assert (pasta / "manifesto.json").read_text() == original
    assert (pasta2 / "e.tsv").read_text() == (pasta / "e.tsv").read_text()


def test_repeticao_recusa_cache_de_outra_versao(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr("mods.evidencias.cache._buscar_padrao", _servicos())
    monkeypatch.setattr(cli, "_enriquecer", lambda *a: [])
    codigo, pasta = _rodar_enriquecido(tmp_path, "um")
    assert codigo == 0
    caminho = pasta / "manifesto.json"
    m = json.loads(caminho.read_text())
    m["versoes"]["open_targets"] = "23.12"
    caminho.write_text(json.dumps(m))
    codigo, pasta2 = _rodar_enriquecido(tmp_path, "dois", ["--manifesto", str(caminho)])
    assert codigo == 2
    assert "open_targets" in capsys.readouterr().err
    assert not (pasta2 / "e.tsv").exists()


def test_verificar_nomeia_o_arquivo_de_entrada_que_mudou():
    anterior = {"versoes": {"ligantes_md5": "a", "alvos_md5": "b", "sintomas_md5": "c"}}
    atual = {"versoes": {"ligantes_md5": "a", "alvos_md5": "b", "sintomas_md5": "X"}}
    msg = " ".join(verificar_manifesto(anterior, atual))
    assert "sintomas" in msg and "ligantes" not in msg and "alvos" not in msg


def _servico_falso():
    """Os tres servicos, no formato real, para um alvo (SCN9A)."""
    scn9a = {"target_chembl_id": "CHEMBL4296", "pref_name": "SCN9A",
             "organism": "Homo sapiens", "target_type": "SINGLE PROTEIN"}
    chamadas = []

    def servico(url, corpo=None):
        chamadas.append(url)
        if corpo and "Versao" in corpo["query"] or url.endswith(("status.json", "database/version")):
            return _servicos()(url, corpo)
        if corpo is not None:
            return {"data": {"disease": {"associatedTargets": {"count": 1, "rows": [
                {"score": 0.6, "datatypeScores": [],
                 "target": {"id": "E", "approvedSymbol": "SCN9A",
                            "proteinIds": [{"id": "Q15858", "source": "uniprot_swissprot"}]}}]}}}}
        if url.startswith("https://reactome.org"):
            if "mapping" in url:
                return [{"stId": "R-X", "displayName": "x"}]
            return [{"databaseName": "UniProt", "identifier": "Q15858"}]
        if "target/CHEMBL4296.json" in url:
            return {**scn9a, "target_components": [{"accession": "Q15858"}]}
        if "accession=Q15858" in url:
            return {"targets": [scn9a], "page_meta": {"total_count": 1}}
        return {"targets": [], "page_meta": {"total_count": 0}}
    return servico, chamadas


def _entrada_real(tmp_path):
    (tmp_path / "alvos.tsv").write_text("alvo_chembl_id\nCHEMBL4296\n", encoding="utf-8")
    (tmp_path / "sint.txt").write_text("dor\n", encoding="utf-8")
    (tmp_path / "lig.tsv").write_text(
        open("tests/dados/ligantes_minimo.tsv", encoding="utf-8").read(), encoding="utf-8")
    (tmp_path / "labels.tsv").write_text(
        "target_idx\ttid\tchembl_id\tpref_name\ttarget_type\torganism\tactive_molecules\tis_ckd\n"
        "0\t1\tCHEMBL4296\tSCN9A\tSINGLE PROTEIN\tHomo sapiens\t50\tTrue\n", encoding="utf-8")


def _cli_real(tmp_path, nome, extra=()):
    pasta = tmp_path / nome
    codigo = main(["--ligantes", str(tmp_path / "lig.tsv"), "--alvos", str(tmp_path / "alvos.tsv"),
                   "--sintomas", str(tmp_path / "sint.txt"), "--saida", str(pasta / "e.tsv"),
                   "--rotulos", str(tmp_path / "labels.tsv"), "--efo", "dor=EFO_1",
                   "--alvos-enriquecidos", str(pasta / "enr.tsv"),
                   "--cache", str(tmp_path / "cache"), *extra])
    return codigo, pasta


def test_enriquecimento_real_repetido_so_do_cache_e_entrada_alterada_recusa(
        tmp_path, monkeypatch, capsys):
    servico, chamadas = _servico_falso()
    monkeypatch.setattr("mods.evidencias.cache._buscar_padrao", servico)
    _entrada_real(tmp_path)
    codigo, vivo = _cli_real(tmp_path, "vivo")
    assert codigo == 0
    enr = (vivo / "enr.tsv").read_text(encoding="utf-8")
    assert "CHEMBL4296" in enr and "\t0.6\t" in enr and len(chamadas) > 3
    m = ler_manifesto(vivo / "manifesto.json")
    assert {"ligantes_md5", "alvos_md5", "sintomas_md5"} <= m["versoes"].keys()

    def proibido(url, corpo=None):
        raise AssertionError(f"repeticao foi a rede: {url}")
    monkeypatch.setattr("mods.evidencias.cache._buscar_padrao", proibido)
    ref = ["--manifesto", str(vivo / "manifesto.json")]
    codigo, rep = _cli_real(tmp_path, "rep", ref)     # socket bloqueado tambem
    assert codigo == 0
    assert (rep / "enr.tsv").read_text(encoding="utf-8") == enr

    # --limite diferente: recusa.
    codigo, _ = _cli_real(tmp_path, "rep2", [*ref, "--limite", "1"])
    assert codigo == 2 and "limite" in capsys.readouterr().err

    # Entrada alterada: recusa nomeando o arquivo.
    (tmp_path / "sint.txt").write_text("dor\nfadiga\n", encoding="utf-8")
    codigo, rep3 = _cli_real(tmp_path, "rep3", ref)
    err = capsys.readouterr().err
    assert codigo == 2 and "sintomas_md5" in err and "ligantes_md5" not in err
    assert not (rep3 / "enr.tsv").exists()


def test_manifesto_sem_hash_de_entrada_recusa_a_repeticao(tmp_path, monkeypatch, capsys):
    servico, _ = _servico_falso()
    monkeypatch.setattr("mods.evidencias.cache._buscar_padrao", servico)
    _entrada_real(tmp_path)
    codigo, vivo = _cli_real(tmp_path, "vivo")
    assert codigo == 0
    caminho = vivo / "manifesto.json"
    m = json.loads(caminho.read_text())
    del m["versoes"]["ligantes_md5"]
    caminho.write_text(json.dumps(m))
    capsys.readouterr()
    codigo, rep = _cli_real(tmp_path, "rep", ["--manifesto", str(caminho)])
    err = capsys.readouterr().err
    assert codigo == 2
    assert "anterior ao registro de MD5" in err and "ligantes_md5" in err
    assert "alvos_md5" not in err
    assert not (rep / "e.tsv").exists() and not (rep / "manifesto.json").exists()
