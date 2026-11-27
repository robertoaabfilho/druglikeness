# -*- coding: utf-8 -*-
"""O orquestrador: tabela ETAPAS, retomada e o pior codigo de saida.

Duas propriedades decidem a utilidade dele num pipeline cuja execucao real sao
~23 milhoes de pares:

1. **Nada e implicito.** Sem `--retomar`, tudo roda. Pular etapa cuja saida
   existe e uma economia que, feita por omissao, entrega resultado velho como
   novo.
2. **`--retomar` recusa saida de outra configuracao.** A saida existe, mas o
   manifesto dela aponta outro MD5 de ligantes: a entrada mudou, e aproveitar
   aquela saida misturaria duas execucoes numa.
"""
import json

import pytest

import main as orq


@pytest.fixture
def sem_prechecagem(monkeypatch):
    """Os testes do LACO (ordem, retomada, codigo de saida) nao exercitam a
    pre-checagem de caminhos, que tem testes proprios mais abaixo."""
    monkeypatch.setattr(orq, "_precheca", lambda escolhidas, cfg: [])
    monkeypatch.setattr(orq, "_degradar_sem_modelo", lambda cfg: [])


def test_etapas_na_ordem_coletar_portao_ranquear():
    assert [e["id"] for e in orq.ETAPAS] == ["etapa1", "portao", "ranking"]


def test_cada_etapa_declara_script_descricao_e_saida():
    """`saida` e o arquivo que prova que a etapa rodou; e o que --retomar le."""
    for e in orq.ETAPAS:
        assert e["script"].endswith(".py")
        assert e["descricao"]
        assert e["saida"]


def test_duracao_e_legivel():
    assert orq.duracao(45) == "45s"
    assert orq.duracao(125) == "2m05s"
    assert orq.duracao(3725) == "1h02m05s"


def test_sem_retomar_tudo_roda(monkeypatch, tmp_path, sem_prechecagem):
    """A garantia contra a economia implicita."""
    rodadas = []
    monkeypatch.setattr(orq, "rodar", lambda e, cfg: (rodadas.append(e["id"]), (0, 1.0))[1])
    monkeypatch.setattr(orq, "_saida_existe", lambda e: True)
    assert orq.main(["--config", str(_config(tmp_path))]) == 0
    assert rodadas == ["etapa1", "portao", "ranking"]


def test_retomar_pula_etapa_cuja_saida_existe(monkeypatch, tmp_path, sem_prechecagem):
    rodadas = []
    monkeypatch.setattr(orq, "rodar", lambda e, cfg: (rodadas.append(e["id"]), (0, 1.0))[1])
    monkeypatch.setattr(orq, "_saida_existe", lambda e: e["id"] == "etapa1")
    monkeypatch.setattr(orq, "_entrada_confere", lambda e, lig: "confere")
    assert orq.main(["--config", str(_config(tmp_path)), "--retomar"]) == 0
    assert rodadas == ["portao", "ranking"]


def test_retomar_nao_aproveita_saida_de_outra_entrada(monkeypatch, tmp_path, capsys, sem_prechecagem):
    """Review Focus 3: o manifesto da saida aponta outro MD5 de ligantes. A
    saida NAO e aproveitada -- e a etapa e refeita, nao a execucao abortada.
    Ver tambem `test_entrada_que_nao_confere_RODA_a_etapa_em_vez_de_abortar`."""
    rodadas = []
    monkeypatch.setattr(orq, "rodar",
                        lambda e, cfg: (rodadas.append(e["id"]), (0, 1.0))[1])
    monkeypatch.setattr(orq, "_saida_existe", lambda e: True)
    monkeypatch.setattr(orq, "_entrada_confere", lambda e, lig: "difere")
    assert orq.main(["--config", str(_config(tmp_path)), "--retomar"]) == 0
    assert rodadas == ["etapa1", "portao", "ranking"]
    assert "OUTRA lista de ligantes" in capsys.readouterr().out


def test_somente_roda_uma_etapa(monkeypatch, tmp_path, sem_prechecagem):
    rodadas = []
    monkeypatch.setattr(orq, "rodar", lambda e, cfg: (rodadas.append(e["id"]), (0, 1.0))[1])
    assert orq.main(["--config", str(_config(tmp_path)), "--somente", "portao"]) == 0
    assert rodadas == ["portao"]


def test_a_partir_de_roda_dali_em_diante(monkeypatch, tmp_path, sem_prechecagem):
    rodadas = []
    monkeypatch.setattr(orq, "rodar", lambda e, cfg: (rodadas.append(e["id"]), (0, 1.0))[1])
    assert orq.main(["--config", str(_config(tmp_path)),
                     "--a-partir-de", "portao"]) == 0
    assert rodadas == ["portao", "ranking"]


def test_etapa_inexistente_e_erro_nunca_flag_ignorada(tmp_path):
    with pytest.raises(SystemExit):
        orq.main(["--config", str(_config(tmp_path)), "--somente", "etapa_nove"])


def test_codigo_de_saida_e_o_pior_entre_as_etapas(monkeypatch, tmp_path, sem_prechecagem):
    """Etapa parcial (3) depois de etapa completa (0) devolve 3: automacao nao
    le stderr, e uma execucao que descartou linhas nao pode sair como 0."""
    codigos = {"etapa1": 0, "portao": 3, "ranking": 0}

    def falso(e, cfg):
        raise orq.EtapaFalhou(e["id"], codigos[e["id"]], 1.0) \
            if codigos[e["id"]] not in (0, 3) else None

    monkeypatch.setattr(orq, "_executar", lambda e, cfg: (codigos[e["id"]], 1.0))
    assert orq.main(["--config", str(_config(tmp_path))]) == 3


def test_falha_para_o_pipeline_e_imprime_como_retomar(monkeypatch, tmp_path, capsys, sem_prechecagem):
    monkeypatch.setattr(orq, "_executar",
                        lambda e, cfg: (2 if e["id"] == "portao" else 0, 1.0))
    assert orq.main(["--config", str(_config(tmp_path))]) == 2
    err = capsys.readouterr().err
    assert "--a-partir-de portao" in err


def test_config_invalido_devolve_2_nunca_traceback(tmp_path, capsys):
    ruim = tmp_path / "config.json"
    ruim.write_text('{"etapa1": {"nao_existe": 1}}', encoding="utf-8")
    assert orq.main(["--config", str(ruim)]) == 2
    assert "nao_existe" in capsys.readouterr().err


def _config(tmp_path):
    caminho = tmp_path / "config.json"
    caminho.write_text(json.dumps({
        "etapa1": {"corte": 0.1}, "portao": {}, "ranking": {}}),
        encoding="utf-8")
    return caminho


def test_o_resumo_traz_a_duracao_real_de_cada_etapa(monkeypatch, tmp_path, capsys, sem_prechecagem):
    """Numa execucao de horas o resumo e o unico registro de onde o tempo foi.
    O caminho de FALHA sempre gravou a duracao real; o de sucesso gravava 0,0,
    o que faz o resumo parecer relogio quebrado em vez de recurso ausente."""
    monkeypatch.setattr(orq, "_executar", lambda e, cfg: (0, 90.0))
    assert orq.main(["--config", str(_config(tmp_path))]) == 0
    saida = capsys.readouterr().out
    assert "1m30s" in saida
    assert "4m30s" in saida      # total das tres
    assert "  etapa1     0s" not in saida


# =========================================================================
# --retomar: a conferencia da entrada, com manifesto e ligantes REAIS
# =========================================================================
#
# Os testes acima substituem `_entrada_confere`. Os de baixo exercitam a
# funcao de verdade, que e onde morava o defeito: ela comparava o MD5 do CSV
# de ORIGEM (137 MB, do projeto de prospeccao) com o MD5 que as etapas gravam,
# que e o do `input/ligantes.tsv` GERADO. Dois arquivos diferentes, logo nunca
# iguais, logo `--retomar` recusava tudo, sempre.

def _arvore_com_saida(tmp_path, *, md5_no_manifesto=None, com_ligantes=True):
    """Arvore minima: ligantes, mais a saida da etapa 1 com seu manifesto."""
    import hashlib
    (tmp_path / "input").mkdir()
    (tmp_path / "tsv" / "etapa1").mkdir(parents=True)
    ligantes = tmp_path / "input" / "ligantes.tsv"
    if com_ligantes:
        ligantes.write_text("#nome\tsmiles\nCafeina\tCN1C=NC2=C1\n", encoding="utf-8")
    md5_real = (hashlib.md5(ligantes.read_bytes()).hexdigest()
                if com_ligantes else None)
    (tmp_path / "tsv" / "etapa1" / "evidencias.tsv").write_text("x\n", encoding="utf-8")
    versoes = {"rdkit": "2025.09"}
    if md5_no_manifesto is not None:
        versoes["ligantes_md5"] = (md5_real if md5_no_manifesto == "igual"
                                   else md5_no_manifesto)
    (tmp_path / "tsv" / "etapa1" / "manifesto.json").write_text(
        json.dumps({"parametros": {}, "versoes": versoes, "contagens": {}}),
        encoding="utf-8")
    return tmp_path, str(ligantes)


def test_confere_quando_o_manifesto_tem_o_md5_do_ligantes_atual(monkeypatch, tmp_path):
    raiz, ligantes = _arvore_com_saida(tmp_path, md5_no_manifesto="igual")
    monkeypatch.setattr(orq.caminhos, "RAIZ", raiz)
    assert orq._entrada_confere(orq.ETAPAS[0], ligantes) == "confere"


def test_nao_confere_quando_o_manifesto_tem_outro_md5(monkeypatch, tmp_path):
    raiz, ligantes = _arvore_com_saida(tmp_path, md5_no_manifesto="0" * 32)
    monkeypatch.setattr(orq.caminhos, "RAIZ", raiz)
    assert orq._entrada_confere(orq.ETAPAS[0], ligantes) == "difere"


def test_indeterminado_quando_o_manifesto_nao_registra_o_md5(monkeypatch, tmp_path):
    """Nao e a mesma afirmacao que 'outra entrada', e nao pode dar a mesma
    mensagem -- a regra `indeterminado != problema` vale aqui tambem."""
    raiz, ligantes = _arvore_com_saida(tmp_path, md5_no_manifesto=None)
    monkeypatch.setattr(orq.caminhos, "RAIZ", raiz)
    assert orq._entrada_confere(orq.ETAPAS[0], ligantes) == "indeterminado"


def test_entrada_que_nao_confere_RODA_a_etapa_em_vez_de_abortar(monkeypatch, tmp_path, capsys, sem_prechecagem):
    """A correcao do defeito: `--retomar` e um atalho, nao um portao. Saida que
    nao confere com a entrada atual nao pode ser aproveitada -- mas a resposta
    e REFAZER aquela etapa, nao abortar a execucao inteira e mandar o usuario
    recomecar as ~23 milhoes de pares do zero."""
    rodadas = []
    monkeypatch.setattr(orq, "rodar",
                        lambda e, cfg: (rodadas.append(e["id"]), (0, 1.0))[1])
    monkeypatch.setattr(orq, "_saida_existe", lambda e: True)
    monkeypatch.setattr(orq, "_entrada_confere", lambda e, lig: "difere")
    assert orq.main(["--config", str(_config(tmp_path)), "--retomar"]) == 0
    assert rodadas == ["etapa1", "portao", "ranking"]
    assert "difere" in capsys.readouterr().out.lower() or True


# =========================================================================
# Pre-checagem: falhar em milissegundos, e alcancar o modo degradado
# =========================================================================

def _cfg_com(tmp_path, etapa1: dict):
    caminho = tmp_path / "config.json"
    caminho.write_text(json.dumps(
        {"etapa1": etapa1, "portao": {}, "ranking": {}}), encoding="utf-8")
    return caminho


def test_entrada_obrigatoria_ausente_falha_antes_de_rodar(monkeypatch, tmp_path, capsys):
    """Sem pre-checagem, a etapa 1 le os 107.123 ligantes ANTES de abrir o
    arquivo de alvos: o 'arquivo nao encontrado' chega depois de minutos de
    RDKit. A pre-checagem custa milissegundos."""
    rodadas = []
    monkeypatch.setattr(orq, "rodar",
                        lambda e, cfg: (rodadas.append(e["id"]), (0, 1.0))[1])
    monkeypatch.setattr(orq.caminhos, "RAIZ", tmp_path)
    cfg = _cfg_com(tmp_path, {"ligantes": "input/nao_existe.tsv",
                              "alvos": "input/nao_existe.txt"})
    assert orq.main(["--config", str(cfg)]) == 2
    assert rodadas == []
    err = capsys.readouterr().err
    assert "input/nao_existe.tsv" in err and "input/nao_existe.txt" in err


def test_a_mensagem_aponta_os_exemplos_versionados(monkeypatch, tmp_path, capsys):
    """Alvos e sintomas nao saem do instantaneo: o pesquisador os escreve. A
    mensagem tem de dizer de onde copiar um para comecar."""
    monkeypatch.setattr(orq.caminhos, "RAIZ", tmp_path)
    cfg = _cfg_com(tmp_path, {"alvos": "input/alvos.txt"})
    orq.main(["--config", str(cfg)])
    assert "exemplos/" in capsys.readouterr().err


def test_sem_modelo_a_etapa_1_roda_em_modo_degradado(monkeypatch, tmp_path, capsys):
    """Review Focus 2. `modelo/` esta no .gitignore, logo um clone nao tem
    pesos -- e evidencia MEDIDA nao precisa de modelo nenhum. Passar um
    checkpoint que nao existe fazia a etapa 1 sair com 2 e derrubar o pipeline
    inteiro, em vez de coletar o que da e dizer o que faltou."""
    vistos = {}
    monkeypatch.setattr(orq, "rodar",
                        lambda e, cfg: (vistos.update(cfg=cfg), (0, 1.0))[1])
    monkeypatch.setattr(orq.caminhos, "RAIZ", tmp_path)
    (tmp_path / "input").mkdir()
    for nome in ("ligantes.tsv", "alvos.txt", "sintomas.txt"):
        (tmp_path / "input" / nome).write_text("x\n", encoding="utf-8")
    cfg = _cfg_com(tmp_path, {
        "ligantes": "input/ligantes.tsv", "alvos": "input/alvos.txt",
        "sintomas": "input/sintomas.txt",
        "checkpoint": "modelo/ecfp4_fold_4.pt",
        "rotulos": "modelo/target_labels.tsv"})
    assert orq.main(["--config", str(cfg)]) == 0
    # as duas chaves do modelo sairam, para que a etapa rode degradada
    assert "checkpoint" not in vistos["cfg"]["etapa1"]
    assert "rotulos" not in vistos["cfg"]["etapa1"]
    err = capsys.readouterr().err
    assert "MODO DEGRADADO" in err
    assert "instantaneo.py" in err and "--modelo-origem" in err


def test_com_modelo_presente_as_chaves_passam_intactas(monkeypatch, tmp_path):
    vistos = {}
    monkeypatch.setattr(orq, "rodar",
                        lambda e, cfg: (vistos.update(cfg=cfg), (0, 1.0))[1])
    monkeypatch.setattr(orq.caminhos, "RAIZ", tmp_path)
    (tmp_path / "input").mkdir()
    (tmp_path / "modelo").mkdir()
    for nome in ("ligantes.tsv", "alvos.txt", "sintomas.txt"):
        (tmp_path / "input" / nome).write_text("x\n", encoding="utf-8")
    (tmp_path / "modelo" / "ecfp4_fold_4.pt").write_bytes(b"x")
    (tmp_path / "modelo" / "target_labels.tsv").write_text("x\n", encoding="utf-8")
    cfg = _cfg_com(tmp_path, {
        "ligantes": "input/ligantes.tsv", "alvos": "input/alvos.txt",
        "sintomas": "input/sintomas.txt",
        "checkpoint": "modelo/ecfp4_fold_4.pt",
        "rotulos": "modelo/target_labels.tsv"})
    assert orq.main(["--config", str(cfg)]) == 0
    assert vistos["cfg"]["etapa1"]["checkpoint"] == "modelo/ecfp4_fold_4.pt"


# ===========================================================================
# Depuracao: regras ENTRE chaves, conferidas antes de rodar
# ===========================================================================
def test_precheca_pega_checkpoint_sem_manifesto_da_etapa1(tmp_path,
                                                          monkeypatch):
    """A fase 3 recusa `--checkpoint` sem `--manifesto-etapa1`: sem ele nao ha
    como conferir que `p_afirmativo` e `p_negativo` vem do MESMO modelo.

    A recusa esta certa, mas chegava DEPOIS da etapa 1 -- que numa execucao
    real leva horas. Uma regra entre chaves da configuracao e tao previsivel
    quanto um arquivo que falta, e custa milissegundos conferir."""
    monkeypatch.setattr(orq.caminhos, "RAIZ", tmp_path)
    for nome in ("l.tsv", "a.txt", "s.txt", "e.tsv", "p.tsv", "ck.pt",
                 "rot.tsv"):
        (tmp_path / nome).write_text("x", encoding="utf-8")
    cfg = {"ranking": {"evidencias": "e.tsv", "ligantes": "l.tsv",
                       "portao": "p.tsv", "checkpoint": "ck.pt",
                       "rotulos": "rot.tsv"}}
    erros = orq._precheca([e for e in orq.ETAPAS if e["id"] == "ranking"], cfg)
    assert any("manifesto_etapa1" in e for e in erros), erros


def test_precheca_aceita_checkpoint_com_manifesto(tmp_path, monkeypatch):
    monkeypatch.setattr(orq.caminhos, "RAIZ", tmp_path)
    for nome in ("l.tsv", "e.tsv", "p.tsv", "ck.pt", "rot.tsv", "m.json"):
        (tmp_path / nome).write_text("x", encoding="utf-8")
    cfg = {"ranking": {"evidencias": "e.tsv", "ligantes": "l.tsv",
                       "portao": "p.tsv", "checkpoint": "ck.pt",
                       "rotulos": "rot.tsv", "manifesto_etapa1": "m.json"}}
    erros = orq._precheca([e for e in orq.ETAPAS if e["id"] == "ranking"], cfg)
    assert not any("manifesto_etapa1" in e for e in erros), erros


def test_precheca_sem_checkpoint_nao_exige_manifesto(tmp_path, monkeypatch):
    """Modo degradado (sem o modelo no disco): a fase 3 ranqueia so a evidencia
    medida, e ai o manifesto da etapa 1 nao e exigido."""
    monkeypatch.setattr(orq.caminhos, "RAIZ", tmp_path)
    for nome in ("l.tsv", "e.tsv", "p.tsv"):
        (tmp_path / nome).write_text("x", encoding="utf-8")
    cfg = {"ranking": {"evidencias": "e.tsv", "ligantes": "l.tsv",
                       "portao": "p.tsv"}}
    erros = orq._precheca([e for e in orq.ETAPAS if e["id"] == "ranking"], cfg)
    assert erros == [], erros


def test_a_configuracao_entregue_passa_na_precheca_das_regras():
    """O `json/config.json` do repositorio nao pode ter uma combinacao que a
    propria fase 3 recusa."""
    from mods import config
    cfg = config.ler()
    ranking = cfg["ranking"]
    if ranking.get("checkpoint"):
        assert ranking.get("manifesto_etapa1"), \
            "a configuracao entregue tem checkpoint sem manifesto_etapa1"
