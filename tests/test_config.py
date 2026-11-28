# -*- coding: utf-8 -*-
"""A configuracao: um lugar so para o parametro de cada etapa.

A regra que atravessa o arquivo: **nome que nao existe e ERRO**, nunca chave
ignorada em silencio. E a mesma regra das flags das CLIs, e existe pelo mesmo
motivo: uma politica que o pesquisador acredita estar valendo e nao esta e pior
que uma mensagem de erro.
"""
import json

import pytest

from mods import config


def _escrever(tmp_path, dados) -> str:
    caminho = tmp_path / "config.json"
    caminho.write_text(json.dumps(dados), encoding="utf-8")
    return str(caminho)


def test_le_o_config_versionado_do_projeto():
    """O arquivo que vem no repositorio tem de ser valido por construcao."""
    cfg = config.ler()
    assert set(cfg) >= {"etapa1", "portao", "ranking"}


def test_secao_inexistente_e_erro_nomeando_a_secao(tmp_path):
    caminho = _escrever(tmp_path, {"etapa_quatro": {}})
    with pytest.raises(config.ErroDeConfig) as e:
        config.ler(caminho)
    assert "etapa_quatro" in str(e.value)


def test_chave_inexistente_e_erro_nomeando_a_chave(tmp_path):
    caminho = _escrever(tmp_path, {"etapa1": {"corte_que_nao_existe": 1}})
    with pytest.raises(config.ErroDeConfig) as e:
        config.ler(caminho)
    assert "corte_que_nao_existe" in str(e.value)


def test_a_mensagem_de_erro_sugere_a_chave_parecida(tmp_path):
    """Erro de digitacao e o caso comum; dizer qual era a chave poupa a busca."""
    caminho = _escrever(tmp_path, {"etapa1": {"corte_associacao_": 0.4}})
    with pytest.raises(config.ErroDeConfig) as e:
        config.ler(caminho)
    assert "corte_associacao" in str(e.value)


def test_valor_de_tipo_errado_e_erro_nomeando_a_chave_e_o_tipo(tmp_path):
    """'corte': 'alto' nao vira 0."""
    caminho = _escrever(tmp_path, {"etapa1": {"corte": "alto"}})
    with pytest.raises(config.ErroDeConfig) as e:
        config.ler(caminho)
    msg = str(e.value)
    assert "corte" in msg and "alto" in msg


def test_argumentos_de_monta_as_flags_da_cli(tmp_path):
    cfg = {"etapa1": {"corte": 0.2, "corte_associacao": 0.5,
                      "sem_rede": True, "limite": None}}
    args = config.argumentos_de("etapa1", cfg)
    assert "--corte" in args and "0.2" in args
    assert "--corte-associacao" in args and "0.5" in args
    # booleano verdadeiro entra como flag sem valor
    assert "--sem-rede" in args
    assert args.count("--sem-rede") == 1
    # ausente nao entra: a CLI tem o seu proprio padrao
    assert "--limite" not in args


def test_booleano_falso_nao_entra_como_flag(tmp_path):
    args = config.argumentos_de("etapa1", {"etapa1": {"sem_rede": False}})
    assert "--sem-rede" not in args


def test_lista_entra_como_flag_repetida():
    """`--efo` e repetivel: um sintoma por ocorrencia."""
    args = config.argumentos_de(
        "etapa1", {"etapa1": {"efo": ["neuropathic pain=EFO_0005762",
                                      "peripheral neuropathy=EFO_0003100"]}})
    assert args.count("--efo") == 2


def test_toda_chave_conhecida_corresponde_a_uma_flag_real_da_cli():
    """A guarda que impede o esquema de descolar das CLIs: cada chave aceitada
    tem de existir como flag no script da etapa. Sao 16 + 17 + 13 = 46."""
    import pathlib
    import re
    for etapa, script in (("etapa1", "coletar_evidencias"),
                          ("portao", "aplicar_portao"),
                          ("ranking", "ranquear")):
        fonte = pathlib.Path("mods", f"{script}.py").read_text(encoding="utf-8")
        flags = set(re.findall(r'"--([a-z0-9-]+)"', fonte))
        for chave in config.CHAVES[etapa]:
            assert chave.replace("_", "-") in flags, f"{etapa}.{chave}"


def test_cpu_e_o_padrao_do_config_versionado():
    """GPU nao reproduz: 36 de 312 colunas diferem entre duas execucoes de GPU,
    0 de 312 na CPU. Nenhuma etapa pede GPU por omissao."""
    cfg = config.ler()
    for etapa in ("portao", "ranking"):
        assert cfg[etapa].get("gpu", False) is False
        assert "--gpu" not in config.argumentos_de(etapa, cfg)


# ===========================================================================
# T4 -- gravar: o painel escreve a configuracao, e as notas tem de sobreviver
# ===========================================================================
def test_gravar_preserva_as_notas(tmp_path):
    """Cada numero do config vem com um `_nota_*` ao lado carregando a medicao
    que o escolheu. Um `json.dump` do dicionario de valores apagaria todas, e
    o proximo a mexer no numero nao saberia por que ele e aquele."""
    caminho = _escrever(tmp_path, {"etapa1": {
        "_nota_corte": "0.1 e o regime de producao; 0 inverte o enriquecimento",
        "corte": 0.1}})
    config.gravar({"etapa1": {"corte": 0.2}}, caminho)
    depois = json.loads(open(caminho, encoding="utf-8").read())
    assert depois["etapa1"]["_nota_corte"].startswith("0.1 e o regime")
    assert depois["etapa1"]["corte"] == 0.2


def test_gravar_preserva_nota_de_secao(tmp_path):
    caminho = _escrever(tmp_path, {"_nota": "lido pelo main.py",
                                   "etapa1": {"corte": 0.1}})
    config.gravar({"etapa1": {"corte": 0.2}}, caminho)
    assert json.loads(open(caminho, encoding="utf-8").read())["_nota"]


def test_gravar_mantem_a_ordem_das_chaves(tmp_path):
    """A ordem do arquivo e a ordem em que o pesquisador le; reordenar por
    alfabeto embaralharia cada numero longe da sua nota."""
    original = {"portao": {"evidencias": "a", "_nota_limite": "n",
                           "limite": 1, "workers": 4}}
    caminho = _escrever(tmp_path, original)
    config.gravar({"portao": {"limite": 2}}, caminho)
    depois = json.loads(open(caminho, encoding="utf-8").read())
    assert list(depois["portao"]) == list(original["portao"])


def test_gravar_acrescenta_chave_nova_no_fim(tmp_path):
    caminho = _escrever(tmp_path, {"portao": {"limite": 1}})
    config.gravar({"portao": {"workers": 4}}, caminho)
    depois = json.loads(open(caminho, encoding="utf-8").read())
    assert list(depois["portao"]) == ["limite", "workers"]


def test_gravar_recusa_chave_desconhecida(tmp_path):
    caminho = _escrever(tmp_path, {"etapa1": {"corte": 0.1}})
    with pytest.raises(config.ErroDeConfig, match="cortte"):
        config.gravar({"etapa1": {"cortte": 0.2}}, caminho)


def test_gravar_recusa_secao_desconhecida(tmp_path):
    caminho = _escrever(tmp_path, {"etapa1": {"corte": 0.1}})
    with pytest.raises(config.ErroDeConfig, match="etapa9"):
        config.gravar({"etapa9": {"corte": 0.2}}, caminho)


def test_gravar_recusa_tipo_errado(tmp_path):
    """`limite: true` nao e um limite. A mesma regra da leitura."""
    caminho = _escrever(tmp_path, {"etapa1": {"limite": 10}})
    with pytest.raises(config.ErroDeConfig):
        config.gravar({"etapa1": {"limite": True}}, caminho)


def test_gravar_nada_deixa_o_arquivo_intacto_quando_recusa(tmp_path):
    """Validar ANTES de escrever: uma chave errada no formulario nao pode
    deixar o arquivo pela metade, porque e dele que as tres etapas leem."""
    caminho = _escrever(tmp_path, {"etapa1": {"corte": 0.1}})
    antes = open(caminho, encoding="utf-8").read()
    with pytest.raises(config.ErroDeConfig):
        config.gravar({"etapa1": {"cortte": 0.2}}, caminho)
    assert open(caminho, encoding="utf-8").read() == antes


def test_gravar_remove_chave_com_valor_none(tmp_path):
    """Apagar o campo no formulario tira a chave do arquivo, e a CLI volta a
    aplicar o proprio padrao -- que e onde os padroes moram. Gravar `null`
    criaria a segunda fonte de verdade que o arquivo antigo acumulou."""
    caminho = _escrever(tmp_path, {"etapa1": {"corte": 0.1, "limite": 10}})
    config.gravar({"etapa1": {"limite": None}}, caminho)
    depois = json.loads(open(caminho, encoding="utf-8").read())
    assert "limite" not in depois["etapa1"] and depois["etapa1"]["corte"] == 0.1


def test_gravar_e_lido_de_volta_por_ler(tmp_path):
    caminho = _escrever(tmp_path, {"etapa1": {"corte": 0.1}})
    config.gravar({"etapa1": {"corte": 0.2}}, caminho)
    assert config.ler(caminho)["etapa1"]["corte"] == 0.2


def test_gravar_devolve_o_caminho_escrito(tmp_path):
    caminho = _escrever(tmp_path, {"etapa1": {"corte": 0.1}})
    assert str(config.gravar({"etapa1": {"corte": 0.2}}, caminho)) == caminho


def test_gravar_em_arquivo_ausente_e_erro(tmp_path):
    """Gravar sobre configuracao que nao existe criaria um arquivo sem nenhuma
    das notas, e o painel teria apagado a documentacao dos numeros."""
    with pytest.raises(config.ErroDeConfig):
        config.gravar({"etapa1": {"corte": 0.2}}, tmp_path / "nao_existe.json")


def test_gravar_nao_escreve_ascii_escapado(tmp_path):
    """O arquivo e editado a mao tambem: acento escapado como \\u00e7 e
    ilegivel para quem abre o config."""
    caminho = _escrever(tmp_path, {"etapa1": {"_nota_corte": "corte de associação",
                                              "corte": 0.1}})
    config.gravar({"etapa1": {"corte": 0.2}}, caminho)
    bruto = open(caminho, encoding="utf-8").read()
    assert "\\u" not in bruto
