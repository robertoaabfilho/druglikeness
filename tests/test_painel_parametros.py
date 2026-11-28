# -*- coding: utf-8 -*-
"""
Correcao C2 -- o formulario de parametros.

O achado: as chaves `(str, list)` (`termos_alerta`, `tecidos`, `reprovar`,
`marcar`) tem a string separada por virgula como forma natural, e o formulario
fazia `", ".join(valor)` sobre uma STRING -- iterando caractere a caractere.
Um clique em "Gravar configuracao" transformava
`"reprovar": "alerta_bula_forte"` em `["a","l","e","r","t","a",...]`, e alerta
de bula e um dos dois unicos criterios que REPROVAM. A politica que o
pesquisador acredita estar valendo deixava de valer.

E o segundo achado, da mesma raiz: o formulario coletava e gravava as tres
secoes inteiras a cada clique, acrescentando `sem_rede: false` e removendo
`efo: []` sem ninguem ter tocado nos campos.

A regra que os dois consertos cumprem: **grava-se o que mudou, e nada mais.**
"""
import json

import pytest

from mods import config
from mods.painel import parametros


# --------------------------------------------------------------------------
# texto <-> valor
# --------------------------------------------------------------------------
def test_para_texto_de_string_nao_separa_caractere():
    """O bug C2, no menor teste que o pega."""
    assert parametros.para_texto("alerta_bula_forte") == "alerta_bula_forte"


def test_para_texto_de_string_com_virgula_e_preservada():
    assert parametros.para_texto("neuropathy,tremor") == "neuropathy,tremor"


def test_para_texto_de_lista_junta_com_virgula():
    assert parametros.para_texto(["a", "b"]) == "a, b"


def test_para_texto_de_none_e_vazio():
    assert parametros.para_texto(None) == ""


def test_de_texto_preserva_o_tipo_que_estava_no_arquivo():
    """Quem tinha string continua com string: reescrever como lista mudaria o
    arquivo de quem nunca pediu."""
    assert parametros.de_texto("a,b", anterior="x") == "a,b"
    assert parametros.de_texto("a, b", anterior=["x"]) == ["a", "b"]


def test_de_texto_sem_valor_anterior_vira_lista():
    assert parametros.de_texto("a, b", anterior=None) == ["a", "b"]


def test_de_texto_vazio_e_none():
    assert parametros.de_texto("  ", anterior=["x"]) is None


# --------------------------------------------------------------------------
# somente o que mudou
# --------------------------------------------------------------------------
def test_somente_mudancas_de_formulario_intocado_e_vazio():
    antes = {"corte": 0.1, "limite": 10, "sem_rede": False, "efo": []}
    assert parametros.somente_mudancas(antes, dict(antes)) == {}


def test_somente_mudancas_ignora_chave_ausente_no_arquivo():
    """Booleano que nao esta no arquivo chega ao formulario como `False`;
    gravar isso criaria a chave que ninguem pediu -- a segunda fonte de
    verdade que o config existe para nao ter."""
    assert parametros.somente_mudancas({}, {"sem_rede": False}) == {}


def test_somente_mudancas_grava_booleano_que_o_usuario_ligou():
    assert parametros.somente_mudancas({}, {"sem_rede": True}) == \
        {"sem_rede": True}


def test_somente_mudancas_nao_remove_lista_vazia_do_arquivo():
    """`efo: []` no arquivo vem ao formulario como texto vazio, que volta
    `None` -- e `None` REMOVE a chave, deixando `_nota_efo` documentando uma
    chave que saiu."""
    assert parametros.somente_mudancas({"efo": []}, {"efo": None}) == {}


def test_somente_mudancas_remove_quando_o_usuario_apaga_valor():
    assert parametros.somente_mudancas({"limite": 10}, {"limite": None}) == \
        {"limite": None}


def test_somente_mudancas_ve_a_mudanca_de_numero():
    assert parametros.somente_mudancas({"corte": 0.1}, {"corte": 0.2}) == \
        {"corte": 0.2}


def test_somente_mudancas_nao_confunde_0_com_ausente():
    """`corte: 0` e o modo denso, que inverte o resultado do enriquecimento:
    e um valor, nao uma ausencia."""
    assert parametros.somente_mudancas({"corte": 0.1}, {"corte": 0}) == \
        {"corte": 0}


def test_somente_mudancas_ignora_notas():
    assert parametros.somente_mudancas(
        {"_nota_corte": "x", "corte": 0.1}, {"corte": 0.1}) == {}


# --------------------------------------------------------------------------
# ida e volta sobre o config de verdade
# --------------------------------------------------------------------------
def test_ida_e_volta_no_config_versionado_nao_muda_o_arquivo(tmp_path):
    """O teste que faltava: abrir o painel, clicar "Gravar", nada muda.

    Roda sobre o `json/config.json` do repositorio, com as quatro chaves
    `(str, list)` preenchidas -- que e o estado em que o bug aparecia."""
    from mods import caminhos
    original = json.loads(caminhos.CONFIG.read_text(encoding="utf-8"))
    original["portao"]["termos_alerta"] = "neuropathy,tremor"
    original["portao"]["reprovar"] = "alerta_bula_forte"
    original["portao"]["tecidos"] = ["nerve", "brain"]
    caminho = tmp_path / "config.json"
    # com a nova-linha final, que e como `config.gravar` termina o arquivo
    caminho.write_text(json.dumps(original, ensure_ascii=False, indent=2) + "\n",
                       encoding="utf-8")
    antes = caminho.read_text(encoding="utf-8")

    cfg = config.ler(caminho)
    mudancas = {}
    for etapa in config.CHAVES:
        secao = cfg.get(etapa, {})
        # o que o formulario devolveria sem ninguem tocar em nada
        formulario = {
            chave: parametros.de_texto(
                parametros.para_texto(secao.get(chave)),
                anterior=secao.get(chave))
            if list in tipos or str in tipos else secao.get(chave)
            for chave, tipos in config.CHAVES[etapa].items()}
        mudancas[etapa] = parametros.somente_mudancas(secao, formulario)
    assert all(not m for m in mudancas.values()), mudancas
    config.gravar(mudancas, caminho)
    assert caminho.read_text(encoding="utf-8") == antes


def test_reprovar_como_string_sobrevive_a_gravacao(tmp_path):
    """O cenario exato do achado: `--reprovar a` quebraria a execucao com
    'criterio inexistente', e `--termos-alerta n` buscaria a letra n nas bulas."""
    caminho = tmp_path / "config.json"
    caminho.write_text(json.dumps(
        {"portao": {"reprovar": "alerta_bula_forte",
                    "termos_alerta": "neuropathy,tremor"}}), encoding="utf-8")
    cfg = config.ler(caminho)
    secao = cfg["portao"]
    formulario = {c: parametros.de_texto(parametros.para_texto(secao.get(c)),
                                         anterior=secao.get(c))
                  for c in ("reprovar", "termos_alerta")}
    config.gravar({"portao": parametros.somente_mudancas(secao, formulario)},
                  caminho)
    depois = json.loads(caminho.read_text(encoding="utf-8"))
    assert depois["portao"]["reprovar"] == "alerta_bula_forte"
    assert depois["portao"]["termos_alerta"] == "neuropathy,tremor"


def test_o_painel_nao_tem_mais_a_juncao_de_string(tmp_path):
    """Guarda: a juncao que iterava a string nao pode voltar ao painel."""
    from mods import caminhos
    fonte = (caminhos.RAIZ / "painel.py").read_text(encoding="utf-8")
    assert '", ".join(str(v) for v in atual)' not in fonte
