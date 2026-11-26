# -*- coding: utf-8 -*-
"""Testes do CLI da etapa 2 (`mods/aplicar_portao.py`).

Os criterios de verificacao do desenho da etapa 2 viram teste aqui. O ADMET-AI
e falsificado por padrao: a biblioteca real carrega dez modelos Chemprop, e o
que se testa e a costura, nao a rede neural. Um teste exercita a ausencia dela,
que e o modo degradado prometido.
"""

import json
from pathlib import Path

import pytest

from mods import aplicar_portao
from mods.portao import admet, bula, esquema

AMITRIPTILINA = "CN(C)CCC=C1c2ccccc2CCc2ccccc21"   # PAINS: styrene_A(13)
ASPIRINA = "CC(=O)Oc1ccccc1C(=O)O"                 # BRENK: phenol_ester
COLUNAS_FALSAS = ["hERG", "hERG_drugbank_approved_percentile",
                  "AMES", "AMES_drugbank_approved_percentile"]


@pytest.fixture
def ligantes(tmp_path):
    caminho = tmp_path / "ligantes.tsv"
    caminho.write_text(
        "nome\tsmiles\n"
        f"amitriptilina\t{AMITRIPTILINA}\n"
        f"aspirina\t{ASPIRINA}\n", encoding="utf-8")
    return str(caminho)


@pytest.fixture
def admet_falso(monkeypatch):
    """ADMET-AI presente, com duas colunas e percentis escolhidos: a
    amitriptilina acima do percentil 90 em hERG, a aspirina abaixo."""
    percentis = {AMITRIPTILINA: 95.0, ASPIRINA: 10.0}

    def prever(smiles, **kwargs):
        return {s: {"hERG": 0.8, "hERG_drugbank_approved_percentile": percentis.get(s, 50.0),
                    "AMES": 0.1, "AMES_drugbank_approved_percentile": 20.0}
                for s in smiles if s in percentis}

    monkeypatch.setattr(admet, "disponivel", lambda: True)
    monkeypatch.setattr(admet, "versao", lambda: "2.0.1-falso")
    monkeypatch.setattr(admet, "colunas", lambda num_workers=0: list(COLUNAS_FALSAS))
    monkeypatch.setattr(admet, "prever", prever)
    return percentis


def _rodar(*argv):
    return aplicar_portao.main(list(argv))


def _ler(pasta, nome):
    return list(esquema.ler_tsv(pasta / nome))


def _por_ligante(linhas):
    return {l["ligante"]: l for l in linhas}


# --------------------------------------------------------------------------
# Saidas, e o que NAO e escrito sem ser pedido
# --------------------------------------------------------------------------

def test_execucao_minima_publica_as_saidas_e_termina_em_zero(
        tmp_path, ligantes, admet_falso):
    saida = tmp_path / "res"
    assert _rodar("--ligantes", ligantes, "--saida", str(saida)) == 0
    assert (saida / "propriedades.tsv").is_file()
    assert (saida / "portao.tsv").is_file()
    assert (saida / "manifesto.json").is_file()
    # Bula e expressao nao foram pedidas: tabela vazia seria pior que ausente,
    # porque "sem alerta" e "nao perguntei" nao sao a mesma afirmacao.
    assert not (saida / "alertas_bula.tsv").exists()
    assert not (saida / "alvos_anotados.tsv").exists()


def test_propriedades_traz_descritor_alerta_nomeado_e_colunas_do_admet(
        tmp_path, ligantes, admet_falso):
    saida = tmp_path / "res"
    _rodar("--ligantes", ligantes, "--saida", str(saida))
    linhas = _por_ligante(_ler(saida, "propriedades.tsv"))
    assert "styrene_A(13)" in linhas["amitriptilina"]["alerta_pains"]
    assert "phenol_ester" in linhas["aspirina"]["alerta_brenk"]
    assert linhas["aspirina"]["qed"]
    # As colunas do ADMET-AI saem com os NOMES ORIGINAIS: e o que as torna
    # citaveis no artigo que cita a ferramenta.
    assert linhas["amitriptilina"]["hERG_drugbank_approved_percentile"] == "95.0"
    assert "RDKit" in linhas["aspirina"]["procedencia"]


# --------------------------------------------------------------------------
# Criterio 2 do desenho: marca, nao reprova
# --------------------------------------------------------------------------

def test_alerta_estrutural_marca_e_o_ligante_segue_aprovado(
        tmp_path, ligantes, admet_falso):
    """Reprovar por PAINS apagaria a amitriptilina, que TRATA dor neuropatica
    e e o positivo de referencia do pipeline."""
    saida = tmp_path / "res"
    _rodar("--ligantes", ligantes, "--saida", str(saida))
    linha = _por_ligante(_ler(saida, "portao.tsv"))["amitriptilina"]
    assert linha["alerta_estrutural_veredito"] == "marca"
    assert "styrene_A(13)" in linha["alerta_estrutural_valor"]
    assert linha["aprovado"] == "sim" and linha["motivos"] == ""


def test_toxicidade_acima_do_percentil_marca_sem_reprovar(
        tmp_path, ligantes, admet_falso):
    saida = tmp_path / "res"
    _rodar("--ligantes", ligantes, "--saida", str(saida))
    linhas = _por_ligante(_ler(saida, "portao.tsv"))
    assert linhas["amitriptilina"]["toxicidade_veredito"] == "marca"
    assert "hERG=95" in linhas["amitriptilina"]["toxicidade_valor"]
    assert linhas["aspirina"]["toxicidade_veredito"] == "passa"
    assert linhas["amitriptilina"]["aprovado"] == "sim"


def test_percentil_configuravel_muda_o_veredito(tmp_path, ligantes, admet_falso):
    saida = tmp_path / "res"
    _rodar("--ligantes", ligantes, "--saida", str(saida),
           "--percentil-toxicidade", "99")
    linhas = _por_ligante(_ler(saida, "portao.tsv"))
    assert linhas["amitriptilina"]["toxicidade_veredito"] == "passa"


def test_criterio_nao_pedido_e_indeterminado_nunca_passa(
        tmp_path, ligantes, admet_falso):
    """Bula nao consultada nao e "sem alerta"."""
    saida = tmp_path / "res"
    _rodar("--ligantes", ligantes, "--saida", str(saida))
    linha = _por_ligante(_ler(saida, "portao.tsv"))["aspirina"]
    assert linha["alerta_bula_veredito"] == "indeterminado"
    assert "termos-alerta" in linha["alerta_bula_valor"]


# --------------------------------------------------------------------------
# Criterio 3 do desenho: a politica, e o resumo que mede o botao
# --------------------------------------------------------------------------

def test_reprovar_promove_o_criterio_e_a_amitriptilina_cai(
        tmp_path, ligantes, admet_falso):
    saida = tmp_path / "res"
    assert _rodar("--ligantes", ligantes, "--saida", str(saida),
                  "--reprovar", "alerta_estrutural") == 0
    linha = _por_ligante(_ler(saida, "portao.tsv"))["amitriptilina"]
    assert linha["aprovado"] == "nao"
    assert "alerta_estrutural" in linha["motivos"]
    assert "alerta_estrutural=reprovar" in linha["politica"]


def test_reprovado_continua_na_saida_com_o_motivo(tmp_path, ligantes, admet_falso):
    """Apagar em silencio e anti-objetivo declarado do desenho."""
    saida = tmp_path / "res"
    _rodar("--ligantes", ligantes, "--saida", str(saida),
           "--reprovar", "alerta_estrutural")
    assert "amitriptilina" in _por_ligante(_ler(saida, "portao.tsv"))
    assert "amitriptilina" in _por_ligante(_ler(saida, "propriedades.tsv"))


def test_resumo_diz_quantos_cada_criterio_reprovaria_se_promovido(
        tmp_path, ligantes, admet_falso, capsys):
    saida = tmp_path / "res"
    _rodar("--ligantes", ligantes, "--saida", str(saida))
    texto = capsys.readouterr().err
    assert "reprovaria" in texto
    assert "alerta_estrutural" in texto


def test_criterio_desconhecido_erra_cedo_e_nao_publica_nada(tmp_path, ligantes):
    saida = tmp_path / "res"
    assert _rodar("--ligantes", ligantes, "--saida", str(saida),
                  "--reprovar", "alerta_estruturall") == 2
    assert not (saida / "portao.tsv").exists()


# --------------------------------------------------------------------------
# Criterio 4 do desenho: degradacao sem admet_ai
# --------------------------------------------------------------------------

def test_sem_admet_ai_termina_em_zero_com_toxicidade_indeterminada(
        tmp_path, ligantes, monkeypatch, capsys):
    """Sem a biblioteca a etapa RODA e diz que rodou sem ADMET. Ausencia nunca
    vira aprovacao nem reprovacao silenciosa."""
    monkeypatch.setattr(admet, "disponivel", lambda: False)
    saida = tmp_path / "res"
    assert _rodar("--ligantes", ligantes, "--saida", str(saida)) == 0
    linha = _por_ligante(_ler(saida, "portao.tsv"))["amitriptilina"]
    assert linha["toxicidade_veredito"] == "indeterminado"
    assert linha["aprovado"] == "sim"
    # As 104 colunas NAO sao fixadas no codigo: sem a biblioteca instalada nao
    # ha como nomea-las, e inventa-las faria o repositorio mentir sobre a
    # versao. O cabecalho sai sem elas, e o manifesto registra a ausencia.
    propriedades = _ler(saida, "propriedades.tsv")[0]
    assert "hERG" not in propriedades
    assert "alerta_pains" in propriedades
    manifesto = json.loads((saida / "manifesto.json").read_text(encoding="utf-8"))
    assert manifesto["versoes"]["admet_ai"] == "ausente"
    assert "sem ADMET" in capsys.readouterr().err


# --------------------------------------------------------------------------
# Manifesto, procedencia e codigos de saida
# --------------------------------------------------------------------------

def test_manifesto_registra_politica_versoes_e_contagens(
        tmp_path, ligantes, admet_falso):
    saida = tmp_path / "res"
    _rodar("--ligantes", ligantes, "--saida", str(saida), "--marcar", "alerta_bula")
    manifesto = json.loads((saida / "manifesto.json").read_text(encoding="utf-8"))
    assert manifesto["parametros"]["politica"]["alerta_bula"] == "marcar"
    assert manifesto["parametros"]["percentil_toxicidade"] == 90
    assert manifesto["versoes"]["admet_ai"] == "2.0.1-falso"
    assert manifesto["versoes"]["ligantes_md5"]
    assert manifesto["versoes"]["pipeline_commit"]
    assert manifesto["contagens"]["ligantes"] == 2
    assert manifesto["contagens"]["aprovados"] == 2
    # O numero que o resumo anuncia: quanto custaria promover o criterio.
    por_criterio = manifesto["contagens"]["por_criterio"]
    assert por_criterio["alerta_estrutural"]["adversos"] == 2


def test_linha_de_entrada_descartada_torna_a_execucao_parcial(
        tmp_path, admet_falso, capsys):
    # Linha sem SMILES nenhum: sem estrutura nao ha o que julgar nem o que
    # relatar, entao ela e PERDIDA -- diferente do SMILES ilegivel, que fica na
    # saida e e reprovado (teste proprio, mais abaixo).
    caminho = tmp_path / "com_erro.tsv"
    caminho.write_text(f"nome\tsmiles\naspirina\t{ASPIRINA}\nsem_nada\t\n",
                       encoding="utf-8")
    saida = tmp_path / "res"
    assert _rodar("--ligantes", str(caminho), "--saida", str(saida)) == 3
    manifesto = json.loads((saida / "manifesto.json").read_text(encoding="utf-8"))
    assert manifesto["contagens"]["execucao_parcial"] is True
    assert manifesto["contagens"]["linhas_descartadas"] == 1
    # A saida existe, mas o usuario precisa saber que esta incompleta.
    assert "PARCIAL" in capsys.readouterr().err


def test_entrada_sem_ligante_utilizavel_para_sem_publicar(tmp_path, admet_falso):
    caminho = tmp_path / "vazio.tsv"
    caminho.write_text("nome\tsmiles\n", encoding="utf-8")
    saida = tmp_path / "res"
    assert _rodar("--ligantes", str(caminho), "--saida", str(saida)) == 1
    assert not (saida / "portao.tsv").exists()


def test_limite_corta_a_entrada(tmp_path, ligantes, admet_falso):
    saida = tmp_path / "res"
    _rodar("--ligantes", ligantes, "--saida", str(saida), "--limite", "1")
    assert len(_ler(saida, "portao.tsv")) == 1


def test_saida_apontando_para_arquivo_e_erro(tmp_path, ligantes, admet_falso):
    arquivo = tmp_path / "ocupado"
    arquivo.write_text("nao sou pasta\n", encoding="utf-8")
    assert _rodar("--ligantes", ligantes, "--saida", str(arquivo)) == 2
    assert arquivo.read_text(encoding="utf-8") == "nao sou pasta\n"


def test_ordem_das_linhas_e_deterministica(tmp_path, ligantes, admet_falso):
    """Necessaria para a repeticao: dois arquivos iguais byte a byte."""
    um, dois = tmp_path / "a", tmp_path / "b"
    _rodar("--ligantes", ligantes, "--saida", str(um))
    _rodar("--ligantes", ligantes, "--saida", str(dois))
    for nome in ("propriedades.tsv", "portao.tsv"):
        assert (um / nome).read_bytes() == (dois / nome).read_bytes()


# --------------------------------------------------------------------------
# Colisao de nome de coluna (achada rodando a etapa em input.exemplo.tsv)
# --------------------------------------------------------------------------

def test_coluna_do_admet_que_colide_com_o_esquema_e_prefixada(
        tmp_path, ligantes, monkeypatch):
    """O ADMET-AI 2.0.1 publica `tpsa`, nome que o esquema do repositorio ja
    usa para a TPSA do RDKit. Sem tratamento o cabecalho saia com `tpsa` duas
    vezes e o valor do RDKit era substituido pelo da outra ferramenta, em
    silencio -- duas colunas respondendo a mesma pergunta sem procedencia."""
    monkeypatch.setattr(admet, "disponivel", lambda: True)
    monkeypatch.setattr(admet, "versao", lambda: "2.0.1-falso")
    monkeypatch.setattr(admet, "colunas", lambda num_workers=0: ["tpsa", "hERG"])
    monkeypatch.setattr(admet, "prever", lambda smiles, **k: {
        s: {"tpsa": 999.0, "hERG": 0.5} for s in smiles})
    saida = tmp_path / "res"
    assert _rodar("--ligantes", ligantes, "--saida", str(saida)) == 0

    cabecalho = (saida / "propriedades.tsv").read_text(
        encoding="utf-8").splitlines()[0].lstrip("#").split("\t")
    assert len(cabecalho) == len(set(cabecalho)), "coluna duplicada no cabecalho"
    assert "admet_tpsa" in cabecalho and "hERG" in cabecalho

    linha = _por_ligante(_ler(saida, "propriedades.tsv"))["aspirina"]
    assert linha["admet_tpsa"] == "999.0"          # o valor do ADMET-AI
    assert linha["tpsa"] not in ("", "999.0")      # e o do RDKit, preservado

    manifesto = json.loads((saida / "manifesto.json").read_text(encoding="utf-8"))
    assert manifesto["contagens"]["colunas_de_admet_renomeadas"] == {
        "tpsa": "admet_tpsa"}


def test_desvio_impossivel_erra_em_vez_de_adivinhar(tmp_path, ligantes, monkeypatch):
    """Se `admet_tpsa` tambem existisse, inventar um terceiro nome seria pior
    que parar."""
    monkeypatch.setattr(admet, "disponivel", lambda: True)
    monkeypatch.setattr(admet, "versao", lambda: "x")
    monkeypatch.setattr(admet, "colunas",
                        lambda num_workers=0: ["tpsa", "admet_tpsa"])
    monkeypatch.setattr(admet, "prever", lambda smiles, **k: {})
    assert _rodar("--ligantes", ligantes, "--saida", str(tmp_path / "res")) == 2
    assert not (tmp_path / "res" / "propriedades.tsv").exists()


# --------------------------------------------------------------------------
# SMILES ilegivel: reprovado com motivo, nunca apagado
# --------------------------------------------------------------------------

def test_smiles_ilegivel_fica_na_saida_reprovado_com_motivo(
        tmp_path, admet_falso, capsys):
    """Anti-objetivo do desenho: apagar em silencio. O ligante nao desaparece
    da tabela de vereditos so porque nao pode ser avaliado."""
    caminho = tmp_path / "ligantes.tsv"
    caminho.write_text(f"nome\tsmiles\naspirina\t{ASPIRINA}\nruim\tC1CC\n",
                       encoding="utf-8")
    saida = tmp_path / "res"
    # Codigo 0, nao 3: o portao DECIDIU sobre o ligante. Parcial e quando algo
    # se perde, e aqui nada se perdeu.
    assert _rodar("--ligantes", str(caminho), "--saida", str(saida)) == 0

    linha = _por_ligante(_ler(saida, "portao.tsv"))["ruim"]
    assert linha["aprovado"] == "nao"
    assert linha["smiles_veredito"] == "reprova"
    assert "RDKit nao le" in linha["motivos"]
    # Os outros criterios nao opinam sobre o que nao tem molecula.
    assert linha["alerta_estrutural_veredito"] == "indeterminado"
    assert linha["toxicidade_veredito"] == "indeterminado"

    propriedades = _por_ligante(_ler(saida, "propriedades.tsv"))["ruim"]
    assert propriedades["smiles"] == "C1CC" and propriedades["qed"] == ""

    manifesto = json.loads((saida / "manifesto.json").read_text(encoding="utf-8"))
    assert manifesto["contagens"]["smiles_ilegiveis_mantidos"] == 1
    assert manifesto["contagens"]["execucao_parcial"] is False
    assert "nao apagados" in capsys.readouterr().err


# --------------------------------------------------------------------------
# Criterio 1 do desenho: a bula no conjunto-controle de 12 farmacos
# --------------------------------------------------------------------------

TERMOS = "neuropathy,neuropathic,neuritis,paresthesia"
CONTROLE = "exemplos/controle.dor_neuropatica.tsv"
TRATAM = ["pregabalin", "duloxetine", "gabapentin", "amitriptyline",
          "carbamazepine", "capsaicin"]
CAUSAM = ["ciprofloxacin", "linezolid", "chloramphenicol", "ethambutol",
          "ethionamide", "metronidazole"]


@pytest.fixture
def cache_de_bulas(monkeypatch, respostas):
    """Cache falso servido pelas respostas gravadas do openFDA. A fixture
    `respostas` FALHA o teste para qualquer URL nao gravada, entao nenhum
    caminho novo passa escondido -- e o bloqueio de socket garante que nada
    cai para a rede."""
    class CacheFalso:
        obter_json = staticmethod(respostas)
        diretorio = Path("openfda-2026-11-25")

    monkeypatch.setattr(
        aplicar_portao, "_abrir_cache",
        lambda args, destino: (CacheFalso(), {"openfda": "2026-11-25"},
                               "/raiz/cache", "openfda-2026-11-25"))


@pytest.fixture
def sem_admet(monkeypatch):
    """O controle da bula nao precisa de ADMET, e dez modelos Chemprop por
    teste custariam mais que todo o resto da suite."""
    monkeypatch.setattr(admet, "disponivel", lambda: False)


def _rodar_controle(saida, *extra):
    return _rodar("--ligantes", CONTROLE, "--saida", str(saida),
                  "--termos-alerta", TERMOS, *extra)


def test_bula_no_controle_tem_zero_falso_positivo_e_pega_quatro_de_seis(
        tmp_path, cache_de_bulas, sem_admet):
    """O numero que justifica `alerta_bula` ser o UNICO criterio que reprova
    por padrao. A regra antiga dava 4/6 de FALSO positivo e perdia o
    ciprofloxacino; esta da 0/6 de falso positivo e o pega."""
    saida = tmp_path / "res"
    assert _rodar_controle(saida) == 0
    linhas = _por_ligante(_ler(saida, "portao.tsv"))

    reprovados_que_tratam = [n for n in TRATAM
                             if linhas[n]["alerta_bula_veredito"] == "reprova"]
    assert reprovados_que_tratam == [], (
        "falso positivo: estes TRATAM dor neuropatica")

    reprovados_que_causam = [n for n in CAUSAM
                             if linhas[n]["alerta_bula_veredito"] == "reprova"]
    assert reprovados_que_causam == ["ciprofloxacin", "linezolid",
                                     "ethambutol", "metronidazole"]

    # Os dois nao pegos sao honestos, nao "sem alerta": cloranfenicol tem UMA
    # bula (so secao fraca) e etionamida nenhuma.
    assert linhas["chloramphenicol"]["alerta_bula_veredito"] == "passa"
    assert "fraca" in linhas["chloramphenicol"]["alerta_bula_valor"]
    assert linhas["ethionamide"]["alerta_bula_veredito"] == "indeterminado"
    assert "sem bula" in linhas["ethionamide"]["alerta_bula_valor"]


def test_ciprofloxacino_e_pego_pela_tarja_preta(tmp_path, cache_de_bulas, sem_admet):
    """O caso que originou o projeto: `limit=1` devolvia uma bula arbitraria
    entre 166 e a que traz a tarja preta nao e a primeira."""
    saida = tmp_path / "res"
    _rodar_controle(saida)
    linha = _por_ligante(_ler(saida, "portao.tsv"))["ciprofloxacin"]
    assert "boxed_warning" in linha["alerta_bula_valor"]
    assert linha["aprovado"] == "nao"
    assert "alerta_bula" in linha["motivos"]


def test_capsaicina_passa_por_indicacao_e_o_peso_fica_registrado(
        tmp_path, cache_de_bulas, sem_admet):
    """A capsaicina TRATA dor neuropatica e nao pode ser reprovada.

    Ela tem 211 bulas no openFDA e UMA so a indica para a condicao (Qutenza).
    Dois mecanismos independentes a salvam hoje, e isso e deliberado:

    * a janela passou a ser ordenada por data, e a Qutenza entra nas 50 mais
      recentes (effective_time 2026-07-07);
    * a consulta dirigida de indicacao roda no acervo inteiro e nao depende de
      a Qutenza ser recente.

    A redundancia e o ponto: a ordenacao salvou este caso por SORTE DE DATA (1
    bula de 211), e uma bula antiga com genericos novos por cima escaparia.
    `test_dirigido_decide_sobre_a_janela`, em test_portao.py, prende a regra
    sem depender destes dados.
    """
    saida = tmp_path / "res"
    _rodar_controle(saida)
    linha = _por_ligante(_ler(saida, "portao.tsv"))["capsaicin"]
    assert linha["alerta_bula_veredito"] == "passa"
    assert "indicado" in linha["alerta_bula_valor"]
    # O peso da evidencia, que uma contagem dentro da janela nao saberia dizer.
    assert "1 de 211 bulas" in linha["alerta_bula_valor"]
    assert linha["aprovado"] == "sim"


def test_peso_da_evidencia_distingue_118_de_166_de_1_de_211(
        tmp_path, cache_de_bulas, sem_admet):
    """O ciprofloxacino e reprovado por 118 das 166 bulas; a capsaicina e
    salva por 1 de 211. Os dois numeros vem do acervo, nao da janela, e e a
    diferenca entre um alerta robusto e um fio de cabelo."""
    saida = tmp_path / "res"
    _rodar_controle(saida)
    linhas = _por_ligante(_ler(saida, "portao.tsv"))
    assert "118 de 166 bulas" in linhas["ciprofloxacin"]["alerta_bula_valor"]
    assert "35 de 35 bulas" in linhas["linezolid"]["alerta_bula_valor"]


def test_alertas_bula_guarda_o_rastro_conferivel(tmp_path, cache_de_bulas, sem_admet):
    """Sem o trecho e a secao, "tem alerta" nao se confere. O rastro e o
    produto desta tabela."""
    saida = tmp_path / "res"
    _rodar_controle(saida)
    linhas = [l for l in _ler(saida, "alertas_bula.tsv")
              if l["ligante"] == "ciprofloxacin"]
    assert linhas, "nenhum achado gravado para o ciprofloxacino"
    forte = [l for l in linhas if l["secao"] == "boxed_warning"]
    assert forte and forte[0]["forca"] == "forte"
    assert "neuropathy" in forte[0]["trecho"].lower()
    assert int(forte[0]["n_bulas_consultadas"]) == 50
    assert forte[0]["campo_do_nome"] == "generic_name"
    assert forte[0]["procedencia"] == "openFDA 2026-11-25"
    # A secao fraca tambem fica registrada, mesmo nao contando como alerta.
    assert any(l["forca"] == "fraca" for l in _ler(saida, "alertas_bula.tsv"))


def test_teto_de_bula_sai_indeterminado_com_o_motivo(
        tmp_path, cache_de_bulas, sem_admet):
    """108.935 consultas a um servico publico com limite de taxa nao sao
    aceitaveis. O que passa do teto nunca vira "sem alerta"."""
    saida = tmp_path / "res"
    assert _rodar_controle(saida, "--limite-bula", "1") == 0
    linhas = _ler(saida, "portao.tsv")
    assert linhas[0]["alerta_bula_veredito"] != "indeterminado"
    acima = [l for l in linhas[1:] if l["alerta_bula_veredito"] == "indeterminado"]
    assert len(acima) == len(linhas) - 1
    assert "--limite-bula 1" in acima[0]["alerta_bula_valor"]
    manifesto = json.loads((saida / "manifesto.json").read_text(encoding="utf-8"))
    assert manifesto["contagens"]["bulas_consultadas"] == 1


def test_termos_alerta_com_sem_rede_e_erro(tmp_path, ligantes):
    assert _rodar("--ligantes", ligantes, "--saida", str(tmp_path / "res"),
                  "--termos-alerta", "neuropathy", "--sem-rede") == 2
    assert not (tmp_path / "res" / "portao.tsv").exists()


def test_falha_de_consulta_num_item_torna_a_execucao_parcial(
        tmp_path, monkeypatch, respostas, sem_admet, capsys):
    """Rede caida em UM ligante nao derruba a execucao nem vira "sem alerta":
    aquele item fica `indeterminado` e a execucao inteira fica PARCIAL."""
    def as_vezes(url, corpo=None):
        if "linezolid" in url:
            raise RuntimeError("502 do openFDA")
        return respostas(url, corpo)

    class CacheFalso:
        obter_json = staticmethod(as_vezes)
        diretorio = Path("p")

    monkeypatch.setattr(aplicar_portao, "_abrir_cache",
                        lambda args, destino: (CacheFalso(), {"openfda": "x"},
                                               "/r", "p"))
    saida = tmp_path / "res"
    assert _rodar_controle(saida) == 3
    linha = _por_ligante(_ler(saida, "portao.tsv"))["linezolid"]
    assert linha["alerta_bula_veredito"] == "indeterminado"
    assert "falhou" in linha["alerta_bula_valor"]
    assert linha["aprovado"] == "sim"       # indeterminado nunca reprova
    manifesto = json.loads((saida / "manifesto.json").read_text(encoding="utf-8"))
    assert manifesto["contagens"]["falhas_de_bula"]
    assert "PARCIAL" in capsys.readouterr().err


def test_manifesto_registra_a_versao_do_openfda_e_os_termos(
        tmp_path, cache_de_bulas, sem_admet):
    saida = tmp_path / "res"
    _rodar_controle(saida)
    manifesto = json.loads((saida / "manifesto.json").read_text(encoding="utf-8"))
    assert manifesto["versoes"]["openfda"] == "2026-11-25"
    assert manifesto["parametros"]["termos_alerta"] == TERMOS.split(",")
    assert manifesto["parametros"]["limite_bula"] == 500
    assert manifesto["parametros"]["cache_pasta"] == "openfda-2026-11-25"


def test_janela_e_pedida_ordenada_por_data(tmp_path, monkeypatch, respostas,
                                           sem_admet):
    """"As 50 mais acessadas" nao existe: a bula do openFDA tem 34 campos e
    nenhum de acesso, popularidade ou contagem de uso. A data e o unico
    critetio de relevancia disponivel, e significa rotulagem vigente -- em vez
    de 50 bulas sorteadas entre 211."""
    pedidas = []

    def registrando(url, corpo=None):
        pedidas.append(url)
        return respostas(url, corpo)

    class CacheFalso:
        obter_json = staticmethod(registrando)
        diretorio = Path("p")

    monkeypatch.setattr(aplicar_portao, "_abrir_cache",
                        lambda args, destino: (CacheFalso(), {"openfda": "x"},
                                               "/r", "p"))
    _rodar_controle(tmp_path / "res", "--limite-bula", "1")
    janelas = [u for u in pedidas if "limit=50" in u]
    assert janelas, "nenhuma janela pedida"
    assert all("sort=effective_time:desc" in u for u in janelas)


def test_tres_requisicoes_por_ligante_janela_forte_e_indicacao(
        tmp_path, monkeypatch, respostas, sem_admet):
    """O custo por ligante e conhecido e fica prendido: a secao fraca NAO tem
    consulta propria, porque nao decide nada -- a janela basta para
    registra-la."""
    pedidas = []

    def registrando(url, corpo=None):
        pedidas.append(url)
        return respostas(url, corpo)

    class CacheFalso:
        obter_json = staticmethod(registrando)
        diretorio = Path("p")

    monkeypatch.setattr(aplicar_portao, "_abrir_cache",
                        lambda args, destino: (CacheFalso(), {"openfda": "x"},
                                               "/r", "p"))
    _rodar_controle(tmp_path / "res", "--limite", "1", "--limite-bula", "1")
    assert len(pedidas) == 3, pedidas
    assert sum("limit=50" in u for u in pedidas) == 1
    assert sum("indications_and_usage" in u for u in pedidas) == 1
    assert sum("boxed_warning" in u for u in pedidas) == 1
