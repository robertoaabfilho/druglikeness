# -*- coding: utf-8 -*-
"""A fase 3 no REGIME DE PRODUCAO da etapa 1.

Por que este arquivo existe
---------------------------
A suite inteira passava verde sobre uma afirmacao falsa: que as faixas `1b` e
`2b` sao inalcancaveis em producao, porque a etapa 1 escreveria uma linha
prevista para todo par composto x alvo e a faixa `2a` nunca ficaria vazia.

Isso vale SO com `--corte 0,0`. O corte e um filtro pos-predicao
(`mods/evidencias/predito.py:152`, `if p < corte: continue`) e o padrao e 0,1.
No corte padrao, dos 60.384 pares do run real sobram 297 linhas, e
`classificar` poe 51 compostos em `2b` contra 19 em `2a` -- `2b` e a faixa
MAIS comum.

Nenhuma fixture do repositorio tinha a densidade nem o corte da etapa 1 real,
entao nenhum teste podia contradizer a afirmacao. Esta tem as duas coisas, e o
teste do piso de `valor` e o que impede que ela seja regenerada no regime
errado sem a suite notar.

Desenho: `docs/superpowers/specs/2026-11-26-fase3-topologia-design.md`.
"""

import collections
import inspect
import os

from mods.metricas import calcular_auc
from mods.evidencias import predito
from mods.evidencias.esquema import COLUNAS_EVIDENCIA
from mods.ranque.faixa import classificar

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CAMINHO = os.path.join(RAIZ, "tests/dados/previsoes_regime_producao.tsv")
ROTULOS = os.path.join(RAIZ, "dados/validacao_neuro_prospector.tsv")
CONTROLE = os.path.join(RAIZ, "tests/dados/controle_previsoes.tsv")

AFIRMATIVOS = {"usuario", "open_targets", "reactome"}

# O padrao de `mods/evidencias/predito.py`. Nao e um parametro deste teste: e o
# regime que a fixture precisa representar.
CORTE_DE_PRODUCAO = 0.1


def _cabecalho(caminho) -> list[str]:
    with open(caminho, encoding="utf-8") as arquivo:
        return arquivo.readline().lstrip("#").rstrip("\n").split("\t")


def linhas_da_fixture() -> list[dict]:
    colunas = _cabecalho(CAMINHO)
    with open(CAMINHO, encoding="utf-8") as arquivo:
        arquivo.readline()
        return [dict(zip(colunas, l.rstrip("\n").split("\t")))
                for l in arquivo if l.strip()]


# --------------------------------------------------------------------------
# A forma da fixture
# --------------------------------------------------------------------------

def test_a_fixture_tem_o_contrato_da_etapa_1():
    """Se a etapa 1 ganhar ou renomear uma coluna, a fixture passa a
    representar um contrato que nao existe mais."""
    assert list(COLUNAS_EVIDENCIA) == _cabecalho(CAMINHO)


def test_a_fixture_tem_a_forma_do_run_real():
    linhas = linhas_da_fixture()
    assert len(linhas) == 297
    assert len({l["ligante_inchikey_parental"] for l in linhas}) == 70


def test_nenhuma_linha_abaixo_do_corte_de_producao():
    """A guarda central deste arquivo. Uma fixture regenerada com `--corte 0`
    traria a matriz densa (60.384 linhas) e todo composto cairia em `2a`."""
    assert min(float(l["valor"]) for l in linhas_da_fixture()) >= CORTE_DE_PRODUCAO


def test_a_fixture_e_so_de_evidencia_prevista():
    assert {l["tipo_evidencia"] for l in linhas_da_fixture()} == {"predito"}


def test_as_tres_procedencias_estao_representadas():
    assert {l["origem_alvo"] for l in linhas_da_fixture()} == {
        "usuario", "open_targets", "reactome"}


# --------------------------------------------------------------------------
# O que a fase 3 faz nesse regime
# --------------------------------------------------------------------------

def _por_composto() -> dict[str, list[dict]]:
    agrupadas: dict[str, list[dict]] = collections.defaultdict(list)
    for linha in linhas_da_fixture():
        agrupadas[linha["ligante_inchikey_parental"]].append(linha)
    return dict(agrupadas)


def _rotulos() -> dict[str, dict]:
    """O conjunto de validacao VERSIONADO, casado por forma parental."""
    colunas = _cabecalho(ROTULOS)
    with open(ROTULOS, encoding="utf-8") as arquivo:
        arquivo.readline()
        linhas = [dict(zip(colunas, l.rstrip("\n").split("\t")))
                  for l in arquivo if l.strip()]
    return {l["inchikey_parental"]: l for l in linhas}


def _faixas_por_composto() -> dict[str, int]:
    """Quantos compostos a fase 3 poe em cada faixa, SEM canal negativo.

    Sem canal de proposito: nestes 222 compostos o canal mede AUC 0,5156 --
    ruido --, e `None` nao e zero (zero afirmaria 'nenhuma atividade
    indesejada', e ausencia de dado nao e afirmacao)."""
    return dict(collections.Counter(
        classificar(linhas, p_negativo=None)["faixa"]
        for linhas in _por_composto().values()))


def _sem_rotulo() -> list[str]:
    rotulos = _rotulos()
    return sorted(k for k in _por_composto() if k not in rotulos)


def _nomes_do_controle() -> set[str]:
    """Os 12 farmacos que DERIVARAM as decisoes da fase 3, lidos da fixture
    que os guarda -- nao da coluna `no_controle_de_desenho`, que e mantida a
    mao e pode ficar atras do controle sem ninguem notar."""
    colunas = _cabecalho(CONTROLE)
    with open(CONTROLE, encoding="utf-8") as arquivo:
        arquivo.readline()
        return {dict(zip(colunas, l.rstrip("\n").split("\t")))["nome"].casefold()
                for l in arquivo if l.strip()}


def _do_controle_de_desenho(nomes: set[str] | None = None) -> list[str]:
    """Compostos da fixture que estao no controle de desenho. `nomes` existe
    para que a guarda possa ser provada contra um caso positivo."""
    nomes = _nomes_do_controle() if nomes is None else nomes
    rotulos = _rotulos()
    return sorted(rotulos[k]["nome"] for k in _por_composto()
                  if rotulos[k]["nome"].casefold() in nomes)


def _contagem_de_grupos() -> dict[str, int]:
    rotulos = _rotulos()
    return dict(collections.Counter(
        rotulos[k]["grupo_esperado"] for k in _por_composto()))


def _melhor_p(linhas, origens: set[str]) -> float | None:
    """`max` do valor previsto entre as linhas cuja procedencia esta em
    `origens`, ou None se nenhuma esta. `max` e nao soma: somar
    probabilidades de alvos diferentes soma grandezas que nao se somam, e com
    274 alvos a soma vira quase literalmente a contagem de promiscuidade --
    medida e rejeitada (AUC 0,583)."""
    valores = [float(l["valor"]) for l in linhas
               if l["origem_alvo"] in origens]
    return max(valores) if valores else None


def _auc_de(escores: dict[str, float]) -> float:
    """Pares concordantes, empate valendo meio. `trata` e o positivo.

    A estatistica de Mann-Whitney mora em `mods/metricas.py`, num lugar so.
    Aqui se faz apenas a separacao entre quem trata e quem causa."""
    rotulos = _rotulos()
    trata = [v for k, v in escores.items()
             if rotulos[k]["grupo_esperado"] == "trata"]
    causa = [v for k, v in escores.items()
             if rotulos[k]["grupo_esperado"] == "causa"]
    return calcular_auc(trata, causa)


def _auc(origens: set[str]) -> float:
    escores = {}
    for inchikey, linhas in _por_composto().items():
        melhor = _melhor_p(linhas, origens)
        if melhor is not None:
            escores[inchikey] = melhor
    return _auc_de(escores)


def _auc_das_metades() -> list[tuple[float, int]]:
    """[(AUC, n) da metade A, (AUC, n) da metade B], todos os alvos
    afirmativos. A divisao e a que o desenho fixou: os 222 compostos
    ROTULADOS (`trata`/`causa`, nunca `ambiguo`) ordenados por
    `inchikey_parental`, alternando. A populacao e o conjunto rotulado, nao as
    244 linhas do arquivo -- dividir o arquivo inteiro daria outras metades."""
    rotulos = _rotulos()
    rotulados = sorted(k for k, v in rotulos.items()
                       if v["grupo_esperado"] in ("trata", "causa"))
    metade_a = {k for i, k in enumerate(rotulados) if i % 2 == 0}
    metades = (metade_a, set(rotulados) - metade_a)
    saida = []
    for metade in metades:
        escores = {}
        for inchikey, linhas in _por_composto().items():
            if inchikey not in metade:
                continue
            escores[inchikey] = _melhor_p(linhas, AFIRMATIVOS)
        saida.append((_auc_de(escores), len(escores)))
    return saida


def _auc_na_base_comum() -> float:
    """Todos os alvos afirmativos, mas so nos compostos que TAMBEM tem linha
    em alvo do usuario -- a unica comparacao com denominador igual."""
    escores = {}
    for inchikey, linhas in _por_composto().items():
        if _melhor_p(linhas, {"usuario"}) is None:
            continue
        escores[inchikey] = _melhor_p(linhas, AFIRMATIVOS)
    return _auc_de(escores)

def test_a_faixa_2b_e_alcancavel_no_corte_padrao():
    """O achado que retirou o desenho de topologia: a premissa de que `2b` e
    inalcancavel vale SO no corte 0, onde a etapa 1 escreve a matriz densa."""
    assert _faixas_por_composto() == {"2b": 51, "2a": 19}


def test_todo_composto_da_fixture_tem_rotulo_versionado():
    """Sem isto a AUC mudaria de denominador em silencio."""
    assert _sem_rotulo() == []


def test_o_corte_de_producao_e_o_padrao_da_etapa_1():
    """O literal 0,1 deste arquivo tem de SER o padrao real da etapa 1.

    Sem este teste, mudar o padrao para 0,05 ou 0,2 deixaria os outros doze
    verdes enquanto a fixture silenciosamente parava de representar producao
    -- o unico buraco da guarda do piso."""
    padrao = inspect.signature(
        predito.linhas_preditas).parameters["corte"].default
    assert padrao == CORTE_DE_PRODUCAO


def test_nenhum_composto_da_fixture_esta_no_controle_de_desenho():
    """Um composto dos 12 que derivaram a fase 3 tornaria a AUC de treino.

    Derivado de `controle_previsoes.tsv`, NAO da coluna
    `no_controle_de_desenho`: aquela coluna e mantida a mao e marca um unico
    composto (capsaicina), que ja sai do conjunto por ser `ambiguo` -- uma
    guarda que lia so a coluna nao podia falhar."""
    assert _do_controle_de_desenho() == []


def test_a_guarda_de_contaminacao_pega_um_composto_do_controle():
    """Prova que a guarda acima morde: EMETINE esta na fixture, entao se ela
    estivesse no controle de desenho seria apontada."""
    assert _do_controle_de_desenho({"emetine"}) == ["EMETINE"]


def test_a_coluna_do_conjunto_versionado_concorda_com_o_controle():
    """A coluna mantida a mao, conferida contra os 12 nomes do controle nas
    244 linhas do arquivo versionado."""
    nomes = _nomes_do_controle()
    rotulos = _rotulos().values()
    marcados = {r["nome"].casefold() for r in rotulos
                if r["no_controle_de_desenho"] == "sim"}
    derivados = {r["nome"].casefold() for r in rotulos
                 if r["nome"].casefold() in nomes}
    assert marcados == derivados == {"capsaicin"}


def test_os_rotulos_da_fatia():
    assert _contagem_de_grupos() == {"trata": 55, "causa": 15}


def test_auc_no_regime_de_producao_so_alvos_do_usuario():
    assert round(_auc({"usuario"}), 4) == 0.6146          # n = 19


def test_auc_no_regime_de_producao_todos_os_alvos_afirmativos():
    """Deixar o enriquecimento competir PIORA no corte padrao: 0,6146 ->
    0,5145. No corte 0 o mesmo arquivo dava 0,4842 -> 0,6007, invertido."""
    assert round(_auc(AFIRMATIVOS), 4) == 0.5145


def test_na_base_comum_de_compostos_o_enriquecimento_fica_abaixo_do_acaso():
    assert round(_auc_na_base_comum(), 4) == 0.3438       # n = 19


def test_as_duas_metades_da_divisao_deterministica_dao_acaso():
    """A divisao que o desenho fixou ANTES de medir: ordenar os 222 compostos
    rotulados por `inchikey_parental` e alternar. Sem gerador aleatorio, sem
    semente -- reproduzivel a partir do arquivo.

    E o que dissolve a tabela inteira: uma metade nao acompanha a outra, e as
    duas ficam no acaso."""
    (auc_a, n_a), (auc_b, n_b) = _auc_das_metades()
    assert (round(auc_a, 4), n_a) == (0.5000, 29)
    assert (round(auc_b, 4), n_b) == (0.5258, 41)
