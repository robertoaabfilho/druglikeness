# -*- coding: utf-8 -*-
"""O carregador do preditor de alvos (PP4CKD).

Tres coisas que estes testes prendem, e que nao sao detalhe:

1. A arquitetura e os indices do `Sequential`, porque e por eles que o
   `state_dict` publicado encaixa sem renomear nada.
2. `weights_only=True` na carga: um `.pt` e um pickle, e carregar sem a trava
   executaria codigo arbitrario do arquivo baixado.
3. Ausencia de pesos vira mensagem com remedio, nunca traceback -- e nunca
   `sys.exit`, que mataria o painel e o pytest em vez de degradar.
"""
import pickle

import pytest

torch = pytest.importorskip("torch")

from mods import modelo
from tests.conftest import N_ALVOS_FALSO


def test_arquitetura_e_a_do_pp4ckd():
    """4096 -> 1000 -> 500 -> N. O Dropout e Identity porque este modulo e de
    inferencia, mas os INDICES do Sequential sao preservados: o state_dict
    original tem chaves '0.weight', '3.weight' e '6.weight'."""
    rede = modelo.construir_rede(4096, N_ALVOS_FALSO)
    assert rede[0].in_features == 4096
    assert rede[0].out_features == 1000
    assert rede[3].in_features == 1000
    assert rede[3].out_features == 500
    assert rede[6].in_features == 500
    assert rede[6].out_features == N_ALVOS_FALSO
    assert isinstance(rede[2], torch.nn.Identity)
    assert isinstance(rede[5], torch.nn.Identity)


def test_ecfp4_e_morgan_raio_2_com_4096_bits():
    cfg = modelo.FINGERPRINTS["ECFP4"]
    assert (cfg["raio"], cfg["bits"]) == (2, 4096)


def test_n_alvos_do_pp4ckd():
    assert modelo.N_ALVOS == 7676


def test_carregar_modelo_le_o_checkpoint_e_deduz_n_alvos_dele(modelo_falso):
    """O numero de alvos vem do PROPRIO checkpoint, nao de uma constante: e
    isso que permite um fold sintetico de 6 alvos e os 7.676 reais."""
    rede, cfg, dispositivo, nome = modelo.carregar_modelo(
        None, checkpoint=modelo_falso, dispositivo="cpu")
    assert rede[6].out_features == N_ALVOS_FALSO
    assert cfg["bits"] == 4096
    assert dispositivo == "cpu"
    assert nome == "ecfp4_falso.pt"


def test_cpu_e_respeitado_quando_pedido(modelo_falso):
    """GPU nao reproduz: 36 de 312 colunas diferem entre duas execucoes de GPU,
    0 de 312 na CPU. Quem pede CPU recebe CPU."""
    _, _, dispositivo, _ = modelo.carregar_modelo(
        None, checkpoint=modelo_falso, dispositivo="cpu")
    assert dispositivo == "cpu"


def test_checkpoint_inexistente_nomeia_o_caminho(tmp_path):
    with pytest.raises(FileNotFoundError) as e:
        modelo.carregar_modelo(None, checkpoint=str(tmp_path / "nao_existe.pt"))
    assert "nao_existe.pt" in str(e.value)


def test_pasta_sem_pesos_diz_o_que_falta_e_onde_baixar(tmp_path):
    """Review Focus 2: `modelo/` ausente ou vazio da mensagem com remedio."""
    with pytest.raises(FileNotFoundError) as e:
        modelo.carregar_modelo(str(tmp_path))
    msg = str(e.value)
    assert "ecfp4_full_model.pt" in msg
    assert "fold" in msg
    assert "zenodo" in msg.lower()


def test_pasta_sem_full_model_cai_para_o_fold_com_aviso(tmp_path, capsys):
    """O tarball do Zenodo esta truncado na origem e `ecfp4_full_model.pt` nao
    sai integro. Cair para o fold e o contorno, e ele e ANUNCIADO."""
    (tmp_path / "ecfp4_fold_4.pt").write_bytes(b"")
    escolhido = modelo.escolher_checkpoint(tmp_path, modelo.FINGERPRINTS["ECFP4"])
    assert escolhido.name == "ecfp4_fold_4.pt"
    assert "fold" in capsys.readouterr().out.lower()


def test_carga_recusa_pickle_que_nao_seja_de_pesos(tmp_path):
    """SEGURANCA: um .pt e um pickle. Com `weights_only=True`, um arquivo que
    carregue objeto arbitrario e RECUSADO em vez de executado."""
    malicioso = tmp_path / "malicioso.pt"
    with open(malicioso, "wb") as f:
        pickle.dump({"6.weight": print}, f, protocol=2)
    with pytest.raises(Exception) as e:
        modelo.carregar_modelo(None, checkpoint=str(malicioso), dispositivo="cpu")
    assert not isinstance(e.value, SystemExit)


def test_fingerprint_de_smiles_ilegivel_e_none():
    assert modelo.fingerprint_de("nao sou smiles",
                                 modelo.FINGERPRINTS["ECFP4"]) is None


def test_fingerprint_tem_o_tamanho_do_fingerprint():
    fp = modelo.fingerprint_de("CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
                               modelo.FINGERPRINTS["ECFP4"])
    assert fp is not None and len(fp) == 4096


def test_rotulos_normalizam_chembl_id_para_maiusculas(tmp_path):
    """Num lugar so, no ponto de leitura. Antes, 'chembl222' numa TSV dava
    `no_modelo = nao` no arquivo enriquecido enquanto a previsao pontuava o
    alvo -- as duas saidas se contradiziam sobre o mesmo alvo."""
    caminho = tmp_path / "rotulos.tsv"
    caminho.write_text(
        "target_idx\tchembl_id\tpref_name\torganism\ttarget_type\tactive_molecules\n"
        "0\tchembl222\tNET\tHomo sapiens\tSINGLE PROTEIN\t3180\n",
        encoding="utf-8")
    assert modelo.carregar_rotulos(str(caminho))[0]["chembl_id"] == "CHEMBL222"


def test_rotulos_trazem_n_ativos_porque_a_confianca_depende_dele(arquivo_de_rotulos):
    """recall@10 vai de 0,127 a 0,477 conforme o alvo sai de <50 para >=5.000
    moleculas ativas. Uma probabilidade de 0,9 sobre alvo com 88 ativos nao
    vale o mesmo que sobre um com 11.049."""
    rotulos = modelo.carregar_rotulos(arquivo_de_rotulos)
    assert rotulos[4]["n_ativos"] == 11049
    assert rotulos[5]["n_ativos"] == 88


def test_rotulos_preservam_organismo_porque_alvo_nao_humano_e_de_primeira_classe(
        arquivo_de_rotulos):
    """Dos 7.676 alvos reais, 3.536 nao sao humanos. O canal negativo da fase 3
    e deduzido justamente por organismo."""
    rotulos = modelo.carregar_rotulos(arquivo_de_rotulos)
    assert rotulos[4]["organismo"] == "Mycobacterium tuberculosis"
    assert rotulos[5]["organismo"] == "Rattus norvegicus"


def test_rotulos_indexados_por_target_idx_sem_buraco(arquivo_de_rotulos):
    """A lista e indexada pela POSICAO no modelo: o indice i tem de ser o alvo
    i da saida da rede, senao a probabilidade de um alvo leva o nome de outro."""
    rotulos = modelo.carregar_rotulos(arquivo_de_rotulos)
    assert len(rotulos) == N_ALVOS_FALSO
    assert rotulos[0]["chembl_id"] == "CHEMBL228"


def test_rotulos_inexistentes_nomeiam_o_caminho(tmp_path):
    with pytest.raises(FileNotFoundError) as e:
        modelo.carregar_rotulos(str(tmp_path / "nao_existe.tsv"))
    assert "nao_existe.tsv" in str(e.value)


def test_modelo_nao_importa_sys_logo_nao_pode_encerrar_o_processo():
    """O modulo antigo fazia `sys.exit(1)` no import quando faltava torch, o
    que mata o processo hospedeiro (painel, pytest, orquestrador) em vez de
    deixar a etapa 1 rodar em modo degradado. Sem `sys` importado, a saida
    abrupta e impossivel -- e a propriedade, nao a ausencia da palavra."""
    import ast
    arvore = ast.parse(open(modelo.__file__, encoding="utf-8").read())
    importados = {
        alias.name.split(".")[0]
        for no in ast.walk(arvore)
        if isinstance(no, ast.Import) for alias in no.names
    } | {
        (no.module or "").split(".")[0]
        for no in ast.walk(arvore) if isinstance(no, ast.ImportFrom)
    }
    assert "sys" not in importados


def test_dependencia_ausente_levanta_importerror_nao_systemexit():
    """A mensagem tem de dizer o remedio E ser capturavel."""
    import mods.modelo as m
    fonte = open(m.__file__, encoding="utf-8").read()
    assert "raise ImportError" in fonte
    assert "modo degradado" in fonte


def test_toda_cli_insere_a_raiz_no_sys_path_uma_vez_so():
    """A convencao do projeto: um script em `mods/` roda tanto como modulo
    (`python -m mods.x`) quanto por caminho (`python mods/x.py`), porque o
    orquestrador o chama por caminho. Isso exige UM insert apontando para a
    RAIZ. Antes havia dois, e o segundo apontava para `mods/` -- resquicio de
    quando o preditor era importado de forma plana."""
    import pathlib
    for cli in ("coletar_evidencias", "aplicar_portao", "ranquear",
                "instantaneo", "validar_modelo"):
        fonte = pathlib.Path("mods", f"{cli}.py").read_text(encoding="utf-8")
        inserts = [l.strip() for l in fonte.splitlines()
                   if "sys.path.insert" in l]
        assert len(inserts) == 1, f"{cli}: {inserts}"
        assert "parent.parent" in inserts[0], f"{cli}: {inserts[0]}"
