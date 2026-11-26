# -*- coding: utf-8 -*-
"""O CLI da fase 3: faixas, as duas tabelas e a publicacao atomica.

Nenhum teste usa rede nem exige checkpoint ou `target_labels.tsv` reais.
"""

from pathlib import Path

from mods.portao.esquema import ler_tsv
from mods.ranquear import (NOME_MANIFESTO, NOME_NAO_RANQUEADOS, NOME_RANKING,
                           main)

COLUNAS_EVID = ["ligante", "ligante_inchikey_parental", "ligante_chembl_id",
                "alvo_chembl_id", "alvo", "organismo", "origem_alvo",
                "tipo_evidencia", "valor", "unidade", "pchembl",
                "n_moleculas_treino", "no_modelo", "procedencia"]
COLUNAS_PORTAO = ["ligante", "ligante_inchikey_parental", "aprovado",
                  "motivos", "politica"]


def _tsv(caminho: Path, colunas, linhas) -> str:
    """Escreve um TSV com `#` na primeira coluna, como o EscritorTSV faz."""
    texto = ["#" + "\t".join(colunas)]
    for linha in linhas:
        texto.append("\t".join(str(linha.get(c, "")) for c in colunas))
    caminho.write_text("\n".join(texto) + "\n", encoding="utf-8")
    return str(caminho)


def _ev(nome, inchikey, **extra):
    linha = {"ligante": nome, "ligante_inchikey_parental": inchikey,
             "alvo_chembl_id": "CHEMBL228", "alvo": "SERT",
             "organismo": "Homo sapiens", "origem_alvo": "usuario",
             "tipo_evidencia": "medido", "pchembl": "8.0", "unidade": "nM",
             "n_moleculas_treino": "900", "no_modelo": "sim",
             "procedencia": "ChEMBL"}
    linha.update(extra)
    return linha


def _evidencias(tmp_path, linhas) -> str:
    return _tsv(tmp_path / "evidencias.tsv", COLUNAS_EVID, linhas)


def _ler(caminho):
    return list(ler_tsv(caminho))


def _saida(tmp_path):
    destino = tmp_path / "fora"
    destino.mkdir()
    return destino


# --------------------------------------------------------------------------
# As tres saidas e a ordem
# --------------------------------------------------------------------------

def test_ordena_por_faixa_e_publica_as_tres_saidas(tmp_path):
    evid = _evidencias(tmp_path, [
        _ev("previsto_enriq", "K" * 27, origem_alvo="reactome",
            tipo_evidencia="predito", valor="0.99", pchembl=""),
        _ev("medido_usuario", "A" * 27, pchembl="8.0"),
        _ev("previsto_usuario", "B" * 27, tipo_evidencia="predito",
            valor="0.50", pchembl=""),
        _ev("medido_enriq", "C" * 27, origem_alvo="open_targets", pchembl="5.0"),
    ])
    destino = _saida(tmp_path)
    assert main(["--evidencias", evid, "--saida", str(destino),
                 "--sem-canal-negativo"]) == 0
    assert {p.name for p in destino.iterdir()} == {
        NOME_RANKING, NOME_NAO_RANQUEADOS, NOME_MANIFESTO}
    assert [l["ligante"] for l in _ler(destino / NOME_RANKING)] == [
        "medido_usuario", "medido_enriq", "previsto_usuario",
        "previsto_enriq"]


def test_a_posicao_comeca_em_1_e_nao_tem_buraco(tmp_path):
    evid = _evidencias(tmp_path, [_ev("a", "A" * 27), _ev("b", "B" * 27)])
    destino = _saida(tmp_path)
    main(["--evidencias", evid, "--saida", str(destino), "--sem-canal-negativo"])
    assert [l["posicao"] for l in _ler(destino / NOME_RANKING)] == ["1", "2"]


# --------------------------------------------------------------------------
# Quem nao ranqueia
# --------------------------------------------------------------------------

def test_faixa_3_sai_em_nao_ranqueados_sem_posicao(tmp_path):
    # alvo que ninguem pediu nem o enriquecimento trouxe: nao e afirmativo
    evid = _evidencias(tmp_path, [_ev("orfao", "O" * 27, origem_alvo="varredura")])
    destino = _saida(tmp_path)
    main(["--evidencias", evid, "--saida", str(destino), "--sem-canal-negativo"])
    fora = _ler(destino / NOME_NAO_RANQUEADOS)
    assert [l["ligante"] for l in fora] == ["orfao"]
    assert fora[0]["motivo"] == "sem evidencia em alvo afirmativo"
    assert "posicao" not in fora[0] and "#posicao" not in fora[0]
    assert not _ler(destino / NOME_RANKING)


def test_reprovado_no_portao_nao_ranqueia_mas_sai_com_o_criterio(tmp_path):
    evid = _evidencias(tmp_path, [_ev("cipro", "R" * 27)])
    portao = _tsv(tmp_path / "portao.tsv", COLUNAS_PORTAO, [
        {"ligante": "cipro", "ligante_inchikey_parental": "R" * 27,
         "aprovado": "nao", "motivos": "alerta_bula=neuropathy@boxed_warning",
         "politica": "alerta_bula=reprovar"}])
    destino = _saida(tmp_path)
    main(["--evidencias", evid, "--saida", str(destino), "--portao", portao,
          "--sem-canal-negativo"])
    fora = _ler(destino / NOME_NAO_RANQUEADOS)
    assert not _ler(destino / NOME_RANKING)
    assert fora[0]["motivo"].startswith("reprovado no portao")
    assert "alerta_bula" in fora[0]["motivo"]


def test_reprovado_e_faixa_3_aparece_uma_vez_so_com_os_dois_motivos(tmp_path):
    evid = _evidencias(tmp_path, [_ev("duplo", "D" * 27, origem_alvo="varredura")])
    portao = _tsv(tmp_path / "portao.tsv", COLUNAS_PORTAO, [
        {"ligante": "duplo", "ligante_inchikey_parental": "D" * 27,
         "aprovado": "nao", "motivos": "smiles=ilegivel", "politica": "x"}])
    destino = _saida(tmp_path)
    main(["--evidencias", evid, "--saida", str(destino), "--portao", portao,
          "--sem-canal-negativo"])
    fora = _ler(destino / NOME_NAO_RANQUEADOS)
    assert len(fora) == 1
    assert "reprovado no portao" in fora[0]["motivo"]
    assert "sem evidencia em alvo afirmativo" in fora[0]["motivo"]


# --------------------------------------------------------------------------
# O portao carregado
# --------------------------------------------------------------------------

def test_marcado_no_portao_ranqueia_carregando_as_marcas(tmp_path):
    # a etapa 2 marca porque reprovar apaga os positivos; o ranking nao desfaz
    evid = _evidencias(tmp_path, [_ev("amitriptilina", "M" * 27)])
    portao = _tsv(tmp_path / "portao.tsv", COLUNAS_PORTAO, [
        {"ligante": "amitriptilina", "ligante_inchikey_parental": "M" * 27,
         "aprovado": "sim", "motivos": "", "politica": "x"}])
    destino = _saida(tmp_path)
    main(["--evidencias", evid, "--saida", str(destino), "--portao", portao,
          "--sem-canal-negativo"])
    linha = _ler(destino / NOME_RANKING)[0]
    assert linha["portao_aprovado"] == "sim" and linha["portao_motivos"] == ""


def test_sem_portao_a_coluna_diz_que_nao_foi_consultado(tmp_path):
    # nunca `aprovado = sim` por omissao
    evid = _evidencias(tmp_path, [_ev("x", "X" * 27)])
    destino = _saida(tmp_path)
    main(["--evidencias", evid, "--saida", str(destino), "--sem-canal-negativo"])
    assert _ler(destino / NOME_RANKING)[0]["portao_aprovado"] == "nao consultado"


def test_composto_ausente_do_portao_nao_vira_aprovado(tmp_path):
    evid = _evidencias(tmp_path, [_ev("fora_do_portao", "F" * 27)])
    portao = _tsv(tmp_path / "portao.tsv", COLUNAS_PORTAO, [])
    destino = _saida(tmp_path)
    main(["--evidencias", evid, "--saida", str(destino), "--portao", portao,
          "--sem-canal-negativo"])
    assert _ler(destino / NOME_RANKING)[0]["portao_aprovado"] == "ausente"


# --------------------------------------------------------------------------
# Repetibilidade e validacao
# --------------------------------------------------------------------------

def test_empate_de_valor_ordena_estavel_pelo_inchikey(tmp_path):
    evid = _evidencias(tmp_path, [_ev("z", "Z" * 27), _ev("a", "A" * 27)])
    saidas = []
    for n in ("um", "dois"):
        destino = tmp_path / n
        destino.mkdir()
        main(["--evidencias", evid, "--saida", str(destino),
              "--sem-canal-negativo"])
        saidas.append((destino / NOME_RANKING).read_bytes())
    assert saidas[0] == saidas[1]
    assert [l["ligante"] for l in _ler(tmp_path / "um" / NOME_RANKING)] == ["a", "z"]


def test_saida_que_aponta_para_arquivo_existente_e_erro_antes_de_escrever(tmp_path):
    evid = _evidencias(tmp_path, [_ev("a", "A" * 27)])
    arquivo = tmp_path / "ocupado"
    arquivo.write_text("nao sou pasta\n", encoding="utf-8")
    assert main(["--evidencias", evid, "--saida", str(arquivo)]) == 2
    assert arquivo.read_text(encoding="utf-8") == "nao sou pasta\n"


def test_evidencias_inexistente_e_erro(tmp_path):
    assert main(["--evidencias", str(tmp_path / "nao_existe.tsv"),
                 "--saida", str(_saida(tmp_path))]) == 2


def test_nenhum_composto_utilizavel_sai_1_e_nao_escreve_nada(tmp_path):
    vazio = _evidencias(tmp_path, [])
    destino = _saida(tmp_path)
    assert main(["--evidencias", vazio, "--saida", str(destino)]) == 1
    assert not list(destino.iterdir())


def test_limite_usa_so_os_n_primeiros_compostos(tmp_path):
    evid = _evidencias(tmp_path, [_ev("a", "A" * 27), _ev("b", "B" * 27),
                                  _ev("c", "C" * 27)])
    destino = _saida(tmp_path)
    main(["--evidencias", evid, "--saida", str(destino), "--limite", "2",
          "--sem-canal-negativo"])
    assert len(_ler(destino / NOME_RANKING)) == 2


def test_manifesto_registra_as_contagens_e_a_procedencia(tmp_path):
    import json
    evid = _evidencias(tmp_path, [_ev("a", "A" * 27),
                                  _ev("o", "O" * 27, origem_alvo="varredura")])
    destino = _saida(tmp_path)
    main(["--evidencias", evid, "--saida", str(destino), "--sem-canal-negativo"])
    m = json.loads((destino / NOME_MANIFESTO).read_text(encoding="utf-8"))
    assert m["parametros"]["canal_negativo"] is False
    assert m["contagens"]["ranqueados"] == 1
    assert m["contagens"]["nao_ranqueados"] == 1


# --------------------------------------------------------------------------
# O canal negativo
# --------------------------------------------------------------------------

import json

import numpy as np
import pytest

from mods import ranquear
from mods.prever_alvos import carregar_rotulos

ROTULOS_FALSOS = carregar_rotulos("tests/dados/rotulos_minimo.tsv")
IDX_NAO_HUMANO_FORTE = 1        # CHEMBL360 M. tuberculosis, 11.049 ativos


def _ligantes(tmp_path, linhas) -> str:
    return _tsv(tmp_path / "ligantes.tsv",
                ["nome", "smiles", "inchikey_parental"], linhas)


def _ckpt(tmp_path) -> str:
    """Um arquivo de pesos de mentira, estavel por tmp_path.

    Existe porque o canal negativo exige `--manifesto-etapa1`, e o guarda de
    MD5 precisa de um arquivo real para hashear. Os pesos nunca sao lidos: os
    testes substituem `_abrir_modelo`."""
    caminho = tmp_path / "falso.pt"
    if not caminho.exists():
        caminho.write_bytes(b"pesos de mentira")
    return str(caminho)


def _man_etapa1(tmp_path) -> str:
    from mods.evidencias.procedencia import md5
    caminho = tmp_path / "man_etapa1.json"
    if not caminho.exists():
        caminho.write_text(
            json.dumps({"versoes": {"checkpoint_md5": md5(_ckpt(tmp_path))}}),
            encoding="utf-8")
    return str(caminho)


def _modelo_falso(p_do_canal=0.9, dispositivo="cpu"):
    """Modelo com `prever` fixo: o canal e testado, nao a rede.

    O `cfg` e o do ECFP4 de verdade, para que `fingerprint_de` continue sendo
    quem decide se o SMILES e legivel -- e esse caminho que o teste do SMILES
    ilegivel exercita."""
    def prever(fingerprints):
        linha = np.zeros(len(ROTULOS_FALSOS), dtype=float)
        linha[IDX_NAO_HUMANO_FORTE] = p_do_canal
        return np.tile(linha, (len(fingerprints), 1))

    return ranquear.Modelo(rotulos=ROTULOS_FALSOS, cfg={"raio": 2, "bits": 2048},
                           dispositivo=dispositivo, prever=prever,
                           checkpoint="falso.pt")


def test_canal_recusa_checkpoint_diferente_do_que_a_etapa_1_registrou(
        tmp_path, capsys):
    # os dois lados de `p x (1 - p_negativo)` tem de vir do MESMO modelo
    a = tmp_path / "a.pt"
    a.write_bytes(b"pesos a")
    b = tmp_path / "b.pt"
    b.write_bytes(b"pesos b")
    from mods.evidencias.procedencia import md5
    manifesto = tmp_path / "manifesto_etapa1.json"
    manifesto.write_text(json.dumps(
        {"versoes": {"checkpoint_md5": md5(str(a))}}), encoding="utf-8")
    evid = _evidencias(tmp_path, [_ev("x", "X" * 27)])
    assert main(["--evidencias", evid, "--saida", str(_saida(tmp_path)),
                 "--checkpoint", str(b), "--rotulos", "tests/dados/rotulos_minimo.tsv",
                 "--ligantes", _ligantes(tmp_path, []),
                 "--manifesto-etapa1", str(manifesto)]) == 2
    assert "md5" in capsys.readouterr().err.lower()


def test_canal_aceita_o_mesmo_checkpoint_que_a_etapa_1_registrou(tmp_path):
    a = tmp_path / "a.pt"
    a.write_bytes(b"pesos a")
    from mods.evidencias.procedencia import md5
    manifesto = tmp_path / "manifesto_etapa1.json"
    manifesto.write_text(json.dumps(
        {"versoes": {"checkpoint_md5": md5(str(a))}}), encoding="utf-8")
    evid = _evidencias(tmp_path, [_ev("x", "X" * 27)])
    lig = _ligantes(tmp_path, [{"nome": "x", "smiles": "CCO",
                                "inchikey_parental": "X" * 27}])
    assert ranquear._erro_de_checkpoint(str(a), str(manifesto)) is None
    assert main(["--evidencias", evid, "--saida", str(_saida(tmp_path)),
                 "--checkpoint", str(a), "--rotulos",
                 "tests/dados/rotulos_minimo.tsv", "--ligantes", lig,
                 "--manifesto-etapa1", str(manifesto),
                 "--sem-canal-negativo"]) == 0


def test_canal_sem_pesos_registra_falso_e_nao_usa_zero(tmp_path):
    evid = _evidencias(tmp_path, [_ev("x", "X" * 27)])
    destino = _saida(tmp_path)
    assert main(["--evidencias", evid, "--saida", str(destino)]) == 0
    m = json.loads((destino / NOME_MANIFESTO).read_text(encoding="utf-8"))
    assert m["parametros"]["canal_negativo"] is False
    assert _ler(destino / NOME_RANKING)[0]["p_negativo"] == "indeterminado"


def test_canal_exige_ligantes_porque_precisa_do_smiles(tmp_path, capsys):
    evid = _evidencias(tmp_path, [_ev("x", "X" * 27)])
    assert main(["--evidencias", evid, "--saida", str(_saida(tmp_path)),
                 "--checkpoint", str(tmp_path / "falso.pt"),
                 "--rotulos", "tests/dados/rotulos_minimo.tsv"]) == 2
    assert "--ligantes" in capsys.readouterr().err


def test_canal_smiles_ilegivel_da_indeterminado_e_ranqueia_pelo_resto(
        tmp_path, monkeypatch):
    monkeypatch.setattr(ranquear, "_abrir_modelo", lambda args: _modelo_falso())
    evid = _evidencias(tmp_path, [
        _ev("ilegivel", "I" * 27, tipo_evidencia="predito", valor="0.5",
            pchembl=""),
        _ev("legivel", "L" * 27, tipo_evidencia="predito", valor="0.5",
            pchembl="")])
    lig = _ligantes(tmp_path, [
        {"nome": "ilegivel", "smiles": "C1CC", "inchikey_parental": "I" * 27},
        {"nome": "legivel", "smiles": "CCO", "inchikey_parental": "L" * 27}])
    destino = _saida(tmp_path)
    assert main(["--evidencias", evid, "--saida", str(destino),
                 "--checkpoint", _ckpt(tmp_path), "--manifesto-etapa1",
                 _man_etapa1(tmp_path), "--rotulos",
                 "tests/dados/rotulos_minimo.tsv", "--ligantes", lig]) == 3
    por_nome = {l["ligante"]: l for l in _ler(destino / NOME_RANKING)}
    assert por_nome["ilegivel"]["p_negativo"] == "indeterminado"
    assert por_nome["ilegivel"]["posicao"]
    assert por_nome["legivel"]["p_negativo"] == "0.9"


def test_canal_multiplica_so_a_faixa_prevista(tmp_path, monkeypatch):
    monkeypatch.setattr(ranquear, "_abrir_modelo", lambda args: _modelo_falso(0.9))
    evid = _evidencias(tmp_path, [
        _ev("medido", "M" * 27, pchembl="8.0"),
        _ev("previsto", "P" * 27, tipo_evidencia="predito", valor="0.90",
            pchembl="")])
    lig = _ligantes(tmp_path, [
        {"nome": "medido", "smiles": "CCO", "inchikey_parental": "M" * 27},
        {"nome": "previsto", "smiles": "CCN", "inchikey_parental": "P" * 27}])
    destino = _saida(tmp_path)
    main(["--evidencias", evid, "--saida", str(destino), "--checkpoint",
          _ckpt(tmp_path), "--manifesto-etapa1", _man_etapa1(tmp_path),
          "--rotulos", "tests/dados/rotulos_minimo.tsv",
          "--ligantes", lig])
    por_nome = {l["ligante"]: l for l in _ler(destino / NOME_RANKING)}
    assert float(por_nome["medido"]["valor"]) == pytest.approx(8.0)
    assert float(por_nome["previsto"]["valor"]) == pytest.approx(0.09)


def test_canal_nomeia_o_alvo_que_produziu_o_maximo(tmp_path, monkeypatch):
    monkeypatch.setattr(ranquear, "_abrir_modelo", lambda args: _modelo_falso())
    evid = _evidencias(tmp_path, [_ev("x", "X" * 27)])
    lig = _ligantes(tmp_path, [{"nome": "x", "smiles": "CCO",
                                "inchikey_parental": "X" * 27}])
    destino = _saida(tmp_path)
    main(["--evidencias", evid, "--saida", str(destino), "--checkpoint",
          _ckpt(tmp_path), "--manifesto-etapa1", _man_etapa1(tmp_path),
          "--rotulos", "tests/dados/rotulos_minimo.tsv",
          "--ligantes", lig])
    alvo = _ler(destino / NOME_RANKING)[0]["alvo_negativo"]
    assert "CHEMBL360" in alvo and "Mycobacterium tuberculosis" in alvo


def test_canal_registra_o_dispositivo_que_rodou_de_fato(tmp_path, monkeypatch):
    monkeypatch.setattr(ranquear, "_abrir_modelo",
                        lambda args: _modelo_falso(dispositivo="cpu"))
    evid = _evidencias(tmp_path, [_ev("x", "X" * 27)])
    lig = _ligantes(tmp_path, [{"nome": "x", "smiles": "CCO",
                                "inchikey_parental": "X" * 27}])
    destino = _saida(tmp_path)
    main(["--evidencias", evid, "--saida", str(destino), "--checkpoint",
          _ckpt(tmp_path), "--manifesto-etapa1", _man_etapa1(tmp_path),
          "--rotulos", "tests/dados/rotulos_minimo.tsv",
          "--ligantes", lig])
    m = json.loads((destino / NOME_MANIFESTO).read_text(encoding="utf-8"))
    assert m["parametros"]["dispositivo"] == "cpu"
    assert m["parametros"]["canal_negativo"] is True


def test_canal_minimo_treino_muda_quem_entra(tmp_path, monkeypatch):
    monkeypatch.setattr(ranquear, "_abrir_modelo", lambda args: _modelo_falso())
    evid = _evidencias(tmp_path, [_ev("x", "X" * 27)])
    lig = _ligantes(tmp_path, [{"nome": "x", "smiles": "CCO",
                                "inchikey_parental": "X" * 27}])
    destino = _saida(tmp_path)
    main(["--evidencias", evid, "--saida", str(destino), "--checkpoint",
          _ckpt(tmp_path), "--manifesto-etapa1", _man_etapa1(tmp_path),
          "--rotulos", "tests/dados/rotulos_minimo.tsv",
          "--ligantes", lig, "--minimo-treino", "20000"])
    # nenhum alvo passa do minimo: mascara vazia, e isso e indeterminado
    assert _ler(destino / NOME_RANKING)[0]["p_negativo"] == "indeterminado"


# --------------------------------------------------------------------------
# Manifesto e repeticao
# --------------------------------------------------------------------------

NOME_MANIFESTO_REPETICAO = "manifesto_repeticao.json"


def _execucao_base(tmp_path, nome="base"):
    evid = _evidencias(tmp_path, [_ev("a", "A" * 27), _ev("b", "B" * 27,
                                                          pchembl="7.0")])
    destino = tmp_path / nome
    destino.mkdir()
    assert main(["--evidencias", evid, "--saida", str(destino),
                 "--sem-canal-negativo"]) == 0
    return evid, destino


def test_repeticao_da_o_mesmo_arquivo_byte_a_byte(tmp_path):
    evid, base = _execucao_base(tmp_path)
    repetida = tmp_path / "repetida"
    repetida.mkdir()
    assert main(["--evidencias", evid, "--saida", str(repetida),
                 "--sem-canal-negativo",
                 "--manifesto", str(base / NOME_MANIFESTO)]) == 0
    assert (base / NOME_RANKING).read_bytes() == \
        (repetida / NOME_RANKING).read_bytes()
    assert (base / NOME_NAO_RANQUEADOS).read_bytes() == \
        (repetida / NOME_NAO_RANQUEADOS).read_bytes()


def test_repeticao_escreve_manifesto_proprio_e_nao_sobrescreve_o_original(
        tmp_path):
    evid, base = _execucao_base(tmp_path)
    antes = (base / NOME_MANIFESTO).read_bytes()
    assert main(["--evidencias", evid, "--saida", str(base),
                 "--sem-canal-negativo",
                 "--manifesto", str(base / NOME_MANIFESTO)]) == 0
    assert (base / NOME_MANIFESTO).read_bytes() == antes
    assert (base / NOME_MANIFESTO_REPETICAO).is_file()


def test_repeticao_recusa_se_o_md5_das_evidencias_mudou(tmp_path, capsys):
    evid, base = _execucao_base(tmp_path)
    Path(evid).write_text(
        Path(evid).read_text(encoding="utf-8").replace("8.0", "9.0"),
        encoding="utf-8")
    repetida = tmp_path / "r"
    repetida.mkdir()
    assert main(["--evidencias", evid, "--saida", str(repetida),
                 "--sem-canal-negativo",
                 "--manifesto", str(base / NOME_MANIFESTO)]) == 2
    assert "evidencias_md5" in capsys.readouterr().err
    assert not list(repetida.iterdir())


def test_repeticao_recusa_se_o_minimo_treino_mudou(tmp_path, capsys):
    evid, base = _execucao_base(tmp_path)
    repetida = tmp_path / "r"
    repetida.mkdir()
    assert main(["--evidencias", evid, "--saida", str(repetida),
                 "--sem-canal-negativo", "--minimo-treino", "10",
                 "--manifesto", str(base / NOME_MANIFESTO)]) == 2
    assert "minimo_treino" in capsys.readouterr().err


def test_repeticao_recusa_se_o_canal_negativo_mudou(tmp_path, capsys,
                                                    monkeypatch):
    monkeypatch.setattr(ranquear, "_abrir_modelo", lambda args: _modelo_falso())
    evid = _evidencias(tmp_path, [_ev("a", "A" * 27)])
    lig = _ligantes(tmp_path, [{"nome": "a", "smiles": "CCO",
                                "inchikey_parental": "A" * 27}])
    base = tmp_path / "base"
    base.mkdir()
    assert main(["--evidencias", evid, "--saida", str(base), "--checkpoint",
                 _ckpt(tmp_path), "--manifesto-etapa1", _man_etapa1(tmp_path),
                 "--rotulos", "tests/dados/rotulos_minimo.tsv",
                 "--ligantes", lig]) == 0
    repetida = tmp_path / "r"
    repetida.mkdir()
    # `--ligantes` segue igual de proposito: sem ele a recusa viria por
    # `ligantes_md5`, e o teste nao provaria nada sobre o canal.
    assert main(["--evidencias", evid, "--saida", str(repetida),
                 "--ligantes", lig, "--sem-canal-negativo",
                 "--manifesto", str(base / NOME_MANIFESTO)]) == 2
    saida = capsys.readouterr().err
    # Desligar o canal muda tambem as chaves de VERSAO (pesos, rotulos, torch),
    # e a primeira passada tropeca nelas antes de chegar a `canal_negativo`.
    # Registrado como esta: `canal_negativo` e redundante na pratica, e o teste
    # diz a verdade sobre QUEM recusa, em vez de afirmar o que nao acontece.
    assert "pp4ckd_md5" in saida and "rotulos_md5" in saida
    assert not list(repetida.iterdir())


def test_manifesto_registra_as_chaves_da_repeticao(tmp_path):
    _, base = _execucao_base(tmp_path)
    m = json.loads((base / NOME_MANIFESTO).read_text(encoding="utf-8"))
    assert set(m["parametros"]) >= {"canal_negativo", "minimo_treino",
                                    "dispositivo", "lote", "limite", "portao"}
    assert set(m["versoes"]) >= {"rdkit", "evidencias_md5", "pp4ckd_md5",
                                 "portao_md5", "alvos_md5"}


def test_manifesto_inexistente_e_erro_antes_de_escrever(tmp_path, capsys):
    evid = _evidencias(tmp_path, [_ev("a", "A" * 27)])
    destino = _saida(tmp_path)
    assert main(["--evidencias", evid, "--saida", str(destino),
                 "--manifesto", str(tmp_path / "nao_existe.json")]) == 2
    assert not list(destino.iterdir())


def test_ligante_sem_evidencia_nenhuma_e_listado_em_vez_de_desaparecer(
        tmp_path):
    """Nada desaparece em silencio.

    Quem nao tem linha em `evidencias.tsv` nunca chegava aqui: a etapa 1 so
    escreve linhas de evidencia, logo um ligante sem evidencia alguma sumia --
    nao ranqueado e nao listado, o que e o caso MAIS extremo de faixa 3."""
    evid = _evidencias(tmp_path, [_ev("com_evidencia", "C" * 27)])
    lig = _ligantes(tmp_path, [
        {"nome": "com_evidencia", "smiles": "CCO",
         "inchikey_parental": "C" * 27},
        {"nome": "sem_evidencia", "smiles": "CCN",
         "inchikey_parental": "S" * 27}])
    destino = _saida(tmp_path)
    assert main(["--evidencias", evid, "--saida", str(destino),
                 "--ligantes", lig, "--sem-canal-negativo"]) == 0
    fora = {l["ligante"]: l for l in _ler(destino / NOME_NAO_RANQUEADOS)}
    assert "sem_evidencia" in fora
    assert fora["sem_evidencia"]["motivo"] == "sem evidencia na etapa 1"
    assert fora["sem_evidencia"]["faixa"] == "3"
    assert [l["ligante"] for l in _ler(destino / NOME_RANKING)] == \
        ["com_evidencia"]


def test_sem_ligantes_o_que_nao_esta_nas_evidencias_nao_e_inventado(tmp_path):
    # Sem a lista de ligantes nao ha como saber quem ficou de fora, e o
    # manifesto nao deve afirmar uma contagem que a execucao nao pode conhecer.
    evid = _evidencias(tmp_path, [_ev("a", "A" * 27)])
    destino = _saida(tmp_path)
    main(["--evidencias", evid, "--saida", str(destino), "--sem-canal-negativo"])
    assert not _ler(destino / NOME_NAO_RANQUEADOS)


def test_ligante_sem_evidencia_e_reprovado_traz_os_dois_motivos(tmp_path):
    """O motivo nao pode depender do CAMINHO pelo qual o composto saiu.

    Quem esta nas evidencias e foi reprovado ganha "reprovado no portao" no
    motivo; quem nao esta nas evidencias ganhava so "sem evidencia", e a
    reprovacao ficava apenas na coluna. Dois leitores da mesma tabela tirariam
    conclusoes diferentes sobre o mesmo fato."""
    evid = _evidencias(tmp_path, [_ev("outro", "O" * 27)])
    lig = _ligantes(tmp_path, [
        {"nome": "outro", "smiles": "CCO", "inchikey_parental": "O" * 27},
        {"nome": "cipro", "smiles": "CCN", "inchikey_parental": "R" * 27}])
    portao = _tsv(tmp_path / "portao.tsv", COLUNAS_PORTAO, [
        {"ligante": "cipro", "ligante_inchikey_parental": "R" * 27,
         "aprovado": "nao", "motivos": "alerta_bula: neuropathy",
         "politica": "x"}])
    destino = _saida(tmp_path)
    main(["--evidencias", evid, "--saida", str(destino), "--ligantes", lig,
          "--portao", portao, "--sem-canal-negativo"])
    fora = {l["ligante"]: l for l in _ler(destino / NOME_NAO_RANQUEADOS)}
    motivo = fora["cipro"]["motivo"]
    assert "reprovado no portao" in motivo
    assert "sem evidencia na etapa 1" in motivo


# --------------------------------------------------------------------------
# Guardas que a revisao de contexto fresco encontrou
# --------------------------------------------------------------------------

def _manifesto_etapa1(tmp_path, checkpoint) -> str:
    from mods.evidencias.procedencia import md5
    caminho = tmp_path / "manifesto_etapa1.json"
    caminho.write_text(json.dumps(
        {"versoes": {"checkpoint_md5": md5(str(checkpoint))}}),
        encoding="utf-8")
    return str(caminho)


def test_rotulos_entram_na_chave_da_repeticao(tmp_path):
    """Os rotulos decidem organismo, n_ativos e pref_name -- a mascara
    inteira. Sem `rotulos_md5` a repeticao aceitava outro espaco de alvos e
    publicava outro escore como se fosse o mesmo."""
    _, base = _execucao_base(tmp_path)
    m = json.loads((base / NOME_MANIFESTO).read_text(encoding="utf-8"))
    assert "rotulos_md5" in m["versoes"]
    assert "rotulos_md5" in ranquear.VERSOES_DA_REPETICAO


def test_repeticao_recusa_se_os_rotulos_mudaram(tmp_path, capsys, monkeypatch):
    monkeypatch.setattr(ranquear, "_abrir_modelo", lambda args: _modelo_falso())
    rot1 = tmp_path / "rot1.tsv"
    rot1.write_text(Path("tests/dados/rotulos_minimo.tsv").read_text(
        encoding="utf-8"), encoding="utf-8")
    evid = _evidencias(tmp_path, [_ev("a", "A" * 27)])
    lig = _ligantes(tmp_path, [{"nome": "a", "smiles": "CCO",
                                "inchikey_parental": "A" * 27}])
    ckpt = tmp_path / "c.pt"
    ckpt.write_bytes(b"pesos")
    man1 = _manifesto_etapa1(tmp_path, ckpt)
    base = tmp_path / "base"
    base.mkdir()
    assert main(["--evidencias", evid, "--saida", str(base), "--checkpoint",
                 str(ckpt), "--rotulos", str(rot1), "--ligantes", lig,
                 "--manifesto-etapa1", man1]) == 0
    rot1.write_text(rot1.read_text(encoding="utf-8").replace("11049", "12"),
                    encoding="utf-8")
    repetida = tmp_path / "r"
    repetida.mkdir()
    assert main(["--evidencias", evid, "--saida", str(repetida), "--checkpoint",
                 str(ckpt), "--rotulos", str(rot1), "--ligantes", lig,
                 "--manifesto-etapa1", man1,
                 "--manifesto", str(base / NOME_MANIFESTO)]) == 2
    assert "rotulos_md5" in capsys.readouterr().err


def test_canal_exige_o_manifesto_da_etapa_1(tmp_path, capsys):
    """Sem ele nao ha guarda nenhuma, e `p` e `p_negativo` podem vir de
    modelos diferentes com codigo 0 e nenhuma mensagem."""
    evid = _evidencias(tmp_path, [_ev("x", "X" * 27)])
    lig = _ligantes(tmp_path, [{"nome": "x", "smiles": "CCO",
                                "inchikey_parental": "X" * 27}])
    assert main(["--evidencias", evid, "--saida", str(_saida(tmp_path)),
                 "--checkpoint", _ckpt(tmp_path), "--rotulos",
                 "tests/dados/rotulos_minimo.tsv", "--ligantes", lig]) == 2
    assert "--manifesto-etapa1" in capsys.readouterr().err


def test_repeticao_recusa_manifesto_sem_os_md5_de_entrada(tmp_path, capsys):
    evid = _evidencias(tmp_path, [_ev("a", "A" * 27)])
    vazio = tmp_path / "manifesto_vazio.json"
    vazio.write_text(json.dumps({"versoes": {}, "parametros": {}}),
                     encoding="utf-8")
    destino = _saida(tmp_path)
    assert main(["--evidencias", evid, "--saida", str(destino),
                 "--sem-canal-negativo", "--manifesto", str(vazio)]) == 2
    assert "evidencias_md5" in capsys.readouterr().err
    assert not list(destino.iterdir())


def test_saida_nao_sobrescreve_o_manifesto_da_etapa_1(tmp_path, capsys):
    """O manifesto da etapa 1 e a chave de repeticao dela: publicar o da fase
    3 em cima dele apaga um artefato citavel."""
    evid = _evidencias(tmp_path, [_ev("a", "A" * 27)])
    destino = tmp_path / "etapa1"
    destino.mkdir()
    man = destino / NOME_MANIFESTO
    man.write_text(json.dumps({"versoes": {}}), encoding="utf-8")
    antes = man.read_bytes()
    assert main(["--evidencias", evid, "--saida", str(destino),
                 "--manifesto-etapa1", str(man), "--sem-canal-negativo"]) == 2
    assert man.read_bytes() == antes
    assert "manifesto" in capsys.readouterr().err.lower()


def test_saida_nao_gravavel_e_erro_e_nao_entrada_invalida(tmp_path, capsys):
    import os
    evid = _evidencias(tmp_path, [_ev("a", "A" * 27)])
    pai = tmp_path / "somente_leitura"
    pai.mkdir()
    os.chmod(pai, 0o500)
    try:
        codigo = main(["--evidencias", evid, "--saida", str(pai / "fora"),
                       "--sem-canal-negativo"])
    finally:
        os.chmod(pai, 0o700)
    assert codigo == 2          # 1 significaria "nenhuma entrada utilizavel"
    assert "saida" in capsys.readouterr().err.lower()


def test_linha_sem_forma_parental_torna_a_execucao_parcial(tmp_path):
    """Descartar linha em silencio com codigo 0 e o defeito que o codigo 3 da
    etapa 1 existe para evitar: automacao nao le stderr."""
    evid = _evidencias(tmp_path, [_ev("a", "A" * 27), _ev("sem_chave", "")])
    destino = _saida(tmp_path)
    assert main(["--evidencias", evid, "--saida", str(destino),
                 "--sem-canal-negativo"]) == 3
    m = json.loads((destino / NOME_MANIFESTO).read_text(encoding="utf-8"))
    assert m["contagens"]["linhas_sem_forma_parental"] == 1
    assert m["parametros"]["execucao_parcial"] is True


def test_smiles_ilegivel_torna_a_execucao_parcial(tmp_path, monkeypatch):
    monkeypatch.setattr(ranquear, "_abrir_modelo", lambda args: _modelo_falso())
    evid = _evidencias(tmp_path, [_ev("ilegivel", "I" * 27)])
    lig = _ligantes(tmp_path, [{"nome": "ilegivel", "smiles": "C1CC",
                                "inchikey_parental": "I" * 27}])
    ckpt = tmp_path / "c.pt"
    ckpt.write_bytes(b"pesos")
    destino = _saida(tmp_path)
    assert main(["--evidencias", evid, "--saida", str(destino), "--checkpoint",
                 str(ckpt), "--rotulos", "tests/dados/rotulos_minimo.tsv",
                 "--ligantes", lig,
                 "--manifesto-etapa1", _manifesto_etapa1(tmp_path, ckpt)]) == 3
    m = json.loads((destino / NOME_MANIFESTO).read_text(encoding="utf-8"))
    assert m["contagens"]["sem_fingerprint"] == 1


def test_alvo_pedido_nao_vira_antialvo_quando_limite_corta_o_composto(
        tmp_path, monkeypatch):
    """O conjunto afirmativo vinha das linhas INCLUIDAS, nao da lista
    declarada: com --limite, o alvo que o usuario pediu virava antialvo e a
    atividade nele passava a penalizar. E a inversao que o requisito 2 proibe.
    """
    monkeypatch.setattr(ranquear, "_abrir_modelo", lambda args: _modelo_falso())
    evid = _evidencias(tmp_path, [
        _ev("a", "A" * 27, alvo_chembl_id="CHEMBL360",
            alvo="Mycobacterium tuberculosis", organismo="M. tuberculosis"),
        _ev("b", "B" * 27)])
    lig = _ligantes(tmp_path, [
        {"nome": "a", "smiles": "CCO", "inchikey_parental": "A" * 27},
        {"nome": "b", "smiles": "CCN", "inchikey_parental": "B" * 27}])
    ckpt = tmp_path / "c.pt"
    ckpt.write_bytes(b"pesos")
    destino = _saida(tmp_path)
    main(["--evidencias", evid, "--saida", str(destino), "--checkpoint",
          str(ckpt), "--rotulos", "tests/dados/rotulos_minimo.tsv",
          "--ligantes", lig, "--limite", "1",
          "--manifesto-etapa1", _manifesto_etapa1(tmp_path, ckpt)])
    # CHEMBL360 foi pedido pelo composto "a", que --limite deixou de fora: ele
    # continua AFIRMATIVO, logo o 0,9 que o modelo da nele NAO penaliza "b".
    # Antes da correcao, p_negativo saia 0,9 e o escore de "b" era cortado a
    # um decimo por atividade no alvo que o usuario pediu.
    linha = _ler(destino / NOME_RANKING)[0]
    assert float(linha["p_negativo"]) == 0.0
    assert "CHEMBL360" not in linha["alvo_negativo"]


def test_largura_do_checkpoint_conferida_contra_os_rotulos(tmp_path, capsys):
    """O tarball do Zenodo esta truncado na origem e o README diz que trocar o
    par e o caso esperado: um descasamento tem de sair com mensagem, nao com
    IndexError no meio da agregacao."""
    assert ranquear._erro_de_largura(7676, 5).startswith("o --checkpoint")
    assert ranquear._erro_de_largura(5, 5) is None


def test_valor_pequeno_nao_e_impresso_como_zero(tmp_path):
    """Com 6 decimais fixos, p < 5e-7 saia `0`: varios compostos ordenados
    entre si mostravam o mesmo `0`, e a coluna que explica a ordem deixava de
    explica-la. A fixture do controle ja tem celulas dessa magnitude."""
    evid = _evidencias(tmp_path, [
        _ev("minusculo", "M" * 27, tipo_evidencia="predito",
            valor="0.00000008", pchembl=""),
        _ev("menor_ainda", "N" * 27, tipo_evidencia="predito",
            valor="0.00000002", pchembl="")])
    destino = _saida(tmp_path)
    main(["--evidencias", evid, "--saida", str(destino), "--sem-canal-negativo"])
    valores = [l["valor"] for l in _ler(destino / NOME_RANKING)]
    assert valores[0] != "0" and valores[1] != "0"
    assert float(valores[0]) > float(valores[1])
