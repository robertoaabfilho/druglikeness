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
            tipo_evidencia="previsto", valor="0.99", pchembl=""),
        _ev("medido_usuario", "A" * 27, pchembl="8.0"),
        _ev("previsto_usuario", "B" * 27, tipo_evidencia="previsto",
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
