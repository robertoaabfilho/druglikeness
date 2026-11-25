from mods.evidencias.esquema import (
    COLUNAS_EVIDENCIA, EscritorEvidencias, ler_evidencias, formatar_atividade,
)


def test_cabecalho_tem_as_colunas_na_ordem(tmp_path):
    destino = tmp_path / "e.tsv"
    with EscritorEvidencias(destino) as esc:
        esc.escrever({"ligante": "quercetina", "alvo_chembl_id": "CHEMBL230"})
    primeira = destino.read_text(encoding="utf-8").splitlines()[0]
    assert primeira == "#" + "\t".join(COLUNAS_EVIDENCIA)


def test_nome_com_virgula_e_tabulacao_sobrevive_a_ida_e_volta(tmp_path):
    # Review Focus 1: 42,8% dos nomes têm vírgula, e tabulação quebraria a linha.
    nome = "1,8-cineol\tsin.: eucaliptol"
    destino = tmp_path / "e.tsv"
    with EscritorEvidencias(destino) as esc:
        esc.escrever({"ligante": nome, "alvo_chembl_id": "CHEMBL230"})
    linhas = list(ler_evidencias(destino))
    assert len(linhas) == 1
    assert linhas[0]["ligante"] == nome
    assert linhas[0]["alvo_chembl_id"] == "CHEMBL230"


def test_separador_decimal_e_ponto():
    assert formatar_atividade("IC50", "=", 120.5, "nM") == "IC50 = 120.5 nM"


def test_campo_ausente_vira_vazio_nao_a_palavra_none(tmp_path):
    destino = tmp_path / "e.tsv"
    with EscritorEvidencias(destino) as esc:
        esc.escrever({"ligante": "x"})
    assert "None" not in destino.read_text(encoding="utf-8")


def test_barra_invertida_literal_tambem_volta_exata(tmp_path):
    nomes = ["a\\tb", "fim\\", "x\r\ny\nz", "#comeca com cerquilha"]
    destino = tmp_path / "e.tsv"
    with EscritorEvidencias(destino) as esc:
        for n in nomes:
            esc.escrever({"ligante": n})
    assert [l["ligante"] for l in ler_evidencias(destino)] == nomes
    assert len(destino.read_text(encoding="utf-8").splitlines()) == 1 + len(nomes)
