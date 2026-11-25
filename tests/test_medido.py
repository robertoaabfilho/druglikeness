# -*- coding: utf-8 -*-
"""Evidencia medida: linhas de ensaio do ChEMBL, sem rede e sem o arquivo de 60 MB."""
from mods.evidencias.chembl import ListaTruncada, RespostaInesperada
from mods.evidencias.medido import linhas_medidas

AMOSTRA = "tests/dados/evidencias_amostra.csv"


def _lig(nome, chave):
    return {"nome": nome, "inchikey": chave, "inchikey_parental": chave}


LIGANTES = [
    _lig("composto-hts", "AAASNKNLMQBKFV-UHFFFAOYSA-N"),
    _lig("ligante-3-ensaios", "AAOVKJBEBIDNHE-UHFFFAOYSA-N"),
    _lig("pregabalin", "AYXYPKUFHZROOJ-ZETCQYMHSA-N"),
    _lig("ligante-limite", "BIXJFIJYBLJTMK-MEBBXXQBSA-N"),
    _lig("ligante-pct", "FJUKDAZEABGEIH-UHFFFAOYSA-N"),
    _lig("ligante-estereo", "AABSZBMLGKTAOK-IPPJXQLGSA-N"),
]
ALVOS = [
    {"alvo_chembl_id": "CHEMBL1", "alvo": "Microtubule-associated protein tau",
     "organismo": "Homo sapiens"},
    {"alvo_chembl_id": "CHEMBL2", "alvo": "Translocator protein",
     "organismo": "Rattus norvegicus"},
    {"alvo_chembl_id": "CHEMBL3", "alvo": "Kallikrein-7", "organismo": "Homo sapiens"},
    {"alvo_chembl_id": "CHEMBL4", "alvo": "Alpha-1A adrenergic receptor",
     "organismo": "Homo sapiens"},
    {"alvo_chembl_id": "CHEMBL5", "alvo": "Histamine H1 receptor",
     "organismo": "Homo sapiens"},
    {"alvo_chembl_id": "CHEMBL6", "alvo": "Voltage-dependent calcium channel",
     "organismo": "Homo sapiens"},
    {"alvo_chembl_id": "CHEMBL7",
     "alvo": "Peroxisome proliferator-activated receptor alpha",
     "organismo": "Homo sapiens"},
]


def test_registro_em_porcentagem_nao_preenche_pchembl():
    # Review Focus 4: 120 nM e 78% nao se ordenam; pchembl e a unica coluna comparavel.
    linhas = list(linhas_medidas(AMOSTRA, LIGANTES, ALVOS))
    pct = [l for l in linhas if l["unidade"] == "%"]
    assert pct and all(l["pchembl"] == "" for l in pct)
    assert all(l["valor"] != "" for l in pct)


def test_triagem_hts_e_marcada():
    linhas = list(linhas_medidas(AMOSTRA, LIGANTES, ALVOS))
    assert any(l["triagem_hts"] == "sim" for l in linhas)
    assert any(l["triagem_hts"] == "nao" for l in linhas)


def test_relacao_de_limite_inferior_entra_marcada():
    linhas = list(linhas_medidas(AMOSTRA, LIGANTES, ALVOS))
    assert any(l["relacao"] == ">" for l in linhas)


def test_uma_linha_por_ensaio_nao_por_par():
    linhas = list(linhas_medidas(AMOSTRA, LIGANTES, ALVOS))
    pares = [(l["ligante"], l["alvo_chembl_id"]) for l in linhas]
    assert len(pares) > len(set(pares))


def test_linha_medida_nunca_parece_previsao():
    linhas = list(linhas_medidas(AMOSTRA, LIGANTES, ALVOS))
    assert linhas
    for l in linhas:
        assert l["tipo_evidencia"] == "medido" and l["fonte"] == "ChEMBL"
        assert l["medida"] != "probabilidade"
    # pchembl so existe quando a fonte o trouxe, em ponto decimal
    assert any(l["pchembl"] == 6.67 for l in linhas)
    assert all("," not in str(l["valor"]) + str(l["pchembl"]) for l in linhas)


def test_alvo_so_casa_pelo_par_nome_organismo():
    ratos = [{"alvo_chembl_id": "CHEMBL9", "alvo": "Translocator protein",
              "organismo": "Homo sapiens"}]
    assert list(linhas_medidas(AMOSTRA, LIGANTES, ratos)) == []


def test_estereoisomero_fica_marcado_e_sem_id_do_ligante():
    linhas = [l for l in linhas_medidas(AMOSTRA, LIGANTES, ALVOS)
              if l["ligante"] == "ligante-estereo"]
    assert linhas and all("estereoisomero" in l["procedencia"] for l in linhas)
    assert all(l["ligante_chembl_id"] == "" for l in linhas)


def test_alvo_so_com_id_usa_rotulos_e_falha_de_um_item_nao_aborta():
    rotulos = [{"chembl_id": "CHEMBL3", "nome": "Kallikrein-7",
                "organismo": "Homo sapiens", "tipo": "SINGLE PROTEIN", "n_ativos": 5}]
    alvos = [{"alvo_chembl_id": "CHEMBL3", "alvo": "", "organismo": ""},
             {"alvo_chembl_id": "CHEMBL404", "alvo": "", "organismo": ""}]
    contagens = {}
    linhas = list(linhas_medidas(AMOSTRA, LIGANTES, alvos, rotulos=rotulos,
                                 contagens=contagens))
    assert linhas and {l["alvo_chembl_id"] for l in linhas} == {"CHEMBL3"}
    assert linhas[0]["tipo_alvo"] == "SINGLE PROTEIN"
    assert any("CHEMBL404" in p for p in contagens["problemas"])


def test_erro_de_resolucao_de_ligante_e_alvo_vira_problema_e_segue():
    def resolver_ligante(chave, parental):
        if chave.startswith("AABSZB"):  # so o ligante de estereo precisa de consulta
            raise ListaTruncada("truncada")
        return "CHEMBL999"

    def nome_do_alvo(chembl_id):
        raise RespostaInesperada("sem 'pref_name'")

    contagens = {}
    alvos = ALVOS + [{"alvo_chembl_id": "CHEMBLX", "alvo": "", "organismo": ""}]
    linhas = list(linhas_medidas(AMOSTRA, LIGANTES, alvos,
                                 resolver_ligante=resolver_ligante,
                                 resolver_alvo=nome_do_alvo, contagens=contagens))
    assert any(l["ligante"] == "pregabalin" for l in linhas)
    assert any(l["ligante"] == "ligante-estereo" for l in linhas)  # linha mantida
    assert len(contagens["problemas"]) >= 2


def test_correspondencia_do_ligante_e_coluna_propria():
    linhas = list(linhas_medidas(AMOSTRA, LIGANTES, ALVOS))
    por = {l["ligante"]: l["correspondencia_ligante"] for l in linhas}
    assert por["ligante-estereo"] == "estereoisomero"
    assert por["pregabalin"] == "exata"


def test_cli_grava_medido_sem_pchembl_misturado(tmp_path):
    from mods.coletar_evidencias import main
    from mods.evidencias.esquema import ler_evidencias
    rot = tmp_path / "rotulos.tsv"
    rot.write_text("target_idx\tchembl_id\tpref_name\torganism\ttarget_type\t"
                   "active_molecules\n0\tCHEMBL3\tKallikrein-7\tHomo sapiens\t"
                   "SINGLE PROTEIN\t5\n", encoding="utf-8")
    alv = tmp_path / "alvos.tsv"
    alv.write_text("alvo_chembl_id\nCHEMBL3\n", encoding="utf-8")
    lig = tmp_path / "lig.tsv"
    lig.write_text("#nome\tsmiles\nx\tCCO\n", encoding="utf-8")
    amostra = tmp_path / "amostra.csv"
    # InChIKey do etanol, trocada pela do ligante da amostra com Kallikrein-7
    from mods.evidencias.entrada import ler_ligantes
    chave = ler_ligantes(lig)[0][0]["inchikey"]
    from pathlib import Path
    cab, *corpo = Path(AMOSTRA).read_text(encoding="utf-8").splitlines()
    amostra.write_text("\n".join([cab] + [
        l.replace("BIXJFIJYBLJTMK-MEBBXXQBSA-N", chave, 1) for l in corpo
        if l.startswith("BIXJFIJYBLJTMK")]) + "\n", encoding="utf-8")
    saida = tmp_path / "e.tsv"
    rc = main(["--ligantes", str(lig), "--alvos", str(alv), "--sintomas",
               "tests/dados/sintomas.txt", "--saida", str(saida), "--sem-rede",
               "--rotulos", str(rot), "--evidencias-chembl", str(amostra)])
    assert rc == 0
    from mods.evidencias.manifesto import ler_manifesto
    par = ler_manifesto(tmp_path / "manifesto.json")["parametros"]
    assert par["rotulos_md5"] and par["evidencias_chembl_md5"]
    linhas = list(ler_evidencias(saida))
    assert linhas and all(l["tipo_evidencia"] == "medido" for l in linhas)
    assert not [l for l in linhas
                if l["pchembl"] != "" and l["medida"] == "probabilidade"]


def _csv_com(tmp_path, cabecalho, linhas=()):
    f = tmp_path / "x.csv"
    f.write_text(cabecalho + "\n" + "\n".join(linhas) + "\n", encoding="utf-8")
    return str(f)


def test_arquivo_errado_falha_nomeando_as_colunas(tmp_path):
    import pytest
    errado = _csv_com(tmp_path, "id,nome,smiles", ["1,a,CCO"])
    with pytest.raises(ValueError, match="alvo.*organismo_alvo"):
        list(linhas_medidas(errado, LIGANTES, ALVOS))


def test_cli_recusa_arquivo_errado_sem_deixar_saida(tmp_path, capsys):
    from mods.coletar_evidencias import main
    errado = _csv_com(tmp_path, "id,nome", ["1,a"])
    saida = tmp_path / "e.tsv"
    rc = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
               "--alvos", "tests/dados/alvos_minimo.tsv", "--sintomas",
               "tests/dados/sintomas.txt", "--saida", str(saida), "--sem-rede",
               "--evidencias-chembl", errado])
    assert rc == 2 and not saida.exists()
    assert "evidencia medida" in capsys.readouterr().err


def test_cli_avisa_quando_leu_linhas_e_nenhuma_casou(tmp_path, capsys):
    from mods.coletar_evidencias import main
    saida = tmp_path / "e.tsv"
    rc = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
               "--alvos", "tests/dados/alvos_minimo.tsv", "--sintomas",
               "tests/dados/sintomas.txt", "--saida", str(saida), "--sem-rede",
               "--evidencias-chembl", AMOSTRA])
    err = capsys.readouterr().err
    from mods.coletar_evidencias import SAIDA_PARCIAL
    assert rc == SAIDA_PARCIAL and "NENHUMA casou" in err


def test_correspondencia_desconhecida_ou_nfd_nao_vira_estereoisomero(tmp_path):
    import csv
    import unicodedata
    with open(AMOSTRA, encoding="utf-8") as f:
        leitor = csv.DictReader(f)
        base = next(l for l in leitor if l["InChIKey"].startswith("BIXJFIJYBL"))
        cab = leitor.fieldnames
    destino = tmp_path / "c.csv"
    with open(destino, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, cab)
        w.writeheader()
        for valor in ("", "estranha",
                      unicodedata.normalize("NFD", "estereoisômero")):
            w.writerow({**base, "correspondencia": valor})
    contagens = {}
    linhas = list(linhas_medidas(str(destino), LIGANTES, ALVOS, contagens=contagens))
    assert [l["correspondencia_ligante"] for l in linhas] == [
        "desconhecida", "desconhecida", "estereoisomero"]
    assert all("estereoisomero da" not in l["procedencia"] for l in linhas[:2])
    assert len(contagens["problemas"]) == 2


def test_inteiros_saem_como_inteiros():
    linhas = list(linhas_medidas(AMOSTRA, LIGANTES, ALVOS))
    assert all(not str(l["valor"]).endswith(".0") for l in linhas)
    assert all(not l["ano"].endswith(".0") and not l["pmid"].endswith(".0")
               for l in linhas)
    assert any(l["valor"] == 213 for l in linhas)


def test_alvo_pedido_duas_vezes_nao_dobra_a_evidencia():
    """Mesmo alvo por ID e por nome: um ensaio tem de dar UMA linha. Numa tabela
    de uma linha por evidencia, a linha repetida pareceria um segundo ensaio."""
    rotulos = [{"chembl_id": "CHEMBL3", "nome": "Kallikrein-7",
                "organismo": "Homo sapiens", "tipo": "SINGLE PROTEIN", "n_ativos": 5}]
    uma_vez = list(linhas_medidas(AMOSTRA, LIGANTES, [
        {"alvo_chembl_id": "CHEMBL3", "alvo": "", "organismo": ""}], rotulos=rotulos))
    contagens = {}
    duas_vezes = list(linhas_medidas(AMOSTRA, LIGANTES, [
        {"alvo_chembl_id": "CHEMBL3", "alvo": "", "organismo": ""},
        {"alvo_chembl_id": "", "alvo": "Kallikrein-7", "organismo": "Homo sapiens"}],
        rotulos=rotulos, contagens=contagens))
    assert uma_vez and len(duas_vezes) == len(uma_vez)
    assert any("mesmo alvo" in p for p in contagens["problemas"])


def test_origem_do_alvo_sai_na_linha_medida():
    alvos = [{**ALVOS[2], "origem_alvo": "open_targets"}]
    linhas = list(linhas_medidas(AMOSTRA, LIGANTES, alvos))
    assert linhas and all(l["origem_alvo"] == "open_targets" for l in linhas)
