# -*- coding: utf-8 -*-
from mods.evidencias.entrada import ler_ligantes, ler_alvos, ler_sintomas


def test_ligante_sem_smiles_e_erro():
    ligantes, problemas = ler_ligantes("tests/dados/ligantes_minimo.tsv")
    assert any(p.gravidade == "erro" and p.campo == "smiles" for p in problemas)


def test_inchikey_parental_e_recalculado_quando_ausente():
    ligantes, _ = ler_ligantes("tests/dados/ligantes_minimo.tsv")
    aspirina = next(l for l in ligantes if l["nome"] == "aspirina")
    assert aspirina["inchikey_parental"] == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"


def test_ligante_repetido_pela_parental_e_deduplicado():
    # Review Focus 2: sal e base livre, nomes diferentes, mesma molécula.
    ligantes, problemas = ler_ligantes("tests/dados/ligantes_minimo.tsv")
    chaves = [l["inchikey_parental"] for l in ligantes]
    assert len(chaves) == len(set(chaves))
    assert any("duplicad" in p.mensagem for p in problemas)


def test_alvo_por_nome_sem_organismo_e_erro():
    alvos, problemas = ler_alvos("tests/dados/alvos_minimo.tsv")
    assert any(p.gravidade == "erro" and p.campo == "organismo" for p in problemas)


def test_sintomas_ignora_comentario():
    assert ler_sintomas("tests/dados/sintomas.txt") == ["neuropathic pain", "peripheral neuropathy"]


def test_deduplicacao_mantem_primeiro_nome_e_nomeia_descartado():
    ligantes, problemas = ler_ligantes("tests/dados/ligantes_minimo.tsv")
    nomes = [l["nome"] for l in ligantes]
    assert "cloridrato de duloxetina" in nomes and "duloxetina base" not in nomes
    aviso = next(p for p in problemas if "duplicad" in p.mensagem)
    assert aviso.gravidade == "aviso" and "duloxetina base" in aviso.mensagem


def test_alvos_validos_passam_e_nao_ha_no_modelo():
    alvos, _ = ler_alvos("tests/dados/alvos_minimo.tsv")
    assert [a["alvo_chembl_id"] for a in alvos] == ["CHEMBL222", ""]
    assert alvos[1]["organismo"] == "Homo sapiens"
    assert all("no_modelo" not in a for a in alvos)


def test_csv_com_virgula_nao_citada_avisa_e_descarta(tmp_path):
    f = tmp_path / "l.csv"
    f.write_text('nome,smiles\naspirina,CC(=O)Oc1ccccc1C(=O)O\nfoo, bar,C\n"ok, citado",CCO\n')
    ligantes, problemas = ler_ligantes(f)
    assert any(p.gravidade == "aviso" and ".csv" in p.mensagem for p in problemas)
    assert any(p.gravidade == "erro" and "virgula" in p.mensagem for p in problemas)
    assert [l["nome"] for l in ligantes] == ["aspirina", "ok, citado"]


def test_coluna_obrigatoria_ausente_e_erro(tmp_path):
    f = tmp_path / "a.tsv"
    f.write_text("alvo\torganismo2\nX\tY\n")
    alvos, problemas = ler_alvos(f)
    assert alvos == [] and any(p.campo == "organismo" for p in problemas)


def test_smiles_ilegivel_pelo_rdkit_e_erro_e_descartado(tmp_path):
    f = tmp_path / "l.tsv"
    f.write_text("nome\tsmiles\nlixo\tXXX\naspirina\tCC(=O)Oc1ccccc1C(=O)O\n")
    ligantes, problemas = ler_ligantes(f)
    erro = next(p for p in problemas if p.gravidade == "erro" and p.campo == "smiles")
    assert "XXX" in erro.mensagem and "ausente" not in erro.mensagem
    assert [l["nome"] for l in ligantes] == ["aspirina"]


def test_inchikey_parental_informado_divergente_avisa_e_usa_recalculado(tmp_path):
    f = tmp_path / "l.tsv"
    errada = "AAAAAAAAAAAAAA-UHFFFAOYSA-N"
    f.write_text(f"nome\tsmiles\tinchikey_parental\naspirina\tCC(=O)Oc1ccccc1C(=O)O\t{errada}\n")
    ligantes, problemas = ler_ligantes(f)
    aviso = next(p for p in problemas if p.gravidade == "aviso" and p.campo == "inchikey_parental")
    assert errada in aviso.mensagem and "BSYNRYMUTXBXSQ-UHFFFAOYSA-N" in aviso.mensagem
    assert ligantes[0]["inchikey_parental"] == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
