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


def test_le_o_exemplo_do_readme_com_cabecalho_com_cerquilha():
    ligantes, problemas = ler_ligantes("input.exemplo.tsv")
    assert len(ligantes) >= 150, (len(ligantes), problemas[:3])


def test_le_o_controle_e_ignora_coluna_extra():
    ligantes, problemas = ler_ligantes("exemplos/controle.dor_neuropatica.tsv")
    assert len(ligantes) == 12 and not [p for p in problemas if p.gravidade == "erro"]
    assert ligantes[0]["grupo_esperado"] == "trata"


def test_le_lista_de_alvos_uma_por_linha(tmp_path):
    alvos, problemas = ler_alvos("exemplos/alvos.dor_neuropatica.txt")
    assert not problemas
    assert [a["alvo_chembl_id"] for a in alvos][:3] == [
        "CHEMBL1919", "CHEMBL3896", "CHEMBL228"]
    f = tmp_path / "a.txt"
    f.write_text("enoyl reductase | Staphylococcus aureus\nsem organismo\n",
                 encoding="utf-8")
    alvos, problemas = ler_alvos(f)
    assert alvos == [{"alvo_chembl_id": "", "alvo": "enoyl reductase",
                      "organismo": "Staphylococcus aureus",
                      "origem_alvo": "usuario"}]
    assert [p.gravidade for p in problemas] == ["erro"]


def test_alvos_em_csv_continuam_lidos_como_tabela(tmp_path):
    f = tmp_path / "alvos.csv"
    f.write_text("alvo_chembl_id,organismo\nCHEMBL222,Homo sapiens\n", encoding="utf-8")
    alvos, problemas = ler_alvos(f)
    assert [a["alvo_chembl_id"] for a in alvos] == ["CHEMBL222"]
    assert not [p for p in problemas if p.gravidade == "erro"]


def test_aspas_no_nome_sobrevivem_ao_tsv(tmp_path):
    # TSV nao tem citacao: com o dialeto padrao do csv, `"cis" isomer of X`
    # era lido como `cis isomer of X` -- o rotulo de identidade reescrito em
    # silencio, sem nenhum problema relatado.
    f = tmp_path / "l.tsv"
    f.write_text('nome\tsmiles\n"cis" isomer of X\tCCO\n', encoding="utf-8")
    ligantes, problemas = ler_ligantes(f)
    assert [l["nome"] for l in ligantes] == ['"cis" isomer of X']
    assert not [p for p in problemas if p.gravidade == "erro"]


def test_aspas_em_alvo_e_em_csv(tmp_path):
    f = tmp_path / "a.tsv"
    f.write_text('alvo\torganismo\n"noradrenaline transporter"\tHomo sapiens\n',
                 encoding="utf-8")
    alvos, _ = ler_alvos(f)
    assert alvos[0]["alvo"] == '"noradrenaline transporter"'
    # no .csv a citacao continua valendo: e dela que a virgula do nome depende
    g = tmp_path / "a.csv"
    g.write_text('alvo,organismo\n"transporter, sodium",Homo sapiens\n',
                 encoding="utf-8")
    alvos, _ = ler_alvos(g)
    assert alvos[0]["alvo"] == "transporter, sodium"


def test_id_do_chembl_e_normalizado_na_leitura(tmp_path):
    f = tmp_path / "a.tsv"
    f.write_text("alvo_chembl_id\nchembl222\n", encoding="utf-8")
    alvos, _ = ler_alvos(f)
    assert [a["alvo_chembl_id"] for a in alvos] == ["CHEMBL222"]


def test_alvo_repetido_e_deduplicado_e_relatado(tmp_path):
    f = tmp_path / "a.tsv"
    f.write_text("alvo_chembl_id\talvo\torganismo\n"
                 "CHEMBL222\t\t\n"
                 "chembl222\t\t\n"
                 "\tNET\tHomo sapiens\n"
                 "\tnet\thomo sapiens\n", encoding="utf-8")
    alvos, problemas = ler_alvos(f)
    assert [(a["alvo_chembl_id"], a["alvo"]) for a in alvos] == [
        ("CHEMBL222", ""), ("", "NET")]
    repetidos = [p for p in problemas if "repetido" in p.mensagem]
    assert len(repetidos) == 2 and all(p.gravidade == "aviso" for p in repetidos)


def test_lista_um_por_linha_tambem_deduplica(tmp_path):
    f = tmp_path / "a.txt"
    f.write_text("CHEMBL222\nchembl222  # o mesmo alvo\n", encoding="utf-8")
    alvos, problemas = ler_alvos(f)
    assert [a["alvo_chembl_id"] for a in alvos] == ["CHEMBL222"]
    assert any("repetido" in p.mensagem for p in problemas)


def test_origem_do_alvo_vem_da_lista_enriquecida(tmp_path):
    # O alvos_enriquecidos.tsv grava a procedencia na coluna `origem`; realimentar
    # esse arquivo como --alvos tem de preservar essa procedencia.
    f = tmp_path / "enr.tsv"
    f.write_text("#alvo_chembl_id\talvo\torganismo\torigem\n"
                 "CHEMBL222\tNET\tHomo sapiens\tusuario\n"
                 "CHEMBL233\tMOR\tHomo sapiens\topen_targets\n"
                 "CHEMBL1\tX\tHomo sapiens\treactome\n", encoding="utf-8")
    alvos, _ = ler_alvos(f)
    assert [a["origem_alvo"] for a in alvos] == ["usuario", "open_targets", "reactome"]


def test_alvo_sem_coluna_de_origem_e_do_usuario(tmp_path):
    f = tmp_path / "a.tsv"
    f.write_text("alvo_chembl_id\nCHEMBL222\n", encoding="utf-8")
    alvos, _ = ler_alvos(f)
    assert alvos[0]["origem_alvo"] == "usuario"


# --------------------------------------------------------------------------
# manter_ilegiveis: a etapa 2 julga o que a etapa 1 descarta
# --------------------------------------------------------------------------

def _tsv(tmp_path, texto):
    caminho = tmp_path / "ligantes.tsv"
    caminho.write_text(texto, encoding="utf-8")
    return caminho


LINHAS_COM_ILEGIVEL = ("nome\tsmiles\n"
                       "aspirina\tCC(=O)Oc1ccccc1C(=O)O\n"
                       "ruim\tC1CC\n")


def test_por_padrao_smiles_ilegivel_e_descartado_com_erro(tmp_path):
    """O que a etapa 1 precisa: sem molecula nao ha fingerprint nem InChIKey
    parental, logo nao ha evidencia que a linha possa preencher."""
    ligantes, problemas = ler_ligantes(_tsv(tmp_path, LINHAS_COM_ILEGIVEL))
    assert [l["nome"] for l in ligantes] == ["aspirina"]
    (problema,) = [p for p in problemas if p.campo == "smiles"]
    assert problema.gravidade == "erro"


def test_manter_ilegiveis_deixa_o_ligante_na_lista_como_aviso(tmp_path):
    """O que a etapa 2 precisa: um ligante que desaparece da tabela de
    vereditos nao foi reprovado, foi apagado."""
    ligantes, problemas = ler_ligantes(
        _tsv(tmp_path, LINHAS_COM_ILEGIVEL), manter_ilegiveis=True)
    assert [l["nome"] for l in ligantes] == ["aspirina", "ruim"]
    ruim = ligantes[1]
    assert ruim["inchikey_parental"] == "" and ruim["smiles"] == "C1CC"
    (problema,) = [p for p in problemas if p.campo == "smiles"]
    # `aviso`, nao `erro`: quem conta erros torna a execucao parcial, e aqui
    # nada se perdeu -- a diferenca entre "nao avaliei este" e "perdi este".
    assert problema.gravidade == "aviso"
    assert "mantido para julgamento" in problema.mensagem


def test_ilegiveis_nao_sao_deduplicados_entre_si(tmp_path):
    """Sem chave parental nao ha como saber se sao a mesma molecula."""
    ligantes, _ = ler_ligantes(_tsv(tmp_path, "nome\tsmiles\na\tC1CC\nb\tC1CCC\n"),
                               manter_ilegiveis=True)
    assert [l["nome"] for l in ligantes] == ["a", "b"]


def test_smiles_ausente_continua_sendo_erro_mesmo_mantendo_ilegiveis(tmp_path):
    """Sem estrutura nenhuma nao ha o que julgar nem o que relatar."""
    ligantes, problemas = ler_ligantes(
        _tsv(tmp_path, "nome\tsmiles\nsem_nada\t\n"), manter_ilegiveis=True)
    assert ligantes == []
    assert [p.gravidade for p in problemas if p.campo == "smiles"] == ["erro"]


def test_limite_conta_o_ilegivel_mantido(tmp_path):
    ligantes, _ = ler_ligantes(_tsv(tmp_path, LINHAS_COM_ILEGIVEL), 2,
                               manter_ilegiveis=True)
    assert len(ligantes) == 2
