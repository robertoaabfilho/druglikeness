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
