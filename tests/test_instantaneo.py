# -*- coding: utf-8 -*-
"""O importador contra o contrato real do pipeline de prospeccao.

Medido no arquivo real (126.240 compostos, 137 MB) ao escrever estes testes:
em 30.000 linhas o importador escolhia a MESMA estrutura e a MESMA fonte que a
prospeccao, e nenhuma coluna do contrato faltava. O que ele fazia era
DESCARTAR quase tudo o que a prospeccao produz -- e nao notar se uma coluna
mudasse de nome.
"""

import csv

import pytest

from mods import instantaneo as imp
from mods import instantaneo as inst

# Uma linha com os campos que o contrato nomeia. O cafeina-como-exemplo e
# proposital: SMILES curto que o RDKit le em qualquer versao.
CAFEINA = "CN1C=NC2=C1C(=O)N(C)C(=O)N2C"
LINHA = {
    "composto_nome": "Cafeina", "nome_pubchem": "Caffeine",
    "InChIKey": "RYYVLZVUVIJVGH-UHFFFAOYSA-N", "smiles": CAFEINA,
    "fonte_smiles": "NPASS", "smiles_npass": CAFEINA, "smiles_coconut": "",
    "smiles_lotus": "", "smiles_pubchem": "", "conferencia_smiles": "",
    "chembl_id": "CHEMBL113", "correspondencia_chembl": "nome exato",
    "pubchem_cid": "2519", "composto_ids": "NPASS123",
    "plantas": "Coffea arabica", "estados_composto": "folha;semente",
    "n_plantas": "7", "atividades_relatadas_composto": "analgesica (PMID 1)",
    "melhor_evidencia_composto": "ensaio in vivo", "estudos_composto": "3",
    "n_artigos_positivos": "2", "bases": "NPASS;LOTUS",
    "atividades_duke": "Analgesic;Stimulant", "referencias_lotus": "doi:10/x",
    "atividades_chembl": "adenosine receptor",
}


def _csv(tmp_path, linhas, *, com_bom=True, colunas=None):
    """Escreve um CSV como a prospeccao escreve: UTF-8 com BOM."""
    colunas = colunas or list(LINHA)
    caminho = tmp_path / "atividades.csv"
    texto = caminho.open("w", encoding="utf-8-sig" if com_bom else "utf-8",
                         newline="")
    with texto as arq:
        escritor = csv.DictWriter(arq, fieldnames=colunas)
        escritor.writeheader()
        for linha in linhas:
            escritor.writerow({c: linha.get(c, "") for c in colunas})
    return str(caminho)


def _ler_saida(caminho):
    with open(caminho, encoding="utf-8", newline="") as arq:
        colunas = arq.readline().rstrip("\n").lstrip("#").split("\t")
        return [dict(zip(colunas, l.rstrip("\n").split("\t"))) for l in arq]


# --------------------------------------------------------------------------
# O contrato
# --------------------------------------------------------------------------

def test_contrato_completo_nao_reclama():
    assert imp.colunas_ausentes(list(LINHA)) == []


def test_coluna_do_contrato_que_falta_e_nomeada():
    colunas = [c for c in LINHA if c != "correspondencia_chembl"]
    assert imp.colunas_ausentes(colunas) == ["correspondencia_chembl"]


def test_bom_no_cabecalho_nao_esconde_a_primeira_coluna():
    """A prospeccao escreve UTF-8 com BOM, e `composto_nome` e a PRIMEIRA
    coluna: lida sem `utf-8-sig`, a chave vira '\\ufeffcomposto_nome' e o nome
    do composto cai em silencio para `nome_pubchem`."""
    # Com o BOM tratado, os dois cabecalhos tem de dar o MESMO resultado
    assert imp.colunas_ausentes(["﻿composto_nome"]) == \
        imp.colunas_ausentes(["composto_nome"])
    assert "composto_nome" not in imp.colunas_ausentes(["﻿composto_nome"])
    # e o contrato completo com BOM na primeira coluna nao reclama de nada
    com_bom = ["﻿" + list(LINHA)[0]] + list(LINHA)[1:]
    assert imp.colunas_ausentes(com_bom) == []


def test_importacao_recusa_arquivo_fora_do_contrato(tmp_path, capsys):
    origem = _csv(tmp_path, [LINHA],
                  colunas=[c for c in LINHA if c != "smiles_npass"])
    pasta = tmp_path / "input"
    assert inst.main(["--ligantes-origem", origem, "--saida", str(pasta)]) == 2
    assert "smiles_npass" in capsys.readouterr().err
    assert not (pasta / "ligantes.tsv").exists()


# --------------------------------------------------------------------------
# O que a prospeccao decidiu vale
# --------------------------------------------------------------------------

def test_usa_o_smiles_que_a_prospeccao_escolheu():
    """A prospeccao ja escolheu, comparando estruturas canonicas em quatro
    bases. Re-escolher aqui abre espaco para divergir dela em silencio se a
    regra dela mudar."""
    linha = dict(LINHA, smiles=CAFEINA, fonte_smiles="COCONUT",
                 smiles_npass="C1=CC=CC=C1")
    smiles, fonte = imp.escolher_smiles(linha)
    assert fonte == "COCONUT"
    from rdkit import Chem
    assert Chem.MolToSmiles(Chem.MolFromSmiles(smiles)) == \
        Chem.MolToSmiles(Chem.MolFromSmiles(CAFEINA))


def test_cai_para_a_ordem_de_fontes_quando_o_escolhido_e_ilegivel():
    linha = dict(LINHA, smiles="nao sou smiles", fonte_smiles="PubChem",
                 smiles_npass=CAFEINA)
    smiles, fonte = imp.escolher_smiles(linha)
    assert fonte == "NPASS" and smiles == CAFEINA


def test_fallback_de_fonte_e_registrado_na_conferencia():
    linha = dict(LINHA, smiles="nao sou smiles", fonte_smiles="PubChem",
                 smiles_npass=CAFEINA)
    convertida = imp.converter_linha(linha)
    assert "PubChem" in convertida["conferencia"]
    assert "NPASS" in convertida["conferencia"]


# --------------------------------------------------------------------------
# Nada do que a prospeccao produz se perde
# --------------------------------------------------------------------------

def test_correspondencia_chembl_e_carregada():
    """O pipeline tem o conceito (`correspondencia_ligante` em
    evidencias.tsv): um `chembl_id` casado por nome nao vale o mesmo que um
    casado por InChIKey, e jogar essa distincao fora obriga a redescobri-la."""
    assert imp.converter_linha(LINHA)["correspondencia_chembl"] == "nome exato"


def test_procedencia_da_prospeccao_e_carregada():
    convertida = imp.converter_linha(LINHA)
    assert convertida["n_plantas"] == "7"
    assert convertida["plantas"] == "Coffea arabica"
    assert convertida["bases"] == "NPASS;LOTUS"
    assert convertida["pubchem_cid"] == "2519"


def test_atividades_sao_carregadas():
    """Sao elas que permitem rotular um conjunto de validacao: foi do
    `atividades_duke` que saiu `dados/validacao_neuro_prospector.tsv`."""
    convertida = imp.converter_linha(LINHA)
    assert convertida["atividades_duke"] == "Analgesic;Stimulant"
    assert convertida["atividades_chembl"] == "adenosine receptor"
    assert convertida["melhor_evidencia_composto"] == "ensaio in vivo"
    assert convertida["n_artigos_positivos"] == "2"


def test_as_colunas_do_pipeline_vem_primeiro(tmp_path):
    """`nome` e `smiles` abrem a tabela: quem abre o arquivo no editor ve a
    identidade antes dos 12 campos de procedencia."""
    origem = _csv(tmp_path, [LINHA])
    pasta = tmp_path / "input"
    assert inst.main(["--ligantes-origem", origem, "--saida", str(pasta)]) == 0
    linhas = _ler_saida(pasta / "ligantes.tsv")
    assert list(linhas[0])[:4] == ["nome", "smiles", "inchikey",
                                   "inchikey_parental"]
    assert linhas[0]["nome"] == "Cafeina"


def test_importacao_le_o_bom_e_nao_usa_o_nome_do_pubchem(tmp_path):
    origem = _csv(tmp_path, [LINHA], com_bom=True)
    pasta = tmp_path / "input"
    inst.main(["--ligantes-origem", origem, "--saida", str(pasta)])
    assert _ler_saida(pasta / "ligantes.tsv")[0]["nome"] == "Cafeina"


def test_o_resultado_continua_legivel_pela_etapa_1(tmp_path):
    from mods.evidencias.entrada import ler_ligantes
    origem = _csv(tmp_path, [LINHA])
    pasta = tmp_path / "input"
    inst.main(["--ligantes-origem", origem, "--saida", str(pasta)])
    ligantes, problemas = ler_ligantes(str(pasta / "ligantes.tsv"))
    assert [p.mensagem for p in problemas] == []
    assert len(ligantes) == 1
    # as colunas extras chegam ao ligante em vez de serem descartadas
    assert ligantes[0]["atividades_duke"] == "Analgesic;Stimulant"


def test_limite_e_respeitado(tmp_path):
    origem = _csv(tmp_path, [LINHA, dict(LINHA, composto_nome="Outro")])
    pasta = tmp_path / "input"
    inst.main(["--ligantes-origem", origem, "--saida", str(pasta),
               "--limite", "1"])
    assert len(_ler_saida(pasta / "ligantes.tsv")) == 1


# =========================================================================
# O instantaneo: procedencia, integridade e copia das entradas de fora
# =========================================================================

def test_md5_e_estavel_e_depende_do_conteudo(tmp_path):
    a = tmp_path / "a.txt"
    b = tmp_path / "b.txt"
    a.write_text("conteudo", encoding="utf-8")
    b.write_text("conteudo", encoding="utf-8")
    assert inst.md5_de(a) == inst.md5_de(b)
    b.write_text("outro", encoding="utf-8")
    assert inst.md5_de(a) != inst.md5_de(b)


def test_linha_truncada_no_fim_do_arquivo_e_relatada_e_descartada(tmp_path):
    """Review Focus 1: a base de origem e escrita por uma execucao EM CURSO (o
    nome dela traz `_parcial`). Uma linha cortada no meio tem de ser contada e
    nomeada, nunca somada ao bolo de 'sem SMILES utilizavel' -- perder composto
    em silencio e o pior resultado possivel de um instantaneo."""
    origem = tmp_path / "parcial.csv"
    cabecalho = ",".join(inst.COLUNAS_CONTRATO)
    completa = {c: "" for c in inst.COLUNAS_CONTRATO}
    completa["composto_nome"] = "Cafeina"
    completa["smiles"] = CAFEINA
    completa["fonte_smiles"] = "NPASS"
    linha_ok = ",".join(completa[c] for c in inst.COLUNAS_CONTRATO)
    # a ultima linha para no meio: menos campos que o cabecalho
    origem.write_text(f"{cabecalho}\n{linha_ok}\nTruncada,algo\n", encoding="utf-8")

    destino = tmp_path / "ligantes.tsv"
    relatorio = inst.converter_ligantes(origem, destino)

    assert relatorio["truncadas"] == 1
    assert relatorio["gravados"] == 1
    corpo = destino.read_text(encoding="utf-8")
    assert "Cafeina" in corpo
    assert "Truncada" not in corpo


def test_arquivo_sem_as_colunas_do_contrato_e_erro_nomeando_as_que_faltam(tmp_path):
    origem = tmp_path / "errado.csv"
    origem.write_text("coluna_a,coluna_b\n1,2\n", encoding="utf-8")
    with pytest.raises(ValueError) as e:
        inst.converter_ligantes(origem, tmp_path / "saida.tsv")
    assert "composto_nome" in str(e.value)


def test_procedencia_registra_origem_data_e_md5(tmp_path):
    """Sem isto a execucao nao diz de qual estado da base de origem ela saiu, e
    o manifesto deixa de ser rastreavel ate a entrada."""
    origem = tmp_path / "origem.csv"
    origem.write_text("conteudo", encoding="utf-8")
    destino = tmp_path / "procedencia.json"

    inst.registrar([{"papel": "ligantes", "origem": origem,
                     "linhas": 3, "gravados": 2}], destino)

    import json
    dados = json.loads(destino.read_text(encoding="utf-8"))
    entrada = dados["entradas"][0]
    assert entrada["papel"] == "ligantes"
    assert entrada["origem"] == str(origem)
    assert entrada["md5"] == inst.md5_de(origem)
    assert entrada["linhas"] == 3
    assert entrada["bytes"] == len("conteudo")
    # data em ISO, para que duas execucoes sejam comparaveis
    assert dados["data"].count("-") >= 2


def test_copiar_entrada_leva_o_arquivo_e_devolve_a_procedencia(tmp_path):
    origem = tmp_path / "evidencias_chembl.csv"
    origem.write_text("molecule_chembl_id,pchembl\nCHEMBL1,6.2\n", encoding="utf-8")
    pasta = tmp_path / "input"
    pasta.mkdir()

    proc = inst.copiar_entrada(origem, pasta, papel="evidencias_chembl")

    assert (pasta / "evidencias_chembl.csv").read_text(encoding="utf-8") == \
        origem.read_text(encoding="utf-8")
    assert proc["md5"] == inst.md5_de(origem)
    assert proc["papel"] == "evidencias_chembl"


def test_copiar_entrada_recusa_sobrescrever_sem_permissao(tmp_path):
    """O instantaneo substitui a entrada da execucao anterior; fazer isso em
    silencio apagaria a lista que o usuario revisou a mao."""
    origem = tmp_path / "alvos.txt"
    origem.write_text("CHEMBL228\n", encoding="utf-8")
    pasta = tmp_path / "input"
    pasta.mkdir()
    (pasta / "alvos.txt").write_text("revisado a mao\n", encoding="utf-8")

    with pytest.raises(FileExistsError):
        inst.copiar_entrada(origem, pasta, papel="alvos")

    proc = inst.copiar_entrada(origem, pasta, papel="alvos", sobrescrever=True)
    assert proc["substituiu"] is True


def test_cli_do_instantaneo_grava_ligantes_e_procedencia(tmp_path):
    """Uma chamada, tres coisas: a tabela convertida, as copias e o registro."""
    origem = tmp_path / "compostos.csv"
    cabecalho = ",".join(inst.COLUNAS_CONTRATO)
    campos = {c: "" for c in inst.COLUNAS_CONTRATO}
    campos["composto_nome"] = "Cafeina"
    campos["smiles"] = CAFEINA
    campos["fonte_smiles"] = "NPASS"
    origem.write_text(
        f"{cabecalho}\n" + ",".join(campos[c] for c in inst.COLUNAS_CONTRATO) + "\n",
        encoding="utf-8")
    evidencias = tmp_path / "evidencias_chembl.csv"
    evidencias.write_text("molecule_chembl_id,pchembl\nCHEMBL1,6.2\n", encoding="utf-8")
    saida = tmp_path / "input"

    codigo = inst.main(["--ligantes-origem", str(origem),
                        "--evidencias-origem", str(evidencias),
                        "--saida", str(saida)])

    assert codigo == 0
    assert (saida / "ligantes.tsv").exists()
    assert (saida / "evidencias_chembl.csv").exists()
    import json
    proc = json.loads((saida / "procedencia.json").read_text(encoding="utf-8"))
    papeis = {e["papel"] for e in proc["entradas"]}
    assert papeis == {"ligantes", "evidencias_chembl"}


def test_cli_devolve_3_quando_a_conversao_descarta_linha(tmp_path):
    """Codigo 3 e execucao PARCIAL: a saida existe, mas linhas foram
    descartadas. Automacao nao le stderr, entao o codigo tem de dizer."""
    origem = tmp_path / "compostos.csv"
    cabecalho = ",".join(inst.COLUNAS_CONTRATO)
    campos = {c: "" for c in inst.COLUNAS_CONTRATO}
    campos["composto_nome"] = "Cafeina"
    campos["smiles"] = CAFEINA
    campos["fonte_smiles"] = "NPASS"
    linha_ok = ",".join(campos[c] for c in inst.COLUNAS_CONTRATO)
    origem.write_text(f"{cabecalho}\n{linha_ok}\nTruncada,x\n", encoding="utf-8")

    codigo = inst.main(["--ligantes-origem", str(origem),
                        "--saida", str(tmp_path / "input")])
    assert codigo == 3


def test_nome_com_aspas_atravessa_intacto_ate_o_leitor_do_pipeline(tmp_path):
    """O defeito: o escritor usava QUOTE_MINIMAL e o leitor da etapa 1 usa
    QUOTE_NONE, por desenho ("TSV nao tem citacao: aspas fazem parte do
    nome"). As duas metades discordavam, e todo nome com `"` chegava torcido
    -- 13 registros no instantaneo de producao. Um deles e real:
    `5"-Methoxyhexahydrocurcumin`."""
    from mods.evidencias.entrada import ler_ligantes
    nome = '5"-Methoxyhexahydrocurcumin'
    origem = tmp_path / "compostos.csv"
    campos = {c: "" for c in inst.COLUNAS_CONTRATO}
    campos["composto_nome"] = nome
    campos["smiles"] = CAFEINA
    campos["fonte_smiles"] = "NPASS"
    import csv as _csv
    with open(origem, "w", encoding="utf-8", newline="") as f:
        w = _csv.DictWriter(f, fieldnames=inst.COLUNAS_CONTRATO)
        w.writeheader()
        w.writerow(campos)

    destino = tmp_path / "ligantes.tsv"
    inst.converter_ligantes(origem, destino)
    ligantes, problemas = ler_ligantes(str(destino))

    assert [p.mensagem for p in problemas] == []
    assert ligantes[0]["nome"] == nome


def test_tabulacao_no_nome_e_trocada_e_RELATADA_nao_quebra_a_linha(tmp_path):
    """Nome com tabulacao, escrito com citacao, seria partido em DUAS linhas
    pelo leitor QUOTE_NONE -- um composto viraria dois. A tabulacao e trocada
    por espaco e o caso vai para `conferencia`: marcado, nunca descartado em
    silencio, que e a regra do projeto."""
    origem = tmp_path / "compostos.csv"
    campos = {c: "" for c in inst.COLUNAS_CONTRATO}
    campos["composto_nome"] = "Com\ttabulacao"
    campos["smiles"] = CAFEINA
    campos["fonte_smiles"] = "NPASS"
    import csv as _csv
    with open(origem, "w", encoding="utf-8", newline="") as f:
        w = _csv.DictWriter(f, fieldnames=inst.COLUNAS_CONTRATO)
        w.writeheader()
        w.writerow(campos)

    destino = tmp_path / "ligantes.tsv"
    relatorio = inst.converter_ligantes(origem, destino)

    linhas = destino.read_text(encoding="utf-8").splitlines()
    assert len(linhas) == 2                      # cabecalho + 1 composto
    assert relatorio["sanitizadas"] == 1
    saida = _ler_saida(destino)
    assert saida[0]["nome"] == "Com tabulacao"
    assert "caractere de controle" in saida[0]["conferencia"].lower()
