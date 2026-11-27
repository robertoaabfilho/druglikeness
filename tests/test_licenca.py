import re
from pathlib import Path
RAIZ = Path(__file__).resolve().parent.parent

def test_licenca_mit_na_raiz():
    texto = (RAIZ / "LICENSE").read_text(encoding="utf-8")
    assert "MIT License" in texto
    assert "Copyright (c) 2026 Roberto Alexandre Filho" in texto
    assert "Permission is hereby granted" in texto
    assert "WITHOUT WARRANTY" in texto

def test_dados_derivados_do_chembl_sao_share_alike():
    texto = " ".join((RAIZ / "dados" / "LICENSE").read_text(encoding="utf-8").split())
    assert "CC BY-SA 3.0" in texto
    assert "ChEMBL" in texto
    assert "10.1093/nar/gkad1004" in texto
    assert "se aplica apenas a estes dados" in texto
    assert "distribuido sob a licenca MIT" in texto

def _requisitos(nome):
    linhas = [
        linha.split("#", 1)[0].strip()
        for linha in (RAIZ / nome).read_text(encoding="utf-8").splitlines()
    ]
    linhas = [linha for linha in linhas if linha]
    nomes = {re.split(r"[<>=!~\s\[]", linha, maxsplit=1)[0].lower() for linha in linhas}
    return linhas, nomes

def test_requirements_lista_as_dependencias():
    linhas, nomes = _requisitos("requirements.txt")
    for pacote in ("rdkit", "requests", "pandas", "streamlit", "pytest", "scipy"):
        assert pacote in nomes
    assert "numpy>=2" in linhas

def test_requirements_e_instalavel_como_esta():
    """`pip install -r requirements.txt` e o primeiro comando do README. O pip
    instala toda linha, comentada como opcional ou nao: psycopg2 (sem -binary)
    exige pg_config e quebrava o arquivo inteiro num clone limpo, e as duas
    dependencias de alguns GB nao podem estar no caminho padrao."""
    _, nomes = _requisitos("requirements.txt")
    assert "psycopg2" not in nomes
    assert not ({"torch", "transformers"} & nomes)

def test_opcionais_lista_so_o_que_a_arvore_realmente_usa():
    """`torch` liga a evidencia PREVISTA e tem consumidor. `transformers` e
    `psycopg2-binary` nao tem nenhum aqui -- lista-los faria alguem baixar
    alguns GB para nada. E `psycopg2` puro nunca entra: ele exige pg_config e
    quebra num clone limpo; se um dia voltar, volta como -binary."""
    _, nomes = _requisitos("requirements-opcionais.txt")
    assert "torch" in nomes
    assert "psycopg2" not in nomes
    assert "psycopg2-binary" not in nomes
    assert "transformers" not in nomes


def test_nenhuma_opcional_listada_fica_sem_consumidor():
    """A guarda geral: cada pacote opcional tem de aparecer em algum .py."""
    import pathlib as _p
    _, nomes = _requisitos("requirements-opcionais.txt")
    fontes = " ".join(
        c.read_text(encoding="utf-8")
        for c in (RAIZ / "mods").rglob("*.py"))
    for pacote in nomes:
        modulo = pacote.replace("-binary", "").replace("-", "_")
        assert modulo in fontes, f"{pacote} nao tem consumidor em mods/"

def test_readme_ensina_a_instalar_as_opcionais():
    texto = (RAIZ / "README.md").read_text(encoding="utf-8")
    assert "requirements-opcionais.txt" in texto

def test_fixtures_derivadas_do_chembl_tambem_sao_share_alike():
    """Dado do ChEMBL versionado FORA de dados/ tambem carrega o share-alike;
    sob a MIT da raiz, estes arquivos estariam relicenciados."""
    texto = " ".join((RAIZ / "tests" / "dados" / "LICENSE").read_text(
        encoding="utf-8").split())
    assert "CC BY-SA 3.0" in texto
    assert "10.1093/nar/gkad1004" in texto
    for caminho in ("tests/dados/evidencias_amostra.csv",
                    "tests/dados/respostas/molecule.json_*.json",
                    "tests/dados/respostas/target_component.json_*.json"):
        assert caminho in texto
    # Reactome e Open Targets sao CC0 e ficam nomeados como tal
    assert "CC0" in texto and "reactome_*.json" in texto
    assert "Functional Interaction" in texto

def test_todo_dado_do_chembl_versionado_tem_licenca_ao_lado():
    """Qualquer pasta nova com derivado do ChEMBL precisa de LICENSE propria."""
    for pasta in (RAIZ / "dados", RAIZ / "tests" / "dados"):
        assert (pasta / "LICENSE").is_file(), pasta


def test_dado_do_prospector_em_dados_nao_e_atribuido_ao_chembl():
    """`dados/LICENSE` dizia que TODO arquivo daquela pasta deriva do ChEMBL.
    O conjunto de validacao nao deriva: ele vem do Duke (USDA, dominio
    publico) e de SMILES de NPASS/COCONUT/LOTUS/PubChem, via prospeccao.
    Atribuir a fonte errada e erro de licenciamento, nao de estilo."""
    texto = " ".join((RAIZ / "dados" / "LICENSE").read_text(
        encoding="utf-8").split())
    assert "validacao_neuro_prospector.tsv" in texto
    assert "Duke" in texto
    assert "NPASS" in texto and "COCONUT" in texto and "LOTUS" in texto
    # e o escopo do share-alike do ChEMBL fica nomeado, nao implicito
    assert "evidencias_controle.tsv" in texto


def test_conjunto_de_validacao_tem_leia_me_ao_lado():
    leia_me = RAIZ / "dados" / "validacao_neuro_prospector.LEIA-ME.md"
    assert leia_me.is_file()
    texto = leia_me.read_text(encoding="utf-8")
    # o que torna o conjunto reproduzivel sem o script: a regra de selecao
    assert "atividades_duke" in texto
    assert "Neurotoxic" in texto and "Analgesic" in texto
    assert "eac618972d4900e0bd684597b49f654c" in texto      # MD5 da origem
    assert "Capsaicin" in texto                              # a contaminacao
