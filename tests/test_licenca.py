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

def test_opcionais_tem_arquivo_proprio_com_a_variante_binaria():
    _, nomes = _requisitos("requirements-opcionais.txt")
    assert {"torch", "transformers", "psycopg2-binary"} <= nomes
    assert "psycopg2" not in nomes

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
