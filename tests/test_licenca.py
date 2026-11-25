from pathlib import Path
RAIZ = Path(__file__).resolve().parent.parent

def test_licenca_mit_na_raiz():
    texto = (RAIZ / "LICENSE").read_text(encoding="utf-8")
    assert "MIT License" in texto

def test_dados_derivados_do_chembl_sao_share_alike():
    texto = (RAIZ / "dados" / "LICENSE").read_text(encoding="utf-8")
    assert "CC BY-SA 3.0" in texto
    assert "ChEMBL" in texto

def test_requirements_lista_as_dependencias():
    texto = (RAIZ / "requirements.txt").read_text(encoding="utf-8").lower()
    for pacote in ("rdkit", "torch", "requests", "pandas", "streamlit", "pytest"):
        assert pacote in texto
