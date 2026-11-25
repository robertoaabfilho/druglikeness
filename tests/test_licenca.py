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

def test_requirements_lista_as_dependencias():
    linhas = [
        linha.split("#", 1)[0].strip()
        for linha in (RAIZ / "requirements.txt").read_text(encoding="utf-8").splitlines()
    ]
    linhas = [linha for linha in linhas if linha]
    nomes = {re.split(r"[<>=!~\s\[]", linha, maxsplit=1)[0].lower() for linha in linhas}
    for pacote in ("rdkit", "torch", "requests", "pandas", "streamlit", "pytest"):
        assert pacote in nomes
    assert "numpy>=2" in linhas
