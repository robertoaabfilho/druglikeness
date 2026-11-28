# -*- coding: utf-8 -*-
"""
quimica.py - identidade quimica: molecula, forma parental, InChIKey.

A **forma parental** e a chave de identidade deste projeto. Ela e obtida por
*Cleanup -> maior fragmento organico -> neutralizacao*, e e a que reconhece que
sal e base livre sao o mesmo composto. Duas consequencias que valem dizer:

- A **deduplicacao de ligantes** e feita por ela, nao pelo InChIKey declarado.
- E o **cruzamento com a base de prospeccao** depende dela: a prospecao aplica
  o mesmo tratamento, e por isso as duas bases se encontram.

As duas colunas existem por um motivo: `inchikey` e o que a molecula declara,
`inchikey_parental` e a identidade. Exemplo real: a neferina tem dois `[NH+]`
no SMILES, entao seu `inchikey` sai `...-P`, enquanto o `...-N` que as bases
declaram so aparece na forma parental.

**Nada aqui levanta excecao por molecula ruim.** SMILES ilegivel vira `None` e
a linha e JULGADA depois (o portao a reprova, com o motivo escrito), nunca
descartada em silencio na leitura.
"""
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

# Os padronizadores do RDKit sao caros de construir e nao tem estado entre
# chamadas: um por processo.
_PADRONIZADORES: dict = {}

# Como o CSV de origem escreve "nao tem valor". Todas sao a mesma coisa, e
# nenhuma pode virar o texto literal "nan" dentro de um rotulo de identidade.
_AUSENTES = ("", "nan", "none", "-", "null")


def _padronizadores() -> dict:
    if not _PADRONIZADORES:
        _PADRONIZADORES["frag"] = rdMolStandardize.LargestFragmentChooser(
            preferOrganic=True)
        _PADRONIZADORES["uncharger"] = rdMolStandardize.Uncharger()
    return _PADRONIZADORES


def texto(valor) -> str:
    """Normaliza um campo de tabela para string limpa ('' quando ausente)."""
    if valor is None:
        return ""
    s = str(valor).strip()
    return "" if s.lower() in _AUSENTES else s


def molecula(smiles: str):
    """Molecula do RDKit, ou None se o SMILES nao for legivel."""
    if not smiles:
        return None
    try:
        return Chem.MolFromSmiles(smiles)
    except Exception:
        return None


def forma_parental(mol):
    """Cleanup -> maior fragmento organico -> neutralizacao. O mesmo tratamento
    que a base de prospeccao aplica, para que as duas se cruzem."""
    if mol is None:
        return None
    try:
        p = _padronizadores()
        m = rdMolStandardize.Cleanup(mol)
        m = p["frag"].choose(m)
        return p["uncharger"].uncharge(m)
    except Exception:
        return None


def inchikey_de(mol) -> str:
    """InChIKey da molecula, ou '' quando nao ha molecula."""
    if mol is None:
        return ""
    try:
        return Chem.MolToInchiKey(mol) or ""
    except Exception:
        return ""
