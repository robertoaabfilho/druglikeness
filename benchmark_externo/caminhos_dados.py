"""De onde vem o dado grande que estes scripts leem.

Nenhum caminho absoluto fica escrito aqui. Cada um vem de variavel de
ambiente, e a falta de uma vira erro com o nome dela -- em vez de um
caminho que existia na maquina de quem rodou primeiro.

    LITPCBA_FULL        pasta LIT-PCBA_full
    LITPCBA_AVE         pasta LIT-PCBA_AVE_unbiased
    CHECKPOINT_PP4CKD   arquivo .pt do checkpoint
    ROTULOS_PP4CKD      target_labels.tsv

A raiz do repositorio e deduzida do proprio arquivo, sem variavel.
"""
import os
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

_DESCRICAO = {
    "LITPCBA_FULL": "pasta do LIT-PCBA completo (LIT-PCBA_full)",
    "LITPCBA_AVE": "pasta do LIT-PCBA AVE-unbiased",
    "CHECKPOINT_PP4CKD": "arquivo .pt do checkpoint do PP4CKD",
    "ROTULOS_PP4CKD": "target_labels.tsv do PP4CKD",
}


def de_ambiente(nome: str) -> Path:
    v = os.environ.get(nome)
    if not v:
        raise SystemExit(
            f"falta a variavel de ambiente {nome} "
            f"({_DESCRICAO.get(nome, 'caminho do dado')}).\n"
            f"Defina-a antes de rodar, por exemplo:\n"
            f"  {nome}=/caminho/para/o/dado python3 -I {Path.cwd().name}/...")
    return Path(v)
