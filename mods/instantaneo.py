#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
instantaneo.py

Converte a tabela de compostos do pipeline de prospecção de plantas
(`resultados/atividades_compostos_parcial.csv`) no formato de entrada deste
pipeline (TSV com cabeçalho), aplicando o MESMO tratamento químico que a
prospecção usa, para que as duas bases possam ser cruzadas depois.

Os dois projetos são repositórios SEPARADOS: este script não importa nada
do código da prospecção. Ele só conhece o nome das colunas que a prospecção
escreve, documentados abaixo, e as reimplementa de forma independente.

COMO A PROSPECÇÃO OBTÉM E TRATA OS DADOS (contrato replicado aqui)
-------------------------------------------------------------------
Nome       `composto_nome`, com `nome_pubchem` como alternativa. Nomes de
           produto natural costumam conter vírgula ("1,8-cineol",
           "3,3'-diindolilmetano") -- 42,8% deles na base real --, por isso
           a saída é TSV e nunca CSV.

SMILES     A prospecção coleta o SMILES em até quatro bases e escolhe o da
           primeira fonte cujo SMILES o RDKit consegue ler, nesta ordem de
           prioridade:

               NPASS > COCONUT > LOTUS > PubChem

           O escolhido fica em `smiles`, a fonte em `fonte_smiles`, e os
           candidatos de cada base continuam em `smiles_npass`,
           `smiles_coconut`, `smiles_lotus` e `smiles_pubchem`.

Conferência
           Quando duas ou mais bases têm SMILES para o mesmo composto, a
           prospecção compara as ESTRUTURAS canônicas (não o texto) e grava
           o veredito em `conferencia_smiles`:

               "mesma molécula"
               "estereoquímica diferente: NPASS/PubChem x COCONUT"
               "estruturas diferentes: ..."
               "SMILES inválido: <fonte>"

           Compostos cujo veredito aponta divergência estrutural entram na
           saída marcados, não descartados -- a decisão é do pesquisador.

InChIKey   `InChIKey` é a chave de identidade declarada pelas bases. Aqui
           recalculamos duas chaves a partir do SMILES:

             inchikey            da molécula como está (preserva sal e carga)
             inchikey_parental   depois de Cleanup -> maior fragmento
                                 orgânico -> neutralização

           A chave parental é a que reconhece que "cloridrato de X" e "X"
           são o mesmo composto, e é a coluna correta para cruzar as duas
           bases. Divergência entre o InChIKey declarado e o recalculado é
           reportada na coluna `conferencia`, sem descartar a linha.

USO
---
    python mods/instantaneo.py \\
        --ligantes-origem CAMINHO/atividades_compostos.csv \\
        --evidencias-origem CAMINHO/evidencias_chembl.csv \\
        --modelo-origem CAMINHO/models

Escreve `input/ligantes.tsv` (convertido), copia a evidencia medida para
`input/` e os pesos para `modelo/`, e grava `input/procedencia.json` com
origem, data, MD5 e contagens de cada entrada. `--limite N` para depois de N
compostos unicos (amostra de teste). `--sobrescrever` e NECESSARIO para
substituir entrada que ja existe: ela pode ter sido revisada a mao.

Codigos de saida: 0 completo, 3 parcial (linha truncada descartada), 2 erro.


    # amostra para teste rápido, priorizando compostos com evidência ChEMBL
    python mods/instantaneo.py --origem ... --saida input.exemplo.tsv \
        --limite 200 --preferir-com-chembl

REQUISITOS
----------
    pip install rdkit
"""

import argparse
import csv
import sys
from pathlib import Path

try:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
except ImportError:
    print("[ERRO] Este script precisa do RDKit: pip install rdkit", file=sys.stderr)
    sys.exit(1)

# A camada quimica mora num lugar so: e dela que sai a forma parental, que e a
# chave de identidade e de cruzamento entre as duas bases.
# Roda como script (`python mods/instantaneo.py`) alem de como modulo: as tres CLIs
# das etapas usam a mesma convencao, e o orquestrador as chama por
# caminho de arquivo.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from mods.quimica import forma_parental, inchikey_de, molecula, texto  # noqa: E402

# Ordem de prioridade das fontes de SMILES, igual à da prospecção.
ORDEM_FONTES_SMILES = [
    ("NPASS", "smiles_npass"),
    ("COCONUT", "smiles_coconut"),
    ("LOTUS", "smiles_lotus"),
    ("PubChem", "smiles_pubchem"),
]

# As colunas que a prospeccao TEM de ter. Lidas uma vez, no cabecalho: se
# uma delas mudar de nome, o importador RECUSA em vez de emitir campos vazios
# -- `linha.get(x)` devolvendo None e indistinguivel de celula vazia, e uma
# base importada com nomes vindos do PubChem porque `composto_nome` sumiu
# seria errada em silencio.
COLUNAS_CONTRATO = [
    "composto_nome", "nome_pubchem", "InChIKey", "smiles", "fonte_smiles",
    "smiles_npass", "smiles_coconut", "smiles_lotus", "smiles_pubchem",
    "conferencia_smiles", "chembl_id", "correspondencia_chembl",
    "pubchem_cid", "composto_ids", "plantas", "estados_composto", "n_plantas",
    "atividades_relatadas_composto", "melhor_evidencia_composto",
    "estudos_composto", "n_artigos_positivos", "bases", "atividades_duke",
    "referencias_lotus", "atividades_chembl",
]

# Identidade e quimica primeiro: quem abre o arquivo ve o composto antes da
# procedencia.
COLUNAS_IDENTIDADE = [
    "nome", "smiles", "inchikey", "inchikey_parental",
    "fonte_smiles", "conferencia", "chembl_id", "correspondencia_chembl",
]
# O que a prospeccao sabe e o importador descartava. Carregar custa colunas
# num TSV e devolve ao pipeline a procedencia botanica e a evidencia de
# atividade -- foi do `atividades_duke` que saiu o conjunto de validacao
# `dados/validacao_neuro_prospector.tsv`. O leitor da etapa 1 e guiado por
# cabecalho e repassa coluna desconhecida ao ligante, entao nada quebra.
COLUNAS_DA_PROSPECCAO = [
    "pubchem_cid", "composto_ids", "plantas", "estados_composto", "n_plantas",
    "bases", "atividades_duke", "atividades_chembl",
    "atividades_relatadas_composto", "melhor_evidencia_composto",
    "estudos_composto", "n_artigos_positivos", "referencias_lotus",
]
COLUNAS_SAIDA = COLUNAS_IDENTIDADE + COLUNAS_DA_PROSPECCAO


def colunas_ausentes(cabecalho) -> list[str]:
    """Colunas do contrato que o arquivo nao tem, na ordem do contrato.

    O BOM e tratado aqui tambem: a prospeccao escreve UTF-8 com BOM e
    `composto_nome` e a PRIMEIRA coluna, entao um leitor sem `utf-8-sig` ve
    `\ufeffcomposto_nome` e perderia o nome do composto em silencio."""
    presentes = {c.lstrip("\ufeff") for c in (cabecalho or [])}
    return [c for c in COLUNAS_CONTRATO if c not in presentes]

def escolher_smiles(linha: dict) -> tuple[str, str]:
    """(smiles, fonte), preferindo o que a PROSPECCAO escolheu.

    Ela ja escolheu comparando estruturas canonicas em quatro bases, e
    registrou a escolha em `smiles` + `fonte_smiles`. Medido no arquivo real
    (30.000 linhas): a ordem por fonte reimplementada aqui concordava com ela
    em 100% dos casos -- mesma estrutura e mesma fonte. Entao reimplementar a
    decisao nao ganhava nada e abria espaco para divergir em silencio se a
    regra dela mudasse.

    A ordem por fonte fica como FALLBACK, para o caso de o SMILES escolhido
    nao ser legivel; quando ela e usada, `converter_linha` registra isso na
    coluna `conferencia`."""
    escolhido = texto(linha.get("smiles"))
    if escolhido and molecula(escolhido) is not None:
        return escolhido, texto(linha.get("fonte_smiles")) or "prospeccao"
    for rotulo, coluna in ORDEM_FONTES_SMILES:
        smi = texto(linha.get(coluna))
        if smi and molecula(smi) is not None:
            return smi, rotulo
    return "", ""


def converter_linha(linha: dict) -> dict | None:
    """Uma linha da prospecção -> uma linha do formato de entrada, ou None
    se o composto não tiver nome ou SMILES utilizável."""
    nome = texto(linha.get("composto_nome")) or texto(linha.get("nome_pubchem"))
    smiles, fonte = escolher_smiles(linha)
    if not nome or not smiles:
        return None

    mol = molecula(smiles)
    pai = forma_parental(mol)
    ik = inchikey_de(mol)
    ik_pai = inchikey_de(pai)

    # O veredito da prospecção é preservado; acrescentamos a conferência da
    # identidade declarada contra a recalculada.
    notas = []
    veredito = texto(linha.get("conferencia_smiles"))
    if veredito:
        notas.append(veredito)
    declarado = texto(linha.get("InChIKey"))
    if declarado and ik and declarado != ik:
        if declarado.split("-")[0] == ik.split("-")[0]:
            notas.append(f"InChIKey declarado difere na estereoquímica/carga: {declarado}")
        else:
            notas.append(f"InChIKey declarado é outro esqueleto: {declarado}")

    # Fallback de fonte: a prospecção escolheu uma base cujo SMILES o RDKit
    # não lê aqui. Isso é divergência entre as duas, e fica escrito.
    escolhida = texto(linha.get("fonte_smiles"))
    if escolhida and fonte and fonte != escolhida:
        notas.append(f"fonte da prospeccao ({escolhida}) tem SMILES ilegivel; "
                     f"usada {fonte}")

    registro = {
        "nome": nome,
        # Grava o SMILES canônico: a prospecção compara estruturas, não texto.
        "smiles": Chem.MolToSmiles(mol),
        "inchikey": ik,
        "inchikey_parental": ik_pai,
        "fonte_smiles": fonte,
        "conferencia": "; ".join(notas),
        "chembl_id": texto(linha.get("chembl_id")),
        # Um `chembl_id` casado por nome nao vale o mesmo que um casado por
        # InChIKey, e o pipeline tem o conceito (`correspondencia_ligante` em
        # evidencias.tsv): jogar a distincao fora obriga a redescobri-la.
        "correspondencia_chembl": texto(linha.get("correspondencia_chembl")),
    }
    for coluna in COLUNAS_DA_PROSPECCAO:
        registro.setdefault(coluna, texto(linha.get(coluna)))
    return registro


# =========================================================================
# O instantaneo: procedencia, integridade e copia das entradas de fora
# =========================================================================
#
# `druglikeness2` nao le nada de outro projeto em tempo de execucao. O que vem
# de fora entra UMA VEZ, como instantaneo datado, com origem e MD5 gravados em
# `input/procedencia.json`.
#
# Isso nao e zelo: a base de origem e escrita por uma execucao EM CURSO -- o
# `_parcial` no nome dela e literal. Lida em tempo de execucao, a lista de
# ligantes mudaria debaixo do pipeline e duas etapas da mesma execucao
# poderiam ver conjuntos diferentes.

def _relativo(caminho) -> str:
    """Caminho relativo a raiz do projeto, quando esta dentro dela.

    `procedencia.json` e versionado: um caminho absoluto ali vazaria a pasta
    pessoal de quem rodou e nao resolveria na maquina de mais ninguem. A
    ORIGEM continua absoluta -- e justamente o que se quer registrar.
    """
    from mods import caminhos
    caminho = Path(caminho).resolve()
    try:
        return str(caminho.relative_to(caminhos.RAIZ))
    except ValueError:
        return str(caminho)


def md5_de(caminho) -> str:
    """MD5 do arquivo, lido em blocos (a base de origem passa de 100 MB)."""
    import hashlib
    h = hashlib.md5()
    with open(caminho, "rb") as f:
        for bloco in iter(lambda: f.read(1 << 20), b""):
            h.update(bloco)
    return h.hexdigest()


_CONTROLE = {"\t": "tabulacao", "\n": "nova linha", "\r": "retorno de carro"}


def _sem_controle(linha: dict) -> list[str]:
    """Troca tabulacao e quebra de linha por espaco, e devolve o que trocou.

    Nao e cosmetica: o leitor da etapa 1 usa `QUOTE_NONE`, porque em TSV as
    aspas fazem parte do nome. Um valor com tabulacao, escrito com citacao,
    seria PARTIDO EM DUAS LINHAS por aquele leitor -- um composto viraria
    dois. Trocar e marcar preserva o composto e conta a verdade; a regra do
    projeto e marcar, nunca descartar em silencio.
    """
    achados = []
    for campo, valor in list(linha.items()):
        if not isinstance(valor, str):
            continue
        presentes = [nome for ch, nome in _CONTROLE.items() if ch in valor]
        if not presentes:
            continue
        for ch in _CONTROLE:
            valor = valor.replace(ch, " ")
        linha[campo] = " ".join(valor.split())
        achados.append(f"{campo}: {', '.join(presentes)}")
    return achados


def converter_ligantes(origem, destino, *, limite: int | None = None) -> dict:
    """Converte a tabela de compostos da base de prospeccao no formato de
    entrada e devolve o relatorio da conversao.

    Diferenca em relacao a uma leitura ingenua: **linha truncada e contada e
    nomeada**, nunca somada ao bolo de "sem SMILES utilizavel". Um instantaneo
    tirado de arquivo em escrita pode terminar no meio de uma linha, e perder
    composto em silencio e o pior resultado possivel aqui.
    """
    origem, destino = Path(origem), Path(destino)
    relatorio = {"lidas": 0, "truncadas": 0, "invalidas": 0,
                 "duplicadas": 0, "sanitizadas": 0, "gravados": 0,
                 "problemas": []}
    vistos: dict[str, dict] = {}
    with open(origem, "r", encoding="utf-8-sig", newline="") as f:
        leitor = csv.DictReader(f)
        ausentes = colunas_ausentes(leitor.fieldnames)
        if ausentes:
            raise ValueError(
                f"{origem} nao tem as colunas do contrato: "
                f"{', '.join(ausentes)}")
        n_campos = len(leitor.fieldnames)
        for linha in leitor:
            relatorio["lidas"] += 1
            # DictReader preenche com None o que falta: linha cortada no meio.
            if None in linha.values() or len(linha) < n_campos:
                relatorio["truncadas"] += 1
                relatorio["problemas"].append(
                    f"linha {relatorio['lidas']}: truncada "
                    f"({sum(v is not None for v in linha.values())} de "
                    f"{n_campos} campos)")
                continue
            convertida = converter_linha(linha)
            if convertida is None:
                relatorio["invalidas"] += 1
                continue
            trocados = _sem_controle(convertida)
            if trocados:
                relatorio["sanitizadas"] += 1
                nota = ("caractere de controle trocado por espaco ("
                        + "; ".join(trocados) + ")")
                convertida["conferencia"] = (
                    f"{convertida['conferencia']}; {nota}"
                    if convertida.get("conferencia") else nota)
            chave = convertida["inchikey_parental"] or convertida["inchikey"]
            if chave in vistos:
                relatorio["duplicadas"] += 1
                continue
            vistos[chave] = convertida
            if limite is not None and len(vistos) >= limite:
                break

    destino.parent.mkdir(parents=True, exist_ok=True)
    with open(destino, "w", encoding="utf-8", newline="") as f:
        f.write("#" + "\t".join(COLUNAS_SAIDA) + "\n")
        # QUOTE_NONE: em TSV as aspas FAZEM PARTE do nome, e e assim que o
        # leitor da etapa 1 le (`mods/evidencias/entrada.py`). Escrever com
        # citacao faria todo nome com `"` chegar torcido do outro lado.
        escritor = csv.DictWriter(f, fieldnames=COLUNAS_SAIDA, delimiter="\t",
                                  lineterminator="\n", extrasaction="ignore",
                                  quoting=csv.QUOTE_NONE, quotechar=None,
                                  escapechar=None)
        for linha in vistos.values():
            escritor.writerow(linha)
    relatorio["gravados"] = len(vistos)
    return relatorio


def copiar_entrada(origem, pasta, *, papel: str,
                   sobrescrever: bool = False) -> dict:
    """Copia um arquivo de fora para `pasta` e devolve sua procedencia.

    Recusa sobrescrever sem permissao explicita: a entrada da execucao
    anterior pode ser uma lista que o pesquisador revisou a mao, e apaga-la em
    silencio e exatamente o erro que o projeto nao comete.
    """
    import shutil
    origem, pasta = Path(origem), Path(pasta)
    if not origem.is_file():
        raise FileNotFoundError(f"entrada nao encontrada: {origem}")
    alvo = pasta / origem.name
    existia = alvo.exists()
    if existia and not sobrescrever:
        raise FileExistsError(
            f"{alvo} ja existe. Passe sobrescrever=True para substituir -- a "
            f"entrada atual pode ter sido revisada a mao.")
    pasta.mkdir(parents=True, exist_ok=True)
    shutil.copy2(origem, alvo)
    return {"papel": papel, "origem": str(origem),
            "destino": _relativo(alvo),
            "md5": md5_de(origem), "bytes": origem.stat().st_size,
            "substituiu": existia}


def registrar(entradas: list[dict], destino) -> Path:
    """Escreve `input/procedencia.json`: de qual estado de cada base de origem
    este instantaneo saiu. E o que liga o manifesto de uma etapa a entrada."""
    import json
    from datetime import datetime, timezone
    destino = Path(destino)
    registro = {"data": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "entradas": []}
    for e in entradas:
        origem = Path(e["origem"])
        item = {k: v for k, v in e.items() if k != "origem"}
        item["origem"] = str(origem)
        item.setdefault("md5", md5_de(origem))
        item.setdefault("bytes", origem.stat().st_size)
        registro["entradas"].append(item)
    destino.parent.mkdir(parents=True, exist_ok=True)
    destino.write_text(json.dumps(registro, ensure_ascii=False, indent=2) + "\n",
                       encoding="utf-8")
    return destino


def _argumentos(argv):
    parser = argparse.ArgumentParser(
        description="Tira o instantaneo das entradas de fora: converte a "
                    "tabela de compostos, copia o que for preciso e grava a "
                    "procedencia (origem, data e MD5) de cada arquivo.")
    parser.add_argument("--ligantes-origem", dest="ligantes_origem",
                        help="CSV de compostos da base de prospeccao.")
    parser.add_argument("--evidencias-origem", dest="evidencias_origem",
                        help="CSV de evidencia medida (ChEMBL).")
    parser.add_argument("--modelo-origem", dest="modelo_origem",
                        help="Pasta com os pesos .pt e o target_labels.tsv.")
    parser.add_argument("--saida", default=None,
                        help="Pasta de entrada do projeto (padrao: input/).")
    parser.add_argument("--modelo-saida", dest="modelo_saida", default=None,
                        help="Pasta do modelo (padrao: modelo/).")
    parser.add_argument("--fingerprint", default="ecfp4",
                        help="Prefixo do checkpoint a copiar (padrao: ecfp4).")
    parser.add_argument("--limite", type=int, default=None,
                        help="Para depois de N compostos unicos. Serve para "
                             "montar amostra de teste, nao para producao.")
    parser.add_argument("--sobrescrever", action="store_true",
                        help="Substitui entradas que ja existem. Sem isto, "
                             "arquivo existente e ERRO -- ele pode ter sido "
                             "revisado a mao.")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    """Codigos de saida: 0 completo, 3 parcial (linha descartada), 2 erro."""
    from mods import caminhos
    args = _argumentos(argv)
    if not any((args.ligantes_origem, args.evidencias_origem,
                args.modelo_origem)):
        print("[ERRO] Informe ao menos uma origem.", file=sys.stderr)
        return 2
    pasta = Path(args.saida) if args.saida else caminhos.INPUT
    pasta_modelo = (Path(args.modelo_saida) if args.modelo_saida
                    else caminhos.MODELO)
    entradas, parcial = [], False
    try:
        if args.ligantes_origem:
            origem = Path(args.ligantes_origem)
            destino = pasta / "ligantes.tsv"
            if destino.exists() and not args.sobrescrever:
                print(f"[ERRO] {destino} ja existe; use --sobrescrever.",
                      file=sys.stderr)
                return 2
            rel = converter_ligantes(origem, destino, limite=args.limite)
            for problema in rel["problemas"][:20]:
                print(f"[AVISO] {problema}", file=sys.stderr)
            print(f"[INFO] Linhas lidas               : {rel['lidas']:,}")
            print(f"[INFO]   truncadas (descartadas)  : {rel['truncadas']:,}")
            print(f"[INFO]   sem nome ou SMILES       : {rel['invalidas']:,}")
            print(f"[INFO]   duplicadas (parental)    : {rel['duplicadas']:,}")
            print(f"[INFO] Gravados em {destino}      : {rel['gravados']:,}")
            # Codigo 3 sai por descarte ANORMAL (linha truncada), nao pelo
            # esperado. Produto natural sem SMILES utilizavel e a regra, nao a
            # excecao -- foram 17.305 de 126.240 na base real, 13,7%. Se isso
            # virasse 3, todo instantaneo sairia 3 e o codigo deixaria de
            # informar o que existe para informar. As contagens de `invalidas`
            # e `duplicadas` vao no relatorio e em procedencia.json.
            if rel["truncadas"]:
                parcial = True
            entradas.append({"papel": "ligantes", "origem": origem,
                             "destino": _relativo(destino),
                             "linhas": rel["lidas"],
                             "gravados": rel["gravados"],
                             "truncadas": rel["truncadas"],
                             "invalidas": rel["invalidas"],
                             "duplicadas": rel["duplicadas"]})
        if args.evidencias_origem:
            entradas.append(copiar_entrada(
                Path(args.evidencias_origem), pasta,
                papel="evidencias_chembl", sobrescrever=args.sobrescrever))
        if args.modelo_origem:
            entradas.extend(_copiar_modelo(
                Path(args.modelo_origem), pasta_modelo,
                prefixo=args.fingerprint, sobrescrever=args.sobrescrever))
    except (FileNotFoundError, FileExistsError, ValueError) as e:
        print(f"[ERRO] {e}", file=sys.stderr)
        return 2

    registro = registrar(entradas, pasta / "procedencia.json")
    print(f"[INFO] Procedencia em {registro}")
    if parcial:
        print("[AVISO] Execucao PARCIAL: linhas descartadas (codigo 3).",
              file=sys.stderr)
    return 3 if parcial else 0


def _copiar_modelo(origem: Path, destino: Path, *, prefixo: str,
                   sobrescrever: bool) -> list[dict]:
    """Leva o checkpoint e o `target_labels.tsv` para a pasta do projeto.

    O par tem de ser do MESMO PP4CKD: a coluna *i* de um nao e o alvo *i* do
    outro, e pontuar assim atribuiria a probabilidade de um alvo ao nome de
    outro. Por isso os dois sao copiados juntos, numa chamada.
    """
    rotulos = origem / "target_labels.tsv"
    if not rotulos.is_file():
        # o tarball publicado separa os pesos dos rotulos; aceita um nivel acima
        acima = origem.parent.parent / "target_labels.tsv"
        rotulos = acima if acima.is_file() else rotulos
    checkpoints = sorted(origem.glob(f"{prefixo}_full_model.pt")) or \
        sorted(origem.glob(f"{prefixo}_fold_*.pt"))
    if not checkpoints:
        raise FileNotFoundError(
            f"Nenhum checkpoint '{prefixo}_*' em {origem}.")
    if not rotulos.is_file():
        raise FileNotFoundError(
            f"target_labels.tsv nao encontrado em {origem} nem um nivel "
            f"acima. Ele e OBRIGATORIO: diz quais alvos existem no modelo.")
    return [copiar_entrada(checkpoints[0], destino, papel="checkpoint",
                           sobrescrever=sobrescrever),
            copiar_entrada(rotulos, destino, papel="rotulos",
                           sobrescrever=sobrescrever)]


if __name__ == "__main__":
    sys.exit(main())
