# -*- coding: utf-8 -*-
"""
dominio.py

Quando a PREVISAO do modelo tem habilidade medida, e quando ela e azar.

O que foi medido
----------------
`benchmark_externo/`, adendo 8: holdout temporal em 145 alvos do painel
ChEMBL 28, 1.769 ativos publicados a partir de 2022 contra 9.780 decoys
bioativos (`actives-as-decoys`). O eixo que separa acerto de azar NAO e a
probabilidade do modelo -- e a distancia estrutural do composto aos ativos
JA CONHECIDOS do alvo, medida como `sim_top5`: a media das 5 maiores
similaridades Tanimoto ECFP4 aos ativos de referencia daquele alvo.

    sim_top5 <  0,20   AUROC 0,4987   azar exato
    sim_top5 >= 0,20   AUROC 0,6819

Por decil, a transicao e abrupta (`resultado_decis_similaridade.tsv`):
quatro dos cinco primeiros decis ficam indistinguiveis de 0,5 e dois ficam
ABAIXO; do sexto em diante, 0,65 a 0,71.

Comparacao com o que o projeto usava como proxy de dominio, no mesmo dado
(`resultado_proxy_dominio.tsv`, criterio: parte o estrato na mediana do proxy
e compara a AUROC do modelo nas duas metades):

| proxy                            | alta  | baixa | delta  |    z |
|----------------------------------|-------|-------|--------|------|
| **media das 5 maiores (sim_top5)** | 0,682 | 0,499 | +0,183 | +12,1 |
| similaridade maxima              | 0,676 | 0,533 | +0,144 |  +9,6 |
| numero de bits do ECFP4          | 0,640 | 0,514 | +0,127 |  +7,1 |
| probabilidade maxima entre cabecas | 0,654 | 0,590 | +0,065 |  +4,2 |

Por que ha um minimo de referencias
-----------------------------------
`sim_top5` e a media das 5 MAIORES similaridades. Com poucas referencias, as
5 maiores sao simplesmente as unicas que existem, e o valor cai por
construcao -- nao porque o composto seja distante do conhecido, mas porque
pouco e conhecido. Medido nos 9 alvos de dor do projeto contra os 107.123
ligantes (`resultado_posicao_biblioteca.tsv`):

| alvo       | referencias | fracao com sim_top5 < 0,20 |
|------------|-------------|-----------------------------|
| CHEMBL233  | 5.790       | 28,5%  |
| CHEMBL228  | 4.434       | 53,5%  |
| CHEMBL4794 | 3.085       | 45,4%  |
| CHEMBL222  | 2.856       | 61,5%  |
| CHEMBL1904 | 1.549       | 75,7%  |
| CHEMBL1919 |   373       | 96,4%  |
| CHEMBL2015 |   113       | 96,5%  |
| CHEMBL1845 |    67       | 98,8%  |
| CHEMBL3896 |    19       | 99,9%  |

O salto esta entre 1.549 e 373 referencias. Abaixo disso, `sim_top5 < 0,20`
nao quer dizer "fora do dominio", quer dizer "nao da para saber" -- e punir
o composto por isso rebaixaria o alvo inteiro por falta de bibliografia, nao
por falta de sinal. Entao esses alvos recebem `indeterminado`.

O valor de 500 e DECLARADO, nao calibrado: o dado so o constrange ao
intervalo aberto (373, 1.549). Escolhido o numero redondo mais baixo que
separa, para nao descartar alvo por folga.

Por que ha faixa de TAMANHO, e por que ela tem dois lados
---------------------------------------------------------
O adendo 7 mediu que o numero de bits ligados do ECFP4 escala com o tamanho
da molecula, e que 440 das 616 cabecas do modelo tem rho > +0,3 com essa
contagem. O adendo 8 mediu que DESCONTAR esse efeito piora a AUROC (0,6322
para 0,5637): tamanho esta confundido com atividade real, porque ligante
maior faz mais contatos. O defeito nao e o tamanho em si -- e a
EXTRAPOLACAO, molecula fora da faixa de tamanho onde a relacao foi aprendida.

Isso contamina tambem o eixo de similaridade, dos dois lados, e os dois
casos apareceram na verificacao desta implementacao:

* molecula PEQUENA tem poucos bits, e poucos bits sobrepoem facil. A cafeina
  (14 atomos pesados) saiu `analogo` no NMDA 2B, cujas referencias vao de 20
  a 31 atomos no intervalo p1-p99.
* molecula GRANDE tem bits para sobrepor com qualquer coisa. E o caso da
  vancomicina marcando 0,92 em BTK no adendo 7.

Entao a faixa e dois lados, por alvo, a partir do tamanho das PROPRIAS
referencias dele. Fora dela a zona e `extrapolado`, que ordena junto de
`fora`: a similaridade existe, mas nao mede o que se quer que meca.

Limites de validade
-------------------
* A calibracao de 0,20 vem de um painel de alvos MADUROS (>= 50 ativos
  publicados depois de 2022). Nao foi medida em alvo pouco explorado.
* Vem do modelo multitarefa do ChEMBL 28, nao do PP4CKD. Para o PP4CKD nao
  existe numero sem vazamento com dado do ChEMBL (adendo 5).
* `0,20` e `0,25` entram no manifesto da execucao. Nao sao constante de
  codigo a ser esquecida: sao parametro do ranqueamento.
"""

# O fingerprint do EIXO DE DOMINIO e fixo, e NAO o do modelo.
#
# 0,20 foi calibrado com ECFP4 raio 2 em 1024 bits (o modelo multitarefa do
# ChEMBL 28). O PP4CKD usa 4096 bits, e dobragem diferente da Tanimoto
# diferente: com menos bits ha mais colisao e a similaridade sai inflada.
# Aplicar o limite num espaco e calibra-lo noutro compararia numeros que nao
# se comparam. Este eixo mede distancia estrutural entre composto e ativo
# conhecido -- nao e a entrada da rede, e nao tem de acompanhar a largura
# dela.
RAIO_DOMINIO = 2
BITS_DOMINIO = 1024

LIMITE_BANDA = 0.20
LIMITE_ANALOGO = 0.25
MIN_REFERENCIAS = 500

ANALOGO = "analogo"
BANDA = "banda"
FORA = "fora"
EXTRAPOLADO = "extrapolado"
INDETERMINADO = "indeterminado"

# Ordem de confiabilidade para o desempate do ranking. Menor ordena antes.
#
# `indeterminado` fica ENTRE os dois medidos, e nao junto de `fora`, porque
# nao saber nao e a mesma coisa que saber que e azar. Rebaixar o
# indeterminado ao fundo puniria o alvo por falta de bibliografia.
#
# `extrapolado` fica junto de `fora`: ali a similaridade NAO e falta de
# informacao, e informacao enganosa -- ver abaixo.
_ORDEM = {ANALOGO: 0, BANDA: 0, INDETERMINADO: 1, FORA: 2, EXTRAPOLADO: 2}

# Percentis do tamanho das referencias que delimitam a faixa de validade.
# p1-p99 em vez de min-max para que UMA referencia atipica nao abra o portao.
PERCENTIL_BAIXO, PERCENTIL_ALTO = 1, 99


def sim_top5(fingerprint, fingerprints_de_referencia, k: int = 5):
    """Media das `k` maiores similaridades Tanimoto, ou None sem referencia."""
    if fingerprint is None or not fingerprints_de_referencia:
        return None
    from rdkit import DataStructs
    sims = DataStructs.BulkTanimotoSimilarity(
        fingerprint, list(fingerprints_de_referencia))
    sims.sort(reverse=True)
    maiores = sims[:k]
    return sum(maiores) / len(maiores)


def zona(similaridade, n_referencias: int, tamanho=None,
         faixa_de_tamanho=None) -> str:
    """Em que zona de confiabilidade a previsao cai.

    `indeterminado` quando faltam referencias para o valor significar algo,
    ou quando nao houve similaridade a calcular. `extrapolado` quando o
    composto esta fora da faixa de tamanho das referencias do alvo.
    """
    if n_referencias < MIN_REFERENCIAS or similaridade is None:
        return INDETERMINADO
    if tamanho is not None and faixa_de_tamanho:
        baixo, alto = faixa_de_tamanho
        if tamanho < baixo or tamanho > alto:
            return EXTRAPOLADO
    if similaridade >= LIMITE_ANALOGO:
        return ANALOGO
    if similaridade >= LIMITE_BANDA:
        return BANDA
    return FORA


def ordem(zona_do_composto: str | None) -> int:
    """Posicao de desempate. Zona ausente ordena como `indeterminado`."""
    if not zona_do_composto:
        return _ORDEM[INDETERMINADO]
    return _ORDEM.get(zona_do_composto, _ORDEM[INDETERMINADO])


def parametros() -> dict:
    """O que vai para o manifesto da execucao."""
    return {"limite_banda": LIMITE_BANDA,
            "limite_analogo": LIMITE_ANALOGO,
            "minimo_referencias": MIN_REFERENCIAS,
            "k_vizinhos": 5,
            "raio": RAIO_DOMINIO,
            "bits": BITS_DOMINIO,
            "percentis_de_tamanho": [PERCENTIL_BAIXO, PERCENTIL_ALTO]}


def carregar_referencias(caminho) -> dict[str, list[str]]:
    """{alvo_chembl_id: [smiles]} do TSV de referencias.

    Erro explicito quando o arquivo nao existe ou nao tem as colunas: a
    alternativa seria devolver vazio, toda zona sairia `indeterminado` e o
    ranqueamento pareceria ter medido o dominio sem ter medido.
    """
    import csv
    from pathlib import Path

    p = Path(caminho)
    if not p.is_file():
        raise OSError(f"referencias de alvo nao encontradas: {p}")
    por_alvo: dict[str, list[str]] = {}
    with open(p, encoding="utf-8", newline="") as fh:
        leitor = csv.reader(fh, delimiter="\t")
        cabecalho = next(leitor, None)
        if not cabecalho:
            raise ValueError(f"{p}: arquivo vazio")
        nomes = [c.lstrip("#").strip() for c in cabecalho]
        try:
            i_alvo = nomes.index("alvo_chembl_id")
            i_smiles = nomes.index("smiles")
        except ValueError:
            raise ValueError(
                f"{p}: precisa das colunas alvo_chembl_id e smiles; "
                f"achei {nomes}") from None
        for linha in leitor:
            if len(linha) <= max(i_alvo, i_smiles):
                continue
            alvo, smiles = linha[i_alvo].strip(), linha[i_smiles].strip()
            if alvo and smiles:
                por_alvo.setdefault(alvo, []).append(smiles)
    if not por_alvo:
        raise ValueError(f"{p}: nenhuma referencia utilizavel")
    return por_alvo


def gerador():
    """O gerador de fingerprint do eixo de dominio, na largura calibrada."""
    from rdkit.Chem import rdFingerprintGenerator
    return rdFingerprintGenerator.GetMorganGenerator(
        radius=RAIO_DOMINIO, fpSize=BITS_DOMINIO)


def fingerprint(smiles: str, gen=None):
    """Fingerprint do eixo de dominio, ou None para SMILES ilegivel."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles or "")
    if mol is None:
        return None
    return (gen or gerador()).GetFingerprint(mol)


def _percentil(ordenados: list[int], q: int) -> int:
    """Percentil por posto mais proximo, sem numpy."""
    if not ordenados:
        return 0
    i = round((q / 100.0) * (len(ordenados) - 1))
    return ordenados[min(max(i, 0), len(ordenados) - 1)]


def preparar_referencias(por_alvo: dict[str, list[str]]) -> dict:
    """{alvo: {"fps": [...], "n": int, "tamanho": (baixo, alto)}}.

    `tamanho` e a faixa de atomos pesados das referencias daquele alvo, entre
    os percentis `PERCENTIL_BAIXO` e `PERCENTIL_ALTO`.
    """
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    gen = gerador()
    saida = {}
    for alvo, lista in por_alvo.items():
        fps, tamanhos = [], []
        for smiles in lista:
            mol = Chem.MolFromSmiles(smiles or "")
            if mol is None:
                continue
            fps.append(gen.GetFingerprint(mol))
            tamanhos.append(mol.GetNumHeavyAtoms())
        tamanhos.sort()
        saida[alvo] = {
            "fps": fps, "n": len(fps),
            "tamanho": (_percentil(tamanhos, PERCENTIL_BAIXO),
                        _percentil(tamanhos, PERCENTIL_ALTO))
            if tamanhos else None}
    return saida


def tamanho_de(smiles: str):
    """Atomos pesados, ou None para SMILES ilegivel."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles or "")
    return None if mol is None else mol.GetNumHeavyAtoms()


def zonas_do_composto(fingerprint, tamanho, referencias: dict) -> dict:
    """{alvo_chembl_id: zona} para um composto, em todos os alvos com
    referencia."""
    return {alvo: zona(sim_top5(fingerprint, r["fps"]), r["n"],
                       tamanho=tamanho, faixa_de_tamanho=r["tamanho"])
            for alvo, r in referencias.items()}
