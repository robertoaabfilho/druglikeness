# Benchmark externo — o componente de alvo predito contra linhas de base triviais

Medido em 2026-11-27 num contêiner limpo. Nada aqui está commitado.

Versões em `versoes.txt`. Pesos: PP4CKD `ecfp4_fold_4.pt`, do tarball do Zenodo
([10.5281/zenodo.23147411](https://zenodo.org/records/23147411)), MD5
`f9380bb9cef3fe56d6d49d3e4a4d57ec` — **idêntico ao publicado, e truncado na
origem**: `tar` para com `Unexpected EOF` e só 23 dos ~77 `.pt` saem íntegros.
Rótulos: `data/target_labels.tsv` de `github.com/jixing475/pp4ckd` (não vem no
tarball). Dataset: LIT-PCBA, conjunto `full`, de
`drugdesign.unistra.fr` (via HTTPS; o link publicado é HTTP e não resolve).

## O que foi medido, e o que saiu

### 1. Poder estatístico do conjunto-controle (`poder_controle.py`)

Distribuição exata de Mann-Whitney para 6 "trata" × 6 "causa" — o desenho de
onde saem as AUCs do README. Só existem 924 arranjos e 37 valores possíveis de
AUC, de 1/36 em 1/36.

| AUC do README | o que é | p exato (uma cauda) |
|---|---|---:|
| 0,750 | só alvos humanos de dor | **0,0898 — não significante** |
| 0,806 | `p_humano × (1 − p_bact)` | **0,0465** |
| 1,000 | `CHEMBL360` isolado | 0,0011 |

**0,806 é exatamente o menor valor que alcança p ≤ 0,05 neste desenho.** Um
composto trocando de lugar com o vizinho o põe fora da significância. E como a
regra foi escolhida por maximizar a AUC entre cinco candidatas *no mesmo
conjunto*, o p nominal de 0,0465 não é o p da escolha: com cinco tentativas, o
p ajustado fica em torno de 0,2.

### 2. Linhas de base triviais no conjunto-controle (`benchmark_controle.py`)

| escore | AUC | p exato |
|---|---:|---:|
| `−p_Mtb` isolado (CHEMBL360) | **1,000** | 0,0011 |
| PP4CKD humano(max) × (1−p_Mtb) × (1−p_S.aureus) | 0,833 | 0,0325 |
| **BASE: −promiscuidade** (soma dos 7.676) | **0,778** | 0,0660 |
| PP4CKD alvos humanos (max / média / soma) | 0,750 | 0,0898 |
| **BASE: Tanimoto ao vizinho mais próximo dos "trata"** | **0,750** | 0,0898 |
| BASE: +promiscuidade | 0,222 | 0,9535 |

As duas linhas de base triviais **empatam ou ganham** do aparato de alvos.
Nenhuma delas usa lista de alvos, modelo ou conhecimento da doença.

Interpretação: no controle de 12 o que separa os grupos não é farmacologia de
alvo, é quimiotipo — os seis "causa" são antibióticos (promíscuos no modelo) e
os seis "trata" são fármacos de SNC. `CHEMBL360` tira 1,000 porque é
literalmente um classificador de antibacteriano.

> **Defeito achado de passagem:** a rede devolve **logits**, não
> probabilidades. `prever_alvos.py` aplica o sigmoide corretamente, mas o termo
> negativo `(1 − p)` só faz sentido com `p ∈ [0,1]`. Sobre logits ele inverte o
> ranking: no primeiro arame deste benchmark, antes da correção, a AUC caiu de
> 0,750 para 0,444. **Conferir se `mods/ranque/negativo.py::agregar` recebe
> probabilidade ou logit** — não pude inspecionar, está no `druglikeness2`.

### 3. LIT-PCBA (`benchmark_litpcba.py`, `resultado_litpcba.tsv`)

15 alvos, 10.033 ativos, 2.798.737 inativos (o distribuído é maior que os
7.844/407.381 que o artigo cita). 15 minutos em 4 CPUs.

`ESR1_ago` e `ESR1_ant` mapeiam para **o mesmo** `CHEMBL206`.

| alvo | PP4CKD AUROC | PP4CKD EF1% | 1-NN AUROC | 1-NN EF1% |
|---|---:|---:|---:|---:|
| ADRB2 | 0,654 | 11,76 | 0,667 | **35,30** |
| ALDH1 | 0,724 | 8,46 | 0,741 | 6,50 |
| ESR1_ago | 0,681 | **15,37** | 0,582 | 0,00 |
| ESR1_ant | 0,496 | 2,97 | 0,715 | **26,74** |
| FEN1 | 0,621 | 7,05 | 0,857 | **45,78** |
| GBA | 0,678 | 20,48 | 0,838 | **58,40** |
| IDH1 | 0,582 | 5,13 | 0,717 | **38,49** |
| KAT2A | 0,639 | 9,79 | 0,661 | **25,26** |
| MAPK1 | 0,649 | 5,85 | 0,702 | **11,04** |
| MTORC1 | 0,559 | 1,03 | 0,737 | **21,63** |
| OPRK1 | **0,952** | **70,84** | 0,860 | 66,69 |
| PKM2 | 0,794 | 25,45 | 0,831 | **36,47** |
| PPARG | 0,643 | 14,92 | 0,742 | **37,31** |
| TP53 | 0,504 | 1,28 | 0,748 | **39,68** |
| VDR | 0,620 | 3,28 | 0,795 | **34,15** |
| **média** | **0,653** | **13,58** | **0,746** | **32,23** |
| **mediana** | 0,643 | 8,46 | 0,741 | 35,30 |

**O 1-NN ganha em 13 dos 15 alvos.**

Subamostrar inativos para 60.000 não distorceu: as colunas "completo" e
"sub60000" do PP4CKD batem na terceira decimal.

## Três ressalvas que impedem comparação direta com a tabela publicada

Os números publicados em LIT-PCBA (Glide-SP 3,41 · Gnina 4,63 · DrugCLIP 5,51 ·
SPRINT 10,78 de EF1%) **não** são comparáveis de frente com os de cima:

1. **O 1-NN daqui é um teto, não um método.** Leave-one-out sobre o conjunto de
   ativos usa *outros ativos do próprio conjunto de teste* como referência —
   informação que um método zero-shot não tem. Ele mede **viés de analogia** do
   LIT-PCBA, que é o seu uso legítimo (é o probe de Chen et al. 2019), não
   desempenho implantável.
2. **O PP4CKD tem sobreposição com o ChEMBL.** É zero-shot quanto aos rótulos
   do LIT-PCBA, mas treinou no ChEMBL 36, e os ativos do LIT-PCBA vêm de
   ensaios do PubChem que se sobrepõem ao ChEMBL. O caso `OPRK1` é o exemplo:
   AUROC 0,952 e EF1% 70,84 sobre 24 ativos, num alvo com 5.626 ativos de
   treino. **Nenhum dos números do PP4CKD deve ser lido como superando o SOTA**;
   a média de 13,58 acima dos 10,78 do SPRINT é quase certamente vazamento mais
   diferença de protocolo.
3. **O conjunto distribuído difere do que o artigo descreve** em contagem de
   ativos e inativos, e EF depende da proporção.

## O que generalizou, e é o achado do benchmark

Contagem de ativos de treino do alvo **não** explica o desempenho aqui
(Spearman rho = −0,238, p = 0,39) — ao contrário do que o
`validar_alvos_medidos.py` encontrou no seu próprio conjunto. O que separa é a
**natureza do ensaio**:

| ensaio do LIT-PCBA | n | EF1% médio | mediana |
|---|---:|---:|---:|
| ligação / enzimático direto | 9 | **18,31** | 9,79 |
| fenotípico ou dependente de direção | 6 | **6,48** | 3,13 |

(Spearman rho = −0,472, p = 0,075.)

O modelo funciona onde a pergunta é "esta molécula se liga a esta proteína" e
falha onde a pergunta é "esta molécula produz este efeito na célula" — o que é
exatamente o que ele foi treinado para fazer, e exatamente o que um pipeline
terapêutico precisa que ele não faz.

### A cegueira de direção, agora com número

`ESR1_ago` e `ESR1_ant` compartilham a cabeça `CHEMBL206`: o escore que o
PP4CKD dá a um composto é **idêntico** nos dois conjuntos. Ainda assim o EF1%
sai 15,37 no agonista e 2,97 no antagonista — a diferença inteira vem da
composição dos conjuntos. Agonista e antagonista de ESR1 são indicações
**opostas**. A limitação que o README do `druglikeness2` declara em prosa
("a direção farmacológica não entra em nada, e não pode") fica medida: o
componente predito não pode ordenar por "trata" quando é a direção da modulação
que separa tratar de piorar.

## Convergência com a literatura e com o que o projeto já mediu

| fonte | número | o que diz |
|---|---|---|
| Rephetio (Himmelstein 2017) | prior só de grau: **97,9% → 54,1%** | prior trivial brilha in-sample e desaba fora |
| este projeto, controle de 12 | **0,806** | na fronteira exata de p = 0,05 |
| este projeto, 222 externos | **0,484** | canal negativo ao nível do azar |
| Chen et al. 2019 (DUD-E) | KNN **K=1** ótimo em 94% dos alvos | similaridade resolve o benchmark |
| PPB2 (Awale & Reymond) | NN(ECfp4) com a melhor recall em CV | vizinho mais próximo é linha de base forte |
| aqui, LIT-PCBA | 1-NN ganha em **13/15** | idem, em benchmark AVE-debiased |

O padrão do seu projeto é o do Rephetio, com os mesmos dois dígitos. Isso lhe dá
nome (**prior de grau / viés de promiscuidade**) e uma mitigação publicada:
incluir o escore trivial como termo explícito do modelo, para que o que se
reporte seja o ganho **acima** dele, não o total.

## Próximos benchmarks, por esforço

1. **Controle "actives-as-decoys"** (1 dia, CPU). Substituir inativos por
   ativos de alvos de função diferente. Modelo honesto cai para AUROC ≈ 0,500;
   referência: Chen et al., DUD-E 0,618 → 0,500 ([DOI](https://doi.org/10.1371/journal.pone.0220113)).
2. **Prior de grau isolado, estilo Rephetio** (1 dia). Escore só com nº de
   alvos do composto no ChEMBL + nº de fármacos já associados à doença, medido
   in-sample e externo. **Se o pipeline não bater claramente esse prior no
   externo, não há resultado.**
3. **repoDB com negativos reais** (2–3 dias). 6.677 pares aprovados e 4.123
   **falhados** ([DOI](https://doi.org/10.1038/sdata.2017.29)); recall@k
   leave-one-drug-out. Referências: Rephetio 85,5% (DrugCentral) e 70,0%
   (ensaios); CANDO 12,8% no top-10.
4. **Split temporal em Papyrus + MoleculeACE** (1 semana). Resposta direta à
   acusação de vazamento via ChEMBL.

## Reproduzir

```bash
pip install rdkit scipy scikit-learn
pip install torch --index-url https://download.pytorch.org/whl/cpu
python3 poder_controle.py
python3 benchmark_controle.py      # a partir da raiz do repositório
python3 benchmark_litpcba.py       # ~15 min em 4 CPUs
```

Caminhos de pesos, rótulos e dataset estão no topo de cada script.
