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
> probabilidades. O carregador do modelo aplica o sigmoide corretamente, mas o termo
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
a validação do próprio projeto encontrou no seu conjunto. O que separa é a
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

---

# Adendo — trocar o PP4CKD por outro modelo aberto?

Medido em 2026-11-28. Pergunta: vale substituir o PP4CKD pelo
`chembl_multitask_model` oficial do EMBL-EBI
(<https://github.com/chembl/chembl_multitask_model>, MIT, ONNX)?

## Cobertura de alvo: a troca apaga o canal negativo

O painel oficial é `SINGLE PROTEIN` por construção.

| conjunto de alvos | PP4CKD (7.676) | ChEMBL mt36 (786) |
|---|---|---|
| LIT-PCBA | 15/15 | 14/15 (falta KAT2A) |
| dor neuropática | 9/9 | **6/9** (falta α2δ-2, NMDA1, NMDA2B) |
| patógeno | 5/5 | **1/5** |

`CHEMBL360` (*M. tuberculosis*, célula inteira) **não existe** no painel do
ChEMBL. É de onde vem o único resultado do projeto com p sólido (AUC 1,000,
p = 0,0011) e é o alicerce do canal negativo, que é definido sobre alvo não
humano. **Substituir apaga os dois.**

## Desempenho: o modelo oficial fica no azar neste benchmark

LIT-PCBA, 14 alvos comuns. PP4CKD no conjunto completo; ChEMBL mt36 na
subamostra de 60.000 inativos (mostrado equivalente). `benchmark_multitask.py`,
`resultado_multitask.tsv`.

| alvo | PP4CKD AUROC | PP4CKD EF1% | mt36 AUROC | mt36 EF1% | 1-NN EF1% |
|---|---:|---:|---:|---:|---:|
| ADRB2 | 0,654 | 11,76 | 0,704 | 0,00 | 35,30 |
| ALDH1 | 0,724 | 8,46 | 0,489 | 0,98 | 6,50 |
| ESR1_ago | 0,681 | 15,37 | 0,533 | 0,00 | 0,00 |
| ESR1_ant | 0,496 | 2,97 | 0,564 | 0,99 | 26,74 |
| FEN1 | 0,621 | 7,05 | 0,439 | 2,71 | 45,78 |
| GBA | 0,678 | 20,48 | 0,517 | 1,81 | 58,40 |
| IDH1 | 0,582 | 5,13 | 0,637 | 2,57 | 38,49 |
| MAPK1 | 0,649 | 5,85 | 0,463 | 0,97 | 11,04 |
| MTORC1 | 0,559 | 1,03 | 0,482 | 3,09 | 21,63 |
| OPRK1 | 0,952 | 70,84 | 0,651 | 0,00 | 66,69 |
| PKM2 | 0,794 | 25,45 | 0,553 | 1,83 | 36,47 |
| PPARG | 0,643 | 14,92 | 0,647 | 3,73 | 37,31 |
| TP53 | 0,504 | 1,28 | 0,388 | 0,00 | 39,68 |
| VDR | 0,620 | 3,28 | 0,512 | 1,70 | 34,15 |
| **média** | **0,654** | **13,85** | **0,541** | **1,46** | **32,73** |
| **mediana** | 0,646 | 7,75 | 0,525 | 1,34 | 35,89 |

O PP4CKD ganha em **13 de 14** alvos. O mt36 tem EF1% ≤ 1 (azar ou pior) em
**7 de 14** e AUROC abaixo de 0,5 em **5 de 14**.

**Controle positivo: o uso do ONNX está correto.** Duloxetina dá p = 0,685 em
SERT e 0,703 em NET contra 0,032 no VDR (não-alvo); imatinibe 0,334 em ABL1,
estradiol 0,540 em ESR1, morfina 0,411 em OPRK1 — todos do lado certo. O
modelo funciona; ele fica no azar *neste benchmark*.

## Por que: o limiar do rótulo de treino, não a linhagem do modelo

| | PP4CKD | ChEMBL multitask |
|---|---|---|
| limiar de "ativo" | ≤ 10 µM **ou** inibição ≥ 50% | quinase ≤ 30 nM · GPCR/NR ≤ 100 nM · canal iônico ≤ 10 µM · resto ≤ 1 µM |

Os ativos do LIT-PCBA são confirmações de HTS dose-resposta do PubChem, na
faixa de µM. A definição de ativo do painel do ChEMBL exclui essa faixa quase
inteira; a do PP4CKD a inclui. `OPRK1` ilustra: GPCR com corte de 100 nM no
mt36 e EF1% = 0,00, contra 70,84 do PP4CKD.

**Consequência para o projeto:** o corpus de ligantes é de produtos naturais
(NPASS/COCONUT/LOTUS), cuja potência típica é de µM — o mesmo regime do
LIT-PCBA e o mesmo regime em que o mt36 falha. O limiar frouxo do PP4CKD é a
escolha **certa** para este uso. A procedência melhor do modelo oficial não
compensa um modelo cujo rótulo de treino não cobre a faixa de potência que
você tria.

## Recomendação

**Acrescentar como segunda opinião, não substituir.** Onde os dois cobrem o
alvo, concordância entre modelos independentes é o argumento de procedência
que o PP4CKD não tem sozinho (tarball truncado, folds não publicados, sem
revisão por pares) — e o `CHEMBL360` continua existindo. A costura é contida:
`mods/modelo.py` expõe quatro pontos de entrada (`carregar_modelo`,
`fingerprint_de`, `carregar_rotulos`, mais os auxiliares de organismo) e é
importado por quatro módulos de produção.

Nenhum dos dois conserta as duas falhas reais, que exigem consertos diferentes:

* **espaço de rótulos fixo** → modelo condicionado em sequência (ConPLex,
  SPRINT). É a única troca que compra capacidade nova; o teto já está medido
  no seu manifesto (61 de 259 proteínas do enriquecimento fora do modelo, 44%).
* **direção do efeito** → exige rótulo de direção, não preditor de ligação
  melhor. A metade barata é o *direction of effect* do Open Targets (8 fontes
  de evidência), que diz qual direção se quer — hoje não codificada em
  lugar nenhum do pipeline.

## Correção a um alerta anterior deste documento

A ressalva de que `mods/ranque/negativo.py::agregar` poderia receber logit em
vez de probabilidade **não se aplica**: `mods/ranquear.py:389` faz
`torch.sigmoid(rede(...))` antes de agregar. O que permanece é um risco de
documentação: `construir_rede` termina em `Linear` e devolve logit, enquanto o
docstring do módulo desenha "→ sigmoide". Quem reusar `construir_rede` fora do
`ranquear.py` erra — foi o que aconteceu no primeiro arame deste benchmark, e
custou uma AUC de 0,750 virar 0,444.

---

# Adendo 2 — fusão de canais e domínio de aplicabilidade, medidos por precisão no topo

Medido em 2026-11-28. Pergunta: investir em computação reduz o gasto de
bancada? A métrica, portanto, é `precisão@20` — não AUC, não EF1% sobre
centenas de milhares de compostos.

Divisão LIT-PCBA **AVE-unbiased**, que corrige um defeito do Adendo 1: ali o
1-NN era leave-one-out e um ativo de teste via outros ativos de teste. Aqui o
canal de similaridade usa `active_T` como referência e é avaliado em
`active_V` + `inactive_V`. 10 alvos com ≥ 15 ativos em V.
`benchmark_fusao.py`, `resultado_fusao.tsv`.

Canais: `p_alvo` (PP4CKD na cabeça do alvo), `sim_T` (maior Tanimoto aos
ativos de treino), `famil` (maior probabilidade entre os 7.676 alvos, proxy de
domínio de aplicabilidade). Fusões: média de postos percentuais, RRF (k = 60),
produto e máximo dos postos.

## As duas alavancas mais baratas da literatura não deram nada

| canal | prec@20 | prec@50 | prec@100 | EF1% | AUROC |
|---|---:|---:|---:|---:|---:|
| `p_alvo` | 0,100 | 0,078 | 0,065 | 6,50 | 0,600 |
| `sim_T` | 0,060 | 0,066 | 0,048 | 5,41 | 0,555 |
| `famil` | 0,010 | 0,012 | 0,011 | 1,45 | 0,562 |
| média de postos | 0,105 | 0,098 | 0,075 | 4,42 | 0,598 |
| RRF | 0,095 | 0,084 | 0,071 | 6,65 | 0,610 |
| produto | **0,110** | **0,100** | 0,075 | 4,42 | 0,583 |
| máximo | 0,080 | 0,076 | 0,069 | 6,77 | 0,610 |

A melhor fusão dá prec@20 = 0,110 contra 0,100 do PP4CKD sozinho. Em 20
compostos × 10 alvos são 200 vagas, logo **+2 acertos no total** — ruído.

Filtro de domínio (descarta a metade de menor `famil`): `p_alvo` vai de 0,100
para 0,110 em prec@20 (+0,010, ruído) e de 6,50 para **8,77 em EF1%** (+35%
relativo, real mas irrelevante para quem ensaia 20 compostos). O canal `famil`
isolado é inútil como ranqueador (prec@20 = 0,010) — serve como filtro, não
como escore.

## O número que responde à pergunta

| alvo | n | ativos | taxa-base | prec@20 | ganho |
|---|---:|---:|---:|---:|---:|
| ALDH1 | 26.640 | 1.343 | 5,041% | 0,65 | 13× |
| ESR1_ant | 819 | 25 | 3,053% | 0,05 | 2× |
| FEN1 | 87.038 | 91 | 0,105% | 0,05 | 48× |
| GBA | 72.240 | 41 | 0,057% | 0,00 | 0× |
| KAT2A | 84.710 | 48 | 0,057% | 0,05 | 88× |
| MAPK1 | 15.327 | 77 | 0,502% | 0,05 | 10× |
| MTORC1 | 8.267 | 24 | 0,290% | 0,00 | 0× |
| PKM2 | 61.143 | 136 | 0,222% | 0,10 | 45× |
| TP53 | 811 | 16 | 1,973% | 0,05 | 3× |
| VDR | 65.256 | 165 | 0,253% | 0,00 | 0× |

`prec@20` ≤ 0,15 em **9 de 10 alvos**. Sem o ALDH1, cuja taxa-base é 50× a dos
demais, a média dos 8 alvos restantes é **0,038 — 0,8 acerto esperado em 20
compostos ensaiados**.

O ganho relativo sobre a taxa-base é real e grande (mediana 6×, até 88×), mas
**zero em 3 dos 10 alvos**. E a taxa-base não explica estatisticamente o
`prec@20` (Spearman rho = 0,361, p = 0,31) — a hipótese de que tudo se reduz a
taxa-base **não se sustenta**; com n = 10 a variância é que domina.

## Consequência para a estratégia

As duas alavancas computacionais mais barata e melhor documentadas da área —
fusão de canais independentes e domínio de aplicabilidade — foram medidas
sobre o modelo real deste pipeline e **não movem a precisão no topo**. Isso é
informação: mais computação *desse tipo* não é onde está o ganho.

O que a evidência acumulada sustenta, em ordem:

1. **Maximizar a fração do top-k com evidência MEDIDA.** A faixa `1a` já
   ordena medido acima de previsto. Quantificar quantos compostos o top-k
   consegue preencher só com `pchembl` do ChEMBL, antes de recorrer a
   previsão, é contagem — não precisa de modelo nem de validação — e é onde a
   precisão de fato existe.
2. **Encolher o conjunto candidato antes de ranquear, não melhorar o
   ranqueador.** Com 107.123 ligantes e taxa-base de 0,1%, o top-20 rende ~1
   acerto. A alavanca é o denominador.
3. **Especialistas de cobertura e direção** continuam justificados — não por
   melhorarem o ranqueamento, mas por ampliarem quais compostos conseguem
   receber evidência de faixa 1 ou 2.

O que a evidência **não** sustenta: acrescentar preditor esperando
ranqueamento melhor (o painel oficial do ChEMBL ficou no azar; fusão deu
ruído), e finetuning antes de existir instrumento de medição.

---

# Adendo 3 — onde o canal predito vale: o estrato de novidade

Medido em 2026-11-28. Pergunta levantada corretamente pelo autor do projeto:
a novidade de um reposicionamento está nas **previsões**, não na evidência
medida — evidência medida contra o alvo da doença é, por definição,
não-novidade. Preencher o top-k com ela otimizaria precisão gastando
exactamente o que se quer descobrir.

Isso torna a pergunta decisiva: **o canal predito enriquece em compostos
estruturalmente distantes dos ativos conhecidos?** Se os acertos estiverem
todos na vizinhança química do que já se conhece, o 1-NN pega os mesmos e o
modelo não reposiciona nada.

## Desenho

LIT-PCBA AVE-unbiased, 10 alvos. Para cada molécula do conjunto V,
`sim_T` = maior Tanimoto aos ativos de **treino**. Estratifica por `sim_T` e,
**dentro** de cada estrato — comparando ativos e inativos igualmente distantes
do treino — mede se `p_alvo` ainda separa. Estratificar assim é o controle de
viés de analogia no espírito de Chen et al. 2019: um modelo que só memoriza
vizinhança cai para AUROC 0,5 no estrato distante.
`benchmark_novidade.py`, `resultado_novidade.tsv`.

## O canal predito sobrevive ao controle; o de similaridade não

| faixa `sim_T` | n | ativos | taxa-base | AUROC `p_alvo` | AUROC `sim_T` |
|---|---:|---:|---:|---:|---:|
| 0,00–0,25 | 175.877 | 241 | 0,137% | **0,552** | 0,402 |
| 0,25–0,35 | 189.770 | 540 | 0,285% | 0,637 | 0,250 |
| 0,35–0,45 | 41.979 | 547 | 1,303% | 0,645 | 0,195 |
| 0,45–0,60 | 12.379 | 424 | 3,425% | 0,631 | 0,321 |
| 0,60–1,01 | 2.246 | 214 | 9,528% | 0,719 | 0,409 |

O `p_alvo` fica acima do acaso em **todos** os estratos, inclusive o mais
distante. O `sim_T` fica **abaixo** de 0,5 dentro de cada faixa — era de
esperar: controlada a vizinhança, a similaridade não carrega mais informação.

## O estrato de novidade, quantificado

Ranqueando **dentro** do estrato `sim_T < 0,25` (que é o uso operacional:
restringir aos compostos novos e ordená-los entre si), por posto percentual
dentro de cada alvo:

```
n = 175.877   ativos = 241   taxa-base = 0,137%
AUROC = 0,5722   EP 0,0192   IC 95% [0,5345 – 0,6100]   z = 3,76
p (uma cauda, H0: AUROC = 0,5) = 8,7 x 10^-5
```

| medida | valor | ganho sobre a taxa-base |
|---|---:|---:|
| EF 0,1% | 12,44× | — |
| EF 0,5% | 4,15× | — |
| precisão@20 | 1/20 = 5,00% | **36×** |
| precisão@50 | 2/50 = 4,00% | 29× |
| precisão@100 | 3/100 = 3,00% | 22× |
| precisão@200 | 4/200 = 2,00% | 15× |

**Tradução para bancada:** no regime de novidade, ensaiar 20 compostos
escolhidos pelo pipeline rende ~1 acerto, contra 0,027 de escolha aleatória.
É a diferença entre testar 20 compostos e testar ~730 para achar o primeiro.
E 5% está dentro da faixa das campanhas prospectivas publicadas (CACHE #1:
3,7%; equipe vencedora: 8,5%).

## Isso reformula a conclusão do Adendo 1

O Adendo 1 concluiu que o 1-NN ganha do PP4CKD em 13 de 15 alvos e que o
modelo não bate a linha de base trivial. **A medida estava certa e respondia à
pergunta errada** para este projeto. A explicação das duas coisas juntas:

* os ativos do LIT-PCBA são majoritariamente análogos de ativos conhecidos,
  então no agregado o 1-NN ganha;
* no estrato distante, o 1-NN está no acaso ou abaixo, e o `p_alvo` está em
  0,572 com p < 10⁻⁴.

Ou seja: **o valor do modelo é precisamente a fração não-análoga**, que é a
única que interessa a reposicionamento. Julgar o canal predito por EF1%
agregado subestima exatamente o uso pretendido.

## Consequência de desenho

Não misturar os dois regimes num ranking só. Duas listas, pela mesma lógica
de sub-faixa `a`/`b` que o projeto já usa para procedência de alvo, estendida
à novidade **química**:

1. **Lista de confirmação** — evidência medida. Alta precisão, novidade zero.
   Serve para calibrar o canal predito nos alvos da doença e como controle
   positivo de um ensaio.
2. **Lista de reposicionamento** — predita, restrita a `sim_T` baixo contra os
   ativos conhecidos do alvo. É onde está a novidade, com ~5% de acerto
   esperado no top-20 e 36× sobre o acaso.

O limiar de `sim_T` passa a ser parâmetro de primeira classe do ranqueamento,
ao lado do corte de probabilidade — e tem de ir para o manifesto, porque muda
o resultado.

Pendência honesta: a fusão de canais e o filtro de domínio foram medidos no
agregado (Adendo 2) e não no estrato de novidade, onde a linha de base é
outra. Remedir ali é barato e ainda não foi feito.

---

# Adendo 4 — descritor QSAR/físico-químico no estrato de novidade

Medido em 2026-11-28. Proposta avaliada: usar parâmetros QSAR/ADMET também na
previsão de similaridade. `benchmark_descritor.py`,
`resultado_descritor.tsv`.

`sim_desc` = proximidade ao ativo de treino mais próximo no espaço de 12
descritores do RDKit (MW, logP, TPSA, HBD, HBA, ligações rotáveis, anéis
aromáticos, fração Csp3, átomos pesados, QED, MolMR, anéis), padronizados por
z-score sobre o conjunto de triagem. O preditor de ADMET por serviço externo
não entra: é API e
rejeita 30–35% das moléculas.

Avaliado **dentro** do estrato de novidade (`sim_T < 0,25` em ECFP4), nos 10
alvos. n = 175.877, 241 ativos, taxa-base 0,137%.

| canal | AUROC | EF 0,1% | prec@20 |
|---|---:|---:|---:|
| `p_alvo` (PP4CKD) | **0,5722** | **12,44×** | **5,00%** |
| `sim_desc` (descritor) | 0,4953 | 0,00 | 0,00% |
| fusão média de postos | 0,5509 | 0,00 | 0,00% |
| fusão produto de postos | 0,5410 | 0,00 | 0,00% |

**O canal de descritor está exatamente no acaso (0,4953) e a fusão PIORA o
modelo** — AUROC cai de 0,5722 para 0,5509 e a precisão@20 cai de 5,00% para
0,00%. Isso também fecha a pendência do Adendo 2: a fusão não ajuda nem no
estrato onde a linha de base é outra.

## O uso que sobrevive: definir novidade, não pontuar

| | n | % |
|---|---:|---:|
| ECFP4-distante **e** descritor-perto | 73.388 | **41,7%** |
| ECFP4-distante **e** descritor-distante | 102.489 | 58,3% |

Dos compostos que o ECFP4 chama de novos, **41,7% são vizinhos em perfil
físico-químico** — a definição de novidade por fingerprint sozinho é grosseira
e admite como novo quase metade de algo que é, em propriedade, mais do mesmo.

E a consequência é precisa: como `sim_desc` tem AUROC 0,4953, ele **não tem
viés de atividade**. Filtrar por ele não melhora nem piora a taxa de acerto —
encolhe o conjunto sem custo de precisão e torna a afirmação de novidade mais
forte. É de graça como refinamento de novidade, e inútil como escore.

## Veredito sobre a proposta

Três usos distintos, com evidências opostas:

1. **Prever eficácia / pontuar atividade** — não sustentado, e agora duas
   vezes: o QED deu 0,384 (abaixo do acaso) nos 222 compostos externos do
   projeto, e aqui o `sim_desc` dá 0,4953 com a fusão degradando o modelo.
2. **Definir o conjunto candidato (novidade)** — sustentado e barato. 41,7%
   do estrato muda de classificação, sem custo de precisão.
3. **Desempate por desenvolvibilidade a evidência de alvo fixa** — defensável,
   não medido aqui. A literatura diz que desejabilidade empata com
   "similaridade ao que já é bom" ([DOI](https://doi.org/10.1039/D5DD00387C)),
   isto é, não ganha do trivial mas não atrapalha como ordenação secundária.

## Onde a sessão chegou

Quatro intervenções medidas sobre o ranqueamento, três negativas:

| intervenção | resultado |
|---|---|
| trocar o PP4CKD pelo painel oficial do ChEMBL | pior em 13 de 14 alvos |
| fundir canais (similaridade + alvo) | ruído no agregado, **degrada** no estrato de novidade |
| filtro de domínio de aplicabilidade | +0,010 em prec@20 (ruído) |
| **restringir o canal predito ao estrato de novidade** | **AUROC 0,5722, p = 8,7×10⁻⁵, prec@20 36× a taxa-base** |

A única coisa que acrescentou foi **não mudar o modelo, e sim mudar a pergunta
feita a ele**.

---

# Adendo 5 — calibração nos alvos de dor: o vazamento, e o holdout temporal que o contorna

Medido em 2026-11-28. `calibrar_dor.py`, `calibrar_temporal.py`,
`resultado_calibracao_dor.tsv`, `resultado_calibracao_temporal.tsv`.

## Primeira tentativa, e o defeito dela

Desenho: ativos = ChEMBL `pchembl >= 6`, ensaio B/F; divisão temporal por ano
do primeiro documento; decoys = os 10.033 ativos do LIT-PCBA (bioativos,
drug-like, alvos sem relação com dor) — o controle *actives-as-decoys* de
Chen et al. 2019.

| alvo | estrato | ativos | taxa-base | AUROC | prec@20 |
|---|---|---:|---:|---:|---:|
| α2δ-1 | tudo | 242 | 2,42% | **1,0000** | 100% |
| α2δ-1 | novidade | 236 | 2,51% | **1,0000** | 100% |
| SERT | tudo | 1.030 | 9,53% | 0,9971 | 100% |
| SERT | novidade | 31 | 0,56% | 0,9701 | 95% |
| NET | tudo | 657 | 6,30% | 0,9970 | 100% |
| NET | novidade | 48 | 0,81% | 0,9884 | 100% |
| TRPV1 | novidade | 8 | 0,14% | 0,9997 | 40% |

**AUROC 1,0 não é desempenho, é memorização.** A divisão temporal protege o
canal de similaridade e **não** o modelo: o PP4CKD treinou em ChEMBL 36 (2025),
então todo ativo publicado depois de 2015 está dentro do treino dele.
(Parou em 4 de 9 alvos: a API do ChEMBL derrubou a conexão; nova tentativa
foi adicionada depois.)

**O contraste mede o vazamento.** Mesmo modelo, mesma estratificação, mesmo
tipo de decoy:

| origem dos ativos | AUROC no estrato de novidade |
|---|---:|
| PubChem HTS (LIT-PCBA, parcialmente fora do ChEMBL) | **0,572** |
| ChEMBL (dentro do corpus de treino) | **0,97 – 1,00** |

Consequência prática: **qualquer validação do canal predito contra
bioatividade do ChEMBL é inútil**, inclusive o recall@10 de 0,24 do
`validar_modelo.py` — que o próprio `validacao.alvos_medidos.md` já declara
medir memorização e generalização juntas. O 0,572 do LIT-PCBA é o número
honesto disponível hoje.

## A correção: modelo de era anterior

Modelo `chembl_28_multitask.onnx` (ChEMBL 28, jan/2021, 616 alvos) contra
ativos cujo **primeiro documento é de 2022 ou depois** — compostos que não
existiam no treino dele. Referência de similaridade: ativos até 2020. Dos 9
alvos de dor, 3 estão no painel de 616: SERT, NET e receptor mu.

| alvo | estrato | n | ativos | taxa-base | AUROC | EF1% | prec@20 | ganho |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| SERT | tudo | 9.999 | 219 | 2,19% | 0,803 | 15,52 | 40% | 18× |
| SERT | novidade | 3.166 | 9 | 0,28% | 0,662 | 0,00 | 0% | 0× |
| NET | tudo | 9.916 | 136 | 1,37% | 0,695 | 4,42 | 5% | 4× |
| NET | novidade | 3.713 | 8 | 0,22% | 0,494 | 12,54 | 0% | 0× |
| MOR | tudo | 10.273 | 493 | 4,80% | 0,935 | 19,42 | 100% | 21× |
| MOR | novidade | 3.283 | 12 | 0,37% | 0,885 | 49,74 | 20% | 55× |

Agregado do estrato de novidade (posto percentual dentro de cada alvo):

```
n = 10.162   ativos = 29   taxa-base = 0,285%
AUROC = 0,7077   EP 0,0544   IC 95% [0,6010 – 0,8144]   z = 3,82
p (uma cauda) = 6,8 x 10^-5
precisao@20  = 3/20  = 15,00%  (ganho 53x)
precisao@50  = 4/50  =  8,00%  (ganho 28x)
precisao@100 = 7/100 =  7,00%  (ganho 25x)
```

## Como ler isso, com a ressalva do mesmo tamanho do resultado

O sinal é real: p = 6,8 × 10⁻⁵ contra o acaso, em holdout temporal sem
vazamento, com o controle de decoy mais severo que existe.

Mas o estrato de novidade tem **29 ativos** em três alvos, e o agregado é
puxado pelo receptor mu. Por alvo, a precisão@20 é 0%, 0% e 20%. O intervalo
de confiança da AUROC vai de 0,60 a 0,81. **15% é estimativa de ponto com
incerteza grande, não número para fechar orçamento.** A faixa defensável é
"entre 0 e 20%, provavelmente na casa dos 10%", e isso é consistente com as
campanhas prospectivas publicadas (CACHE #1: 3,7%; equipe vencedora: 8,5%) —
sem permitir afirmar que seja melhor.

## Dois achados laterais que valem mais que o número principal

**1. A hipótese do limiar de potência se confirmou.** O Adendo 4 atribuiu o
fracasso do painel oficial do ChEMBL no LIT-PCBA a desencontro de limiar: ele
define ativo em ≤ 30 nM a ≤ 1 µM, e os ativos do LIT-PCBA são confirmações de
HTS em faixa de µM. Aqui os ativos são `pchembl >= 6` (≤ 1 µM), dentro do
regime dele — e **o mesmo modelo que ficou no azar no LIT-PCBA chega a AUROC
0,935 no receptor mu**. A previsão se sustentou.

Isso reabre a questão do modelo com uma resposta mais precisa que "não troque":
**o modelo certo depende do regime de potência que se tria.** Para alvo maduro
e candidato na casa de nM–µM, o painel oficial funciona e traz procedência
melhor, ONNX e nenhuma dependência de torch. Para produto natural em faixa de
µM, o limiar frouxo do PP4CKD é o adequado. É parâmetro de desenho, não
disputa de qualidade.

**2. Novidade química em alvo maduro é rara.** Dos 848 ativos medidos pela
primeira vez a partir de 2022 nesses três alvos, só **29 (3,4%)** são
estruturalmente distantes do que já se conhecia em 2020. O campo publica
majoritariamente análogos. Isso limita por cima quanto reposicionamento há
para achar em alvo bem estudado — e é argumento a favor de apontar o pipeline
a alvos menos explorados, onde a fração nova é maior.

---

# Adendo 6 — o mesmo holdout temporal em 145 alvos: o intervalo fecha, e para baixo

O Adendo 5 chegou a AUROC 0,7077 no estrato de novidade com **29 ativos**, e eu
disse na hora que o intervalo de confiança [0,601, 0,814] era largo demais para
decidir nada. Este adendo repete a receita idêntica em escala.

## Desenho

Igual ao Adendo 5, só maior. Mesmo modelo de era anterior
(`chembl_28_multitask.onnx`, 616 alvos, treino até 2020), mesmos decoys (os
ativos do LIT-PCBA sem medida no alvo), mesmo corte de novidade (Tanimoto
ECFP4 máximo aos ativos de referência < 0,25), mesma divisão temporal
(referência ≤ 2020, teste ≥ 2022).

Duas decisões de amostragem, declaradas antes de rodar:

- **Amostra aleatória de 150 entre os 322 alvos elegíveis** (≥ 50 ativos
  pós-2022), semente 20261128 — não os 150 maiores. Tomar os maiores
  enviesaria para alvo muito estudado e repetiria o defeito do Adendo 5, em que
  o agregado inteiro vinha do receptor mu.
- **Guarda de tamanho**: alvo com mais de 15 mil atividades é pulado. Na
  prática nenhum foi (0 pulados), porque a amostra aleatória não cai nos
  gigantes.

Sondagem prévia (`sondar_alvos.py`, uma requisição por alvo): dos 616 do
painel, **361 têm ≥ 30 ativos pós-2022**, somando 122.118 ativos recentes.
Coletados 145 alvos com dado suficiente nos dois lados do corte.

## O resultado principal

| | Adendo 5 (3 alvos) | Adendo 6 (145 alvos) |
|---|---|---|
| ativos no estrato de novidade | 29 | **1.769** |
| AUROC | 0,7077 | **0,6322** |
| IC95% | [0,601, 0,814] | **[0,618, 0,646]** |
| erro padrão | 0,0544 | 0,0071 |
| z | 3,82 | 18,50 |
| p (uma cauda) | 6,8e-05 | 9,8e-77 |

O intervalo encolheu de 21 pontos para 2,8, e o ponto caiu para perto do piso
do intervalo antigo. **0,7077 era otimismo de n pequeno. O número é 0,632.**

Segue acima do azar com folga enorme, e segue valendo o contraste do Adendo 3 —
no mesmo estrato o canal de similaridade fica abaixo do azar. Mas 0,63 é fraco
em termos absolutos, e a média esconde o que importa.

## Precisão no topo, medida por alvo

O agregado usa posto percentual dentro de cada alvo. Com 145 alvos, o topo de
cada um empata em 1,0 exato, e o desempate pessimista de `metricas_vs` empurra
todo decoy empatado para a frente — o que **zera precisão@k para k < 145 por
construção**. O `precisao@20 = 0%` que sai no agregado é limite pessimista, não
resultado; ficou rotulado assim no log para não enganar ninguém depois. A
AUROC não é afetada (usa postos médios).

A medida correta é dentro de cada alvo — e é a pergunta real: dos *k* que eu
ensaiaria neste alvo, quantos são ativos de verdade? Somada sobre os alvos:

| | análogo + novo (145 alvos) | só novidade (88 alvos) |
|---|---|---|
| precisão@5 | **52,6%** (23×) | **4,3%** (12×) |
| precisão@10 | 49,6% (22×) | 3,4% (9×) |
| precisão@20 | 45,4% (20×) | 2,8% (8×) |
| alvos com ≥ 1 acerto no top-20 | 122/145 | **19/88** |
| taxa-base média | 2,2% | 0,36% |
| AUROC mediana por alvo | 0,850 | 0,638 |
| EF1% mediana por alvo | 14,07 | **0,00** |

## Como ler

**No regime análogo o pipeline funciona de verdade.** Top-5 de um alvo: metade
acerta. 20× acima da taxa-base, consistente em 145 alvos independentes, sob
holdout temporal e com decoy bioativo (`actives-as-decoys`), que é o controle
mais severo da literatura. Isso não é artefato.

**No estrato de novidade o ganho é real mas pequeno em termos absolutos.** 4,3%
no top-5 significa ~23 ensaios por acerto estruturalmente novo. O ganho de 12×
é sólido; 12× sobre 0,36% ainda é 4%.

**E a dispersão importa mais que a média.** Só **19 dos 88 alvos** têm algum
acerto no top-20 — em 69 alvos, os 20 primeiros candidatos novos são todos
falsos. A mediana de EF1% é **zero**: para metade dos alvos, o top 1% não
contém nenhum ativo novo. O retorno não está espalhado, está concentrado em
cerca de um quinto dos alvos.

## Consequência para a decisão

Se novidade estrutural é o objetivo (reposicionamento), o custo esperado é ~23
ensaios por acerto, não 2. Duas saídas, não excludentes:

1. **Usar o regime análogo para selecionar e a novidade como filtro posterior**,
   não como critério de ranqueamento. 45% no top-20 paga a bancada; dali se
   escolhe o mais distante estruturalmente entre os que pontuaram bem.
2. **Prever qual alvo está no quinto bom.** Os 19 alvos que acertam não devem
   ser sorteio — devem ter assinatura (volume de treino, diversidade dos ativos
   de referência, classe de proteína). É mensurável com o TSV deste adendo, sem
   rede e sem bancada. Se houver sinal, o pipeline passa a saber de antemão
   quando não confiar na própria previsão de novidade, o que vale mais que subir
   0,63 para 0,65.

## Limitação declarada

Isto mede o modelo do ChEMBL 28 sobre um painel de 616 alvos **maduros**, não o
PP4CKD. Para o PP4CKD não existe número sem vazamento usando dado do ChEMBL
(Adendo 5). E o painel ser de alvo maduro significa que a pergunta sobre alvo
pouco explorado — onde a fração de novidade é maior — **este experimento não
responde**.

Arquivos: `sondar_alvos.py`, `sondagem_alvos.json`, `holdout_escalado.py`,
`resultado_holdout_escalado.tsv` (233 linhas: 145 alvos × estrato).

---

# Adendo 7 — um teste de falha: 221 antibióticos contra alvos imunes

Pergunta que originou isto: usar uma lista de ~200 antibióticos para buscar
aplicação no mimetismo molecular que leva à síndrome de Guillain-Barré. É o
adendo mais informativo do conjunto, porque mede **como o pipeline falha quando
a pergunta está fora do alcance dele** — e porque três das minhas previsões
sobre essa falha estavam erradas.

## Por que a pergunta não é para este pipeline

O mecanismo da GBS pós-*Campylobacter* é: o lipooligossacarídeo bacteriano imita
gangliosídeos (GM1, GD1a, GQ1b) → IgG anti-gangliosídeo → ativação de
complemento na membrana nodal. **O objeto do mimetismo é um glicolipídeo e o
agente do dano é um anticorpo.** Nenhum dos dois é proteína com bolsão de
ligação, e nenhuma das 616 cabeças (nem das 7.676 do PP4CKD) corresponde a
epítopo de gangliosídeo ou a idiotipo. Não há coluna para indexar.

Os nós a jusante que são proteína e têm química de molécula pequena: C5aR1
(precedente: avacopan), fator D, FcRn, BTK, S1PR1 (fingolimod funciona na
neurite autoimune experimental), MMP-9. Conferido no painel: **C5aR1, S1PR1,
MMP-9, BTK, JAK1, JAK2 e TNF-α estão presentes; TLR4, que seria o sensor do
próprio LOS e o mais próximo da mimetização, não está.**

Há ainda a janela temporal: a GBS começa 1–3 semanas depois da enterite, com a
resposta de anticorpo já montada. Erradicar a bactéria nesse ponto não desfaz
nada. A única hipótese defensável na lista é **efeito imunomodulador fora do
alvo antibacteriano** — minociclina e doxiciclina (inibição de MMP, supressão
microglial, literatura em neurite autoimune experimental) e macrolídeos.

## Desenho

221 antibióticos da classe ATC **J01** do ChEMBL contra 600 ativos do LIT-PCBA
amostrados como referência de "bioativo qualquer". Para cada composto:
similaridade ao domínio (Tanimoto ECFP4 máximo aos 201.196 ativos de referência
≤ 2020 dos 145 alvos do Adendo 6), familiaridade (`p_max`, a probabilidade
máxima entre as 616 cabeças — o proxy de domínio que o projeto usa),
probabilidade nas sete cabeças imunes, e descritores.

## Três previsões minhas, erradas

**Previ que o topo seria fluoroquinolona e tetraciclina, "aromáticos planos e
lipofílicos".** O topo é glicopeptídeo e cefalosporina:
quinupristina/dalfopristina 1,000, oritavancina 0,999, telavancina 0,999,
mezlocilina 0,998, ceftobiprol 0,998, vancomicina 0,995. Ciprofloxacina sai em
**188/221**, metronidazol em 217/221.

**Previ que o ranque seria explicado por lipofilicidade.** O que explica é
**tamanho**: Spearman de `p_max` com peso molecular **+0,407**, com logP só
+0,224. Peso médio dos 20 primeiros **770 Da**, dos 20 últimos **374 Da**.
Molécula maior liga mais bits do ECFP4, mais bits disparam mais cabeças, `p_max`
sobe. É artefato de contagem de bits.

**Previ que a lista cairia no estrato de novidade.** Medido: 7% abaixo de 0,25
(mediana 0,407). Mas o desenho não é comparável e eu não deveria ter
apresentado isso como teste: aqui o Tanimoto é contra a **união** de 201.196
referências, enquanto o Adendo 6 usa as referências **de cada alvo** (mediana
~1.300). Contra 201 mil referências quase nada cai abaixo de 0,25. O domínio
segue em aberto.

## O que o resultado mostra de verdade, e é pior

**1. O canal de familiaridade está saturado.** Mediana de `p_max` dos
antibióticos: **0,968**. Dos 600 bioativos quaisquer: **0,968**. Idênticas. A
biblioteca inteira ocupa 0,77–1,00. O proxy de domínio de aplicabilidade não
distingue um glicopeptídeo de um ligante de receptor humano. Isso explica
retroativamente por que o domínio de aplicabilidade deu ruído no Adendo 2: não
tem variância para dar.

**2. O filtro de novidade não vê o viés.** Spearman de `p_max` com similaridade
ao domínio: **+0,051**. O viés de tamanho passa por baixo dos dois filtros do
pipeline simultaneamente.

**3. As previsões nos alvos imunes são absurdas e sistematicamente infladas.**
Vancomicina com p = 0,92 em BTK; oritavancina 0,94; quinupristina 0,996 em
S1PR1 e 0,95 em JAK1. Glicopeptídeo não inibe quinase — não atravessa membrana.
E a inflação é sistemática (mediana antibiótico vs bioativo qualquer): C5aR1
**0,48 vs 0,17**, BTK **0,21 vs 0,07**, JAK1 **0,34 vs 0,17**, JAK2 0,35 vs
0,17. A triagem **pareceria produtiva**.

**4. O ranque aponta para o lado errado.** Minociclina, o único candidato com
literatura em neurite autoimune experimental, sai em **200/221**. Doxiciclina em
118. Os macrolídeos sobem (azitromicina 59, claritromicina 111) por serem
grandes, não por imunomodulação.

## Veredito e o que ele corrige no pipeline

Rodar essa lista produziria um top-20 plausível à vista e sem conteúdo: ordenado
por peso molecular, inflado justamente nos alvos imunes que se escolheria, e
invisível aos dois filtros do projeto. Não é resultado fraco — é resultado
enganoso, que custa bancada.

Duas correções que valem para **qualquer** uso do pipeline, e que só
apareceram por causa desta pergunta:

1. **Trocar o proxy de domínio.** `p_max` está saturado. Substitutos: a
   similaridade por alvo (que o Adendo 6 já calcula) ou `p_max` normalizado pela
   contagem de bits ligados.
2. **Corrigir ou declarar o viés de tamanho.** Medir a correlação de cada cabeça
   com o número de bits ligados e usá-la como linha de base a subtrair. Cabeça
   explicada por tamanho não deveria entrar no ranque.

Arquivos: `sondar_antibioticos.py`, `lista_j01.json`,
`resultado_antibioticos.tsv` (821 linhas).

---

# Adendo 8 — as três correções: uma funcionou, uma me refutou, uma achou regra prospectiva

O Adendo 7 terminou com duas correções propostas e o Adendo 6 com uma pergunta
aberta. Este adendo mede as três no mesmo dado (`corrigir.py`, 145 alvos do
cache do Adendo 6, 36.439 moléculas pontuadas uma vez). **Duas das três deram o
oposto do que eu havia proposto**, e são justamente as que mudam o desenho.

## Correção 1 — o proxy de domínio: funcionou, e o vencedor não era o proposto

Um proxy de domínio de aplicabilidade não deve prever atividade; deve prever
**confiança**. Logo o critério certo não é a AUROC do proxy, é: parte as
moléculas do estrato de novidade na mediana do proxy e compara a AUROC **do
modelo** nas duas metades. Proxy que serve tem AUROC maior na metade alta.

| proxy | AUROC metade alta | metade baixa | delta | z |
|---|---|---|---|---|
| **média das 5 maiores similaridades** | **0,6819** | **0,4987** | **+0,1832** | **+12,12** |
| similaridade máxima | 0,6763 | 0,5326 | +0,1437 | +9,62 |
| nº de bits do ECFP4 | 0,6402 | 0,5137 | +0,1265 | +7,14 |
| `p_max` (o atual) | 0,6545 | 0,5898 | +0,0647 | +4,20 |
| nº de cabeças > 0,9 | 0,6345 | 0,5847 | +0,0497 | +2,90 |
| `p_max` corrigido por tamanho | 0,6387 | 0,6228 | +0,0159 | +1,10 |

**Substituto: a média das 5 maiores similaridades ECFP4 aos ativos de
referência do alvo.** Separa quase três vezes melhor que `p_max`, e a metade de
baixo fica **no azar exato (0,4987)** — metade do estrato de novidade não tem
informação nenhuma, e dá para identificá-la antes de ensaiar.

Ressalva declarada: o estrato foi *selecionado* por similaridade máxima < 0,25,
então a similaridade é variável truncada dentro dele. O achado não é circular, e
a leitura operacional é que **o corte de 0,25 é generoso demais**.

### Onde o sinal morre, por decil

| decil | sim_top5 média | AUROC | EP | ativos |
|---|---|---|---|---|
| 1 | 0,1445 | 0,4399 | 0,0304 | 82 |
| 2 | 0,1713 | 0,4565 | 0,0341 | 67 |
| 3 | 0,1777 | 0,5228 | 0,0299 | 96 |
| 4 | 0,1835 | 0,4898 | 0,0248 | 134 |
| 5 | 0,1950 | 0,5458 | 0,0241 | 152 |
| 6 | 0,2016 | 0,6451 | 0,0238 | 159 |
| 7 | 0,2073 | **0,7050** | 0,0214 | 189 |
| 8 | 0,2130 | 0,6934 | 0,0201 | 217 |
| 9 | 0,2198 | 0,6931 | 0,0183 | 261 |
| 10 | 0,2298 | 0,6556 | 0,0148 | 412 |

A transição é abrupta entre o 5º e o 7º decil. Abaixo de **sim_top5 ≈ 0,20** o
modelo está no azar: quatro dos cinco primeiros decis têm AUROC indistinguível
de 0,5 e dois estão abaixo. Acima, 0,65–0,71.

**A janela produtiva de novidade é estreita: `0,20 <= sim_top5 < 0,25`, com
similaridade máxima < 0,25.** Fora dela, para cima é análogo; para baixo é
ruído. O `sim_T` deixa de ser filtro binário e passa a ser **banda**, com os
dois limites registrados no manifesto.

## Correção 2 — descontar o viés de tamanho: não funciona, e me refuta

O viés existe e é grande. Medindo a inclinação de cada uma das 616 cabeças
contra o número de bits ligados do ECFP4:

- **440 de 616 cabeças** com ρ > +0,3; **274** com ρ > +0,5
- mediana de ρ: **+0,471**
- o tamanho explica **18,8%** da variância do logit na cabeça mediana (p90: 37,7%)

Mas descontá-lo (resíduo do logit após regressão sobre o número de bits)
**piora**:

| escore | AUROC | IC95% |
|---|---|---|
| cru | **0,6322** | [0,6182, 0,6462] |
| corrigido | **0,5637** | [0,5498, 0,5776] |

Melhorou em só **23 de 88** alvos. Por classe, mediana por alvo: enzima 0,650 →
0,599; membrana 0,625 → 0,577; **epigenética 0,698 → 0,390**.

**Por que eu estava errado.** Tamanho não é artefato puro de contagem de bits.
Ligante maior faz mais contatos e tende genuinamente a ser mais potente — o
tamanho está **confundido com atividade real**, e removê-lo remove sinal
verdadeiro. O "defeito 2" do Adendo 7 estava meio errado.

O que está errado no caso dos antibióticos é outra coisa: **extrapolação**. Um
glicopeptídeo de 770 Da está fora da faixa de tamanho em que a relação
tamanho–atividade foi aprendida. E note a ironia: `n_bits` funciona como proxy
de confiança (delta +0,1265, z +7,14) na direção *mais bits, mais confiável* —
que é exatamente por que os antibióticos não foram sinalizados por ele.

**Correção certa: não descontar tamanho, e marcar tamanho FORA DE FAIXA como
fora do domínio** — verificação de dois lados contra a faixa de tamanho dos
ativos de treino do alvo, não desconto monotônico.

## Correção 3 — quem são os 19 alvos que acertam

Primeiro, uma invalidação própria: os preditores óbvios — AUROC, EF1% e BEDROC
no estrato análogo — **são retrospectivos**. Medem-se no mesmo conjunto de teste
pós-2022, então a correlação alta deles com a AUROC da novidade (ρ +0,482 para a
AUROC análoga) não serve de regra de decisão. Ficaram fora.

Entre os prospectivos, nulos: nº de referências ≤ 2020 (ρ +0,002 com acertar),
diversidade das referências (+0,031), tamanho do estrato (+0,015).

O que funciona é a **forma da distribuição de escore do modelo sobre a própria
biblioteca**, disponível antes de qualquer ensaio:

| preditor | ρ vs acertar@20 | ρ vs AUROC |
|---|---|---|
| razão p90/p50 da coluna | **+0,264** | +0,199 |
| mediana da coluna | **−0,204** | −0,233 |
| p90 da coluna | −0,154 | −0,216 |
| máximo da coluna | −0,032 | −0,140 |

Escore alto em tudo é ruim — assinatura de saturação. Opinião concentrada é bom.
Regra composta nos dois, cortes na mediana:

| | alvos | acertam top-20 | prec@20 | ativos novos (mediana) |
|---|---|---|---|---|
| **regra satisfeita** | 38 | **13 (34%)** | 3,9% | 10 |
| só razão alta | 6 | 0 (0%) | 0,0% | 12 |
| regra falha | 44 | 6 (14%) | 2,3% | 10 |

**Fisher exato, uma cauda: p = 0,0295.** E a mediana de ativos novos é 10 nos
dois grupos — não é o efeito de oportunidade, é seleção de alvo.

Por classe de proteína (n pequeno, registrado como indício):

| classe | alvos | acertam | AUROC mediana |
|---|---|---|---|
| enzyme | 48 | 11 (23%) | 0,650 |
| membrane | 21 | 4 (19%) | 0,625 |
| epigenetic | 7 | 3 (43%) | 0,698 |
| transcription | 6 | **0 (0%)** | **0,454** |
| transporter | 3 | 1 (33%) | 0,612 |

Fator de transcrição no azar é mecanisticamente coerente: sem bolsão definido.

### Correção à leitura do Adendo 6

O Adendo 6 disse "o retorno está concentrado em ~1/5 dos alvos", o que sugere
modelo errático. Estratificando pelo número de ativos novos **disponíveis** no
alvo:

| ativos novos | alvos | acertam top-20 | AUROC mediana |
|---|---|---|---|
| 3–4 | 11 | 1 (9%) | 0,642 |
| 5–9 | 32 | 3 (9%) | 0,576 |
| 10–19 | 23 | 7 (30%) | 0,676 |
| 20+ | 22 | 8 (36%) | 0,678 |

O "19 de 88" é em boa parte **oportunidade, não habilidade**: com 3 ativos num
universo de ~5.500, não acertar no top-20 é o esperado mesmo para um modelo bom.
No agregado, 6,4 acertos esperados por azar contra 50 observados — **7,8×**.

## O que entra no pipeline

1. **Trocar `p_max` pela média das 5 maiores similaridades** como proxy de
   domínio, e transformar a novidade em **banda `0,20 <= sim_top5 < 0,25`**, com
   os dois limites no manifesto.
2. **Não descontar tamanho.** Adicionar verificação de faixa de tamanho em dois
   lados contra os ativos de treino do alvo; fora da faixa, sai do ranque.
3. **Selecionar alvo antes de ensaiar** pela razão p90/p50 e pela mediana da
   coluna, e excluir fator de transcrição.

Arquivos: `corrigir.py`, `resultado_proxy_dominio.tsv`,
`resultado_vies_tamanho.tsv` (as 616 inclinações por cabeça),
`resultado_decis_similaridade.tsv`, `resultado_preditores_alvo.tsv`,
`classes_alvo.json`.
