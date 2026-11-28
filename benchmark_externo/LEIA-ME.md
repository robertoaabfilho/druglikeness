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
