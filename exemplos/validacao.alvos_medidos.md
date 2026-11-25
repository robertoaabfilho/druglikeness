# recall@10 contra alvos medidos

Medida de `mods/prever_alvos.py` com `mods/validar_alvos_medidos.py`, sobre
13.914 compostos do pipeline de prospecção que têm bioatividade medida no
ChEMBL (amostra de 6.000 previstos). Checkpoint `ecfp4_fold_4.pt`, ECFP4.

Reproduzir:

```bash
python mods/validar_alvos_medidos.py \
    --evidencias ~/prospector/.../evidencias_chembl.csv \
    --compostos  ~/prospector/.../atividades_compostos_parcial.csv \
    --rotulos    ~/modelos_baixados/pp4ckd/target_labels.tsv \
    --checkpoint ~/modelos_baixados/pp4ckd/models/ecfp4_fold_4.pt
```

## O vazamento NÃO foi eliminado — todo número aqui é limite superior

Era o plano, e não deu. Dois motivos:

1. A divisão de folds (`data/cv_splits_10fold.npz`) **não foi publicada** — não
   está no repositório do PP4CKD nem no tarball do Zenodo. Então "excluir os
   compostos do fold 4" não é possível: a divisão não é conhecida nem
   reproduzível sem reconstruir a extração do ChEMBL 36 inteira.
2. O caminho alternativo era o limiar de treino: o PP4CKD só usou como rótulo
   positivo pares com `pchembl >= 5,0` (relação `=`, `<`, `<=`, `~` ou nula) ou
   inibição >= 50%. Pares com `pchembl < 5` nunca foram ensinados. Mas depois de
   remover as linhas de triagem HTS sobram **21 pares**, amostra pequena demais
   para estimar nada.

Os números abaixo portanto medem memorização e generalização juntas. O que
continua informativo, mesmo assim, é **como o desempenho varia** — e a variação
é monotônica, o que é a assinatura de uma medida real.

## Resultados

| estrato | pares | rec@1 | rec@10 | hit@10 |
|---|---:|---:|---:|---:|
| todos os pares com pchembl | 32.239 | 0,049 | 0,238 | 0,621 |
| **qualidade do ensaio** | | | | |
| triagem HTS/PubChem (ruído) | 6.727 | 0,003 | 0,111 | 0,255 |
| ensaio dirigido (não-HTS) | 25.512 | 0,061 | 0,272 | 0,675 |
| não-HTS, ligação (tipo B) | 8.425 | 0,088 | 0,342 | 0,585 |
| **treino do alvo** (não-HTS, pchembl≥6) | | | | |
| alvo com < 50 moléculas ativas | 2.178 | 0,033 | 0,127 | 0,454 |
| alvo com 50–500 | 4.254 | 0,049 | 0,147 | 0,523 |
| alvo com 500–5.000 | 4.452 | 0,067 | 0,303 | 0,563 |
| alvo com ≥ 5.000 | 2.364 | 0,088 | **0,477** | 0,661 |
| **organismo** (não-HTS, pchembl≥6) | | | | |
| alvo humano | 9.846 | 0,041 | 0,209 | 0,626 |
| alvo não humano | 3.402 | 0,113 | **0,387** | 0,654 |

Referência do artigo (ECFP4, fold 0, conjunto retido): rec@1 = 0,541,
rec@10 = 0,856.

## Três leituras

**O recall@10 cru de 0,24 não é comparável ao 0,856 do artigo**, e a diferença
é de denominador, não só de modelo. O regime do artigo tem ~1,3 rótulos por
composto (dedutível de `recall@1 = 0,541` com `precision@1 = 0,726`); aqui a
mediana é 2 e a cauda chega a **846** alvos medidos num único composto, por
painéis de triagem. Com mais de 10 rótulos, o recall@10 está limitado por
construção. Restringindo a compostos com ≤2 rótulos, sobe de 0,21 para 0,449.
O resto da diferença é esperado: estes são produtos naturais de
NPASS/COCONUT/LOTUS, estruturalmente mais distantes do núcleo medicinal do
ChEMBL do que as moléculas de uma divisão aleatória do próprio ChEMBL.

**O treino do alvo domina tudo**: 0,127 → 0,147 → 0,303 → 0,477 conforme o
alvo sai de <50 para ≥5.000 moléculas ativas. É um fator de 3,8. É por isso que
`prever_alvos.py` emite a coluna `moleculas_treino` e avisa abaixo de 50 —
sem ela, p = 0,9 num alvo de 8 ativos é indistinguível de p = 0,9 num alvo de
4.000.

**Alvo não humano vai melhor, não pior**: rec@10 0,387 contra 0,209 do humano.
Não é surpresa depois de inspecionar: boa parte dos alvos não humanos aqui são
ensaios fenotípicos de organismo inteiro com muitíssimo treino
(*M. tuberculosis* 11.049 ativos, *P. falciparum* 41.726) e quimiotipos
distintos — antibióticos e antiparasitários são separáveis. Isso sustenta a
decisão de tratar organismo não humano como caso de primeira classe.

## A medida que corresponde ao uso real

Recall sobre 7.676 alvos é a pergunta errada para este pipeline: você **fornece
uma lista curta de alvos**. A medida que não depende da densidade de rótulos é a
posição do alvo medido no ranking (ensaio dirigido, pchembl ≥ 6, 13.248 pares):

| percentil | posição do alvo medido | melhor alvo do composto |
|---:|---:|---:|
| 10 | 2 | 1 |
| 25 | 10 | 1 |
| **50** | **79** de 7.676 (top 1,0%) | **3** |
| 75 | 458 | 18 |
| 90 | 1.159 | 101 |

Traduzindo: numa lista fornecida de 10 alvos, o alvo verdadeiro ficaria entre os
3 primeiros em **94,9%** dos pares; com 20 alvos, 87,0%; com 50, 72,7%.

Essa tradução é uma aproximação otimista: ela supõe que os outros alvos da lista
se comportam como alvos sorteados. Uma lista inteira de alvos plausíveis **para
a mesma doença** é mais difícil, porque são alvos correlacionados e muitos serão
genuinamente prováveis. Trate 94,9% como teto, não como expectativa.
