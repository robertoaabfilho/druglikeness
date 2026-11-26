# Fase 3 — requisitos sobre o uso dos alvos no ranqueamento

Nota de requisitos, nao desenho. Levantada em 2026-11-25, durante a execucao
da etapa 2, para nao ser redescoberta. O desenho da fase 3 parte daqui.

## O buraco que motiva a nota

As etapas de previsao e validacao tratam alvo nao humano como caso de primeira
classe; a etapa de ranqueamento de hoje e **cega a alvo**. Medido:
`organismo`, `Homo`, `9606` e `specie` nao aparecem em `ranking_multibase.py`,
`comparar_metodos_ranking.py` nem `ranking_embeddings.py` -- zero ocorrencias.
Essa etapa ranqueia composto por semelhanca fisico-quimica com as drogas de
referencia do sintoma; nao ha coluna de alvo, logo nao ha organismo.

Consequencia: o alvo e bem tratado em todo lugar onde existe como entidade, e
invisivel exatamente na etapa que produz o resultado final.

## Ganho medido que o pipeline NAO usa

Conjunto-controle de 12 farmacos (6 que tratam dor neuropatica, 6 que a
causam), ECFP4, `ecfp4_fold_4.pt`. Escore do composto = maior probabilidade
entre os alvos pedidos.

| escore | AUC |
|---|---|
| so os 9 alvos humanos | 0,750 |
| `p_humano x (1 - p_patogeno)` | **0,806** |
| so `CHEMBL360` (*M. tuberculosis*) | **1,000** |

Com o alvo bacteriano no calculo, o ciprofloxacino cai do 4o para o 7o lugar.
`p_bacteriano` nao aparece em nenhum `.py`: e medicao documentada no README e
em `exemplos/LEIA-ME.md`, e nada a consome.

**Ressalva de tamanho:** n = 12. Serve para detectar o sinal, nao para
calibrar peso. Constante ajustada nesses 12 esta sobreajustada; os pesos da
fase 3 sao declarados, nao aprendidos ali.

## Requisito 1 — alvo esparso e ausencia de evidencia, nao evidencia de inatividade

Medido: quase toda a perda de AUC do conjunto humano vem da gabapentina
(0,0001) e da pregabalina (0,02), e a causa e que a subunidade alfa-2-delta-2
tem **30 moleculas ativas no treino** -- sendo ela o alvo real das duas. Alvo
nao humano e mais esparso ainda: mediana de 26 moleculas ativas contra 67 dos
humanos, 67% abaixo de 50, 798 abaixo de 10.

Logo: probabilidade BAIXA em alvo pouco treinado nao entra como negativa; e
descartada. A coluna `n_moleculas_treino` ja existe em `evidencias.tsv` e e o
que decide. E a mesma regra do portao da etapa 2 (`indeterminado` nunca
reprova), aplicada a aritmetica.

## Requisito 2 — o sinal do alvo vem da PROCEDENCIA, nao do organismo

O organismo e necessario e insuficiente.

*Necessario*: sem ele, ortologos se confundem -- `noradrenaline transporter` e
`CHEMBL222` (humano), `CHEMBL2370` (camundongo) e `CHEMBL304` (rato). Por isso
a etapa 1 o exige por nome (`entrada.py:275`) e o resolve por ID, recusando o
alvo que fique sem ele (`medido.py:125`).

*Insuficiente*: ele nao carrega intencao. Dois contraexemplos dos dados deste
repositorio:

* `CHEMBL360` -- mesmo ID, mesmo organismo, sinal oposto: afirmativo numa
  busca antituberculose, negativo numa busca de analgesico.
* `CHEMBL304` (NET de rato) -- nao humano e **afirmativo**: e o ortologo do
  alvo pretendido. Uma regra "nao humano = negativo" o marcaria como atividade
  indesejada.

Entao:

* **afirmativo = o alvo que o usuario pediu.** A lista de alvos E a
  declaracao de intencao. `origem_alvo` (etapa 1) ja distingue alvo pedido de
  alvo trazido pelo enriquecimento. Pedir e declarar -- regra que o
  `prever_alvos.py` ja aplica ("pedir por ID vence o filtro de organismo").
* **negativo = atividade em alvo que ninguem pediu**, achada na varredura dos
  7.676. Esse canal nao consulta a doenca.
* **a doenca segue fazendo o que ja faz**: EFO -> alvos humanos associados,
  para o enriquecimento. Nada alem.

## Requisito 3 — alvo do usuario tem de ser ranqueavel sem escore de associacao

Caso que produz o requisito: combater o mimetismo molecular de *Campylobacter
jejuni*, listando as enzimas da via de biossintese do LOS da bacteria
(`cst-II`, `neuA/B/C`) como alvos a inibir.

O Open Targets so tem alvos humanos, logo essas enzimas **nunca** chegam pelo
enriquecimento. Se o escore exigir `escore_associacao`, esse caso inteiro
pontua zero. Portanto alvo de procedencia "usuario" entra no escore com peso
proprio, sem depender de associacao.

O lado humano do mimetismo (atividade do composto nas sialiltransferases do
hospedeiro) e preocupacao legitima e cai no canal negativo, porque ninguem
pediu esses alvos.

## Regra DESCARTADA, e por que

Foi considerada e medida contra a ontologia: *"doenca infecciosa abre o canal
afirmativo para alvo nao humano; o resto fecha"*, detectando a doenca
infecciosa por `MONDO_0005550` nos ancestrais ou nas areas terapeuticas do EFO
que o usuario ja fixa com `--efo`.

Consultado na API do Open Targets em 2026-11-25:

| doenca | classificacao |
|---|---|
| tuberculose (`MONDO_0018076`) | `MONDO_0005550` infectious disease |
| malaria (`MONDO_0005136`) | hematologic disorder + `MONDO_0005550` |
| dor neuropatica (`EFO_0005762`) | `phenotype`; nenhum ancestral infeccioso |
| **Guillain-Barre (`MONDO_0016218`)** | **nervous system + immune system disorder; SEM `MONDO_0005550`** |

O Guillain-Barre mata a regra: a sindrome e causada pelo mimetismo molecular
do LOS de *C. jejuni*, mas nao e classificada como infecciosa. A regra
fecharia o canal afirmativo e as enzimas bacterianas que o usuario listou para
inibir entrariam como atividade indesejada -- o oposto da intencao declarada.

Segundo motivo, independente: a malaria mostra que saber que a doenca e
infecciosa nao diz QUAL patogeno conta -- atividade em *M. tuberculosis*
continua fora do alvo numa busca antimalarica. E derivar o organismo do nome
da doenca e o padrao que o desenho da etapa 2 ja proibiu, por medicao (um
padrao de nome para "nervoso" casou `cortex of kidney`).

A procedencia resolve os dois casos sem dicionario e sem ontologia.

## Em aberto para o desenho da fase 3

* Como agregar varios alvos num escore por composto (maximo, soma ponderada,
  ou classe de evidencia primeiro e valor depois).
* Se medido e previsto podem entrar na mesma conta ou so na mesma ordenacao
  (a etapa 1 proibe misturar as colunas; a fase 3 decide a aritmetica).
* Onde o canal negativo e calculado: a varredura dos 7.676 alvos por composto
  custa uma previsao completa, e a base tem 108.935 compostos.
* Validar `MONDO_0005550` -- ou o que o substituir -- em mais que quatro
  doencas, se alguma regra voltar a depender de classificacao de doenca.
