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

## Adendo 2026-11-26 — o canal negativo DEDUZIDO, medido

Pergunta que motivou a medicao: o usuario declara so o que E alvo, e o
pipeline deduz o antialvo. Da? Medido sobre o controle de 12, ECFP4,
`ecfp4_fold_4.pt`, varrendo os 7.676 alvos do modelo.

### A varredura nao custa previsao a mais

O forward do PP4CKD devolve os 7.676 de uma vez (`torch.sigmoid(rede(x))`). A
varredura e o MESMO calculo que pontuar os alvos pedidos, so nao descartando o
resto do vetor. A pergunta aberta "onde o canal negativo e calculado, ja que a
varredura custa uma previsao completa" partia de premissa errada. O que custaria
e GUARDAR a matriz (108.935 x 7.676 ~ 3,3 GB em float32), e nao e preciso: o
agregado por composto e O(1).

Dos 7.667 alvos nao pedidos, o requisito 1 descarta 4.122 por terem menos de 50
moleculas de treino. Sobram 3.545.

### Como deduzir decide o resultado

| canal negativo deduzido | AUC (trata > causa) |
|---|---|
| nenhum -- so os 9 alvos pedidos (linha de base) | 0,750 |
| **nao pedido + NAO HUMANO** (1.169 alvos) | **0,806** |
| nao pedido + humano (2.376 alvos) | 0,694 -- PIORA |
| contagem de promiscuidade (quantos acima de 0,5) | 0,583 -- PIORA |
| painel curado de seguranca (hERG, CYP3A4/2D6/2C9, MAO-A/B) | 0,750 -- nada |
| `max` sobre todos os 3.545 nao pedidos | 0,778, saturado |

Quatro desenhos eliminados por medicao:

1. **O lado nao humano deduzido empata com o declarado.** 0,806 e exatamente o
   numero que `p_humano x (1 - p_patogeno)` deu com o alvo do patogeno
   DECLARADO. Deduzir chega no mesmo lugar -- a deducao e viavel.
2. **O lado humano deduzido esta errado neste conjunto.** Sozinho e invertido
   da AUC 0,389, abaixo de 0,5: atividade alta em alvo humano nao pedido esta
   associada a TRATAR. Motivo com nome: amitriptilina, 0,910. As drogas que
   tratam dor neuropatica sao as sujas de SNC. Punir isso repete o erro do
   PAINS na etapa 2 -- apaga o positivo de referencia. Efeito fora do alvo em
   humano e assunto do portao (toxicidade/ADMET), nao do ranqueamento.
3. **Contar promiscuidade e o pior desenho:** AUC 0,583 e ZERA a amitriptilina
   (33 alvos nao pedidos acima de 0,5).
4. **Painel de seguranca nao ordena:** AUC 0,750, identica a linha de base;
   sozinho, 0,500 (ruido). Quem lidera o painel e a amitriptilina (0,526) --
   verdade farmacologica conhecida e irrelevante para esta pergunta.

E `max` sobre todos os nao pedidos SATURA: 0,91 a 0,97 para quase todo o
conjunto. Maximo sobre milhares de probabilidades e quase constante; o 0,778 e
artefato de misturar os dois lados, nao sinal.

### A regra que sobrevive a medicao

    canal negativo = alvo NAO pedido
                   & organismo != Homo sapiens
                   & n_moleculas_treino >= 50          (requisito 1)
                   & nao e ortologo de um alvo pedido  (requisito 2)

Ortologo = mesmo `pref_name` de um alvo pedido. Custa 9 alvos e nao muda a AUC
(0,806 com e sem): rede de seguranca de graca que cobre o contraexemplo do
CHEMBL304 (NET de rato) do requisito 2.

Isto NAO e a regra "nao humano = negativo" que este spec derrubou: o filtro de
PROCEDENCIA vem primeiro. No caso do *C. jejuni*, as enzimas que o usuario
listou sao pedidas, logo saem do canal antes do filtro de organismo. A
procedencia continua decidindo o sinal; o organismo so afia.

### Ressalva de validade, que e seria

n = 12 e os SEIS causadores sao todos antibioticos. A atividade nao humana os
separa por construcao do conjunto. O 0,806 pode ser propriedade deste controle,
nao regra geral. Logo: peso DECLARADO e desligavel, nunca constante ajustada
nesses 12, e nenhuma afirmacao de generalidade antes de um conjunto
independente.

### Decisao de interface que isto fecha

Nao entra coluna de sinal na entrada de alvos. A ideia de `sinal_alvo`
(`afirmativo`/`antialvo`) foi considerada e dispensada: o antialvo e deduzido.
Fica registrado, porem, o que a deducao NAO alcanca -- o antialvo que o usuario
conhece e que e HUMANO (as sialiltransferases do hospedeiro no caso do
mimetismo do *C. jejuni*), porque o lado humano do canal esta medido como
prejudicial. Se esse caso voltar, a coluna e a saida, e a medicao acima e o
motivo de ela ter sido dispensada por ora.

E um eixo que a entrada NUNCA deve declarar: silenciamento vs modulacao. A
evidencia do pipeline nao o carrega. Medido no `evidencias_controle.tsv`, a
coluna `medida` traz EC50 18, Ki 15, IC50 10, AC50 5, Potency 1 -- tudo
colapsado em `pchembl`, a unica coluna comparavel, e Ki e afinidade de ligacao,
neutra por definicao. O canal previsto do PP4CKD e P(ativo), sem direcao
alguma, e e o unico canal que existe para todo par sem ensaio no ChEMBL. Uma
coluna de direcao seria registrada e ignorada pela aritmetica.
