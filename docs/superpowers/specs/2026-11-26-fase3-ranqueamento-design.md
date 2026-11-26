# Fase 3 — Ranquear usando os alvos (desenho)

Requisitos: `2026-11-25-fase3-alvos-no-ranqueamento-requisitos.md`, incluindo o
adendo de 2026-11-26 que mediu o canal negativo deduzido. Este documento parte
de lá e decide a aritmetica.

## Intencao

As etapas 1 e 2 entregam evidencia e vereditos com procedencia, e **nada
ordena**. O ranqueamento que existe hoje (`ranking_multibase.py`) e cego a
alvo: `organismo`, `Homo`, `9606` e `specie` tem zero ocorrencias nele, no
`comparar_metodos_ranking.py` e no `ranking_embeddings.py`. Ele ordena por
semelhanca fisico-quimica com as drogas de referencia do sintoma.

O preco disso esta medido no controle de 12 farmacos: os dois QED mais altos
sao ciprofloxacino (0,893) e linezolida (0,892), **ambos do grupo que causa**
dor neuropatica, enquanto os que a tratam ficam entre 0,539 e 0,814. A
drogabilidade classica ordenou os causadores na frente. Nenhum ajuste de peso
corrige isso, porque a informacao que separa os dois grupos -- em que alvo o
composto age -- nao entra na conta.

Esta fase entrega `mods/ranquear.py`: uma ordem de compostos que usa os alvos.

## Objetivos

1. Ordenar compostos usando a evidencia de alvo da etapa 1, com o canal
   negativo deduzido do adendo.
2. Nunca misturar evidencia medida com prevista num mesmo numero.
3. Nao ordenar quem nao foi avaliado -- e dizer quem e.
4. Zero constante aprendida ou arbitrada por mim.
5. Repetir: manifesto, MD5 de entrada, versoes, saidas publicadas juntas.

## Nao objetivos

* Aprender pesos. n = 12 detecta sinal, nao calibra constante (requisitos).
* Faixa permeavel (previsao forte subindo de faixa): decidida como extensao a
  medir depois, com conjunto independente.
* Remover `ranking_multibase.py`, `farmaco_completo.py` e
  `checar_alertas_adversos.py`. Eles saem quando esta fase estiver validada e
  nada mais os chamar; ate la ficam marcados como substituidos no README.
* Qualquer regra derivada de classificacao de doenca (`MONDO_0005550` e
  parentes). Morreu no Guillain-Barre, nos requisitos.

## Anti-objetivos (fariam disto um fracasso mesmo funcionando)

* Um composto **sem evidencia nenhuma** receber um numero e aparecer ordenado
  como se tivesse sido avaliado. E o que acontece hoje.
* Uma constante minha no codigo decidindo a ordem.
* Punir o composto promiscuo: medido, AUC 0,583 e a amitriptilina zerada.
* `indeterminado` virando penalidade. Ausencia de dado nao e evidencia.

## Restricoes

* `pchembl` e a UNICA coluna comparavel entre linhas de `evidencias.tsv`.
  `valor` + `unidade` nao se comparam entre unidades, e `valor` de linha
  prevista e probabilidade, nao potencia.
* O canal negativo precisa do PP4CKD (checkpoint + `target_labels.tsv`), com o
  guarda de largura que a etapa 1 ja aplica.
* A base real tem 108.935 compostos: nada que exija a matriz completa
  (108.935 x 7.676 ~ 3,3 GB em float32) em memoria ou em disco.
* `MINIMO_TREINO = 50` ja existe em `mods/evidencias/resumo.py`. Reusar, nao
  reinventar.

## Abordagens consideradas

**A. Maximo e multiplicacao.** `max(p, pedidos) x (1 - max(p, negativo))`. Uma
linha, e e a forma que mediu AUC 0,806. Descartada como estrutura: `max`
ignora quantos alvos o composto acerta, e sobretudo ela **mistura medido com
previsto** num numero quando os dois existem para o par.

**B. Classe de evidencia primeiro, valor depois (lexicografica).** Faixas por
classe de evidencia; o valor ordena so dentro da faixa. E a regra da etapa 1
levada para a aritmetica em vez de abandonada nela, e e honesta sobre quem nao
tem evidencia. **Escolhida.**

**C. Soma ponderada declarada.** `Sigma wi.pi - Sigma vj.pj`. Expressa que um
alvo importa mais que outro, mas os pesos seriam meu arbitrio escrito no
codigo, ja que n = 12 nao os calibra. Trocar medicao por opiniao numerada.
Descartada.

Escolha do parceiro humano: **B por fora, A por dentro** -- as faixas dao a
estrutura, e dentro da faixa prevista vale o produto que mediu 0,806. E, entre
estrita e permeavel, **estrita**.

## O escore: uma chave lexicografica, sem constante inventada

Por composto, a ordem e dada pela chave

    (classe de evidencia, procedencia do alvo, valor dentro da faixa)

em ordem decrescente de forca:

| faixa | quem cai nela | valor que ordena DENTRO dela |
|---|---|---|
| **1a** | evidencia MEDIDA em alvo de `origem_alvo = usuario` | `max pchembl` sobre esses alvos |
| **1b** | evidencia MEDIDA em alvo do enriquecimento | `max pchembl` sobre esses alvos |
| **2a** | so PREVISTA, em alvo do usuario | `max p x (1 - max p_negativo)` |
| **2b** | so PREVISTA, em alvo do enriquecimento | `max p x (1 - max p_negativo)` |
| **3** | nenhuma evidencia em alvo afirmativo | **nao ranqueia** |

Tres decisoes dentro disto, cada uma com o seu motivo:

**Por que o canal negativo nao multiplica na faixa 1.** `pchembl` e medido;
`p_negativo` e previsto. Multiplica-los poria medido e previsto no mesmo
numero -- o que a etapa 1 proibe. E o 0,806 foi medido inteiramente sobre
valores previstos: e um numero da faixa 2, e aplica-lo a `pchembl` seria
estender a outra grandeza uma regra medida noutra. Na faixa 1 o
`p_negativo` sai como COLUNA, visivel e nao ordenante.

**Por que a procedencia e chave interna e nao peso.** O requisito 2 diz que
pedir e declarar: o alvo do usuario vale mais que o do enriquecimento. Dizer
isso com um peso exigiria uma constante; dizer com sub-faixa nao exige
nenhuma. Alvo do enriquecimento nunca cai no canal negativo porque o Open
Targets so tem alvos humanos, e o canal negativo e so nao humano -- os dois
lados sao consistentes sem regra extra.

**Consequencia da faixa estrita:** uma evidencia medida FRACA num alvo do
enriquecimento (1b) vem antes de uma previsao de 0,99 no alvo que o usuario
pediu (2a). Decorre de "medido sempre antes de previsto".

Apresentada ao parceiro humano como possivel defeito, e mantida por decisao
dele, com a razao dele: *"significa que o enriquecimento pode encontrar um alvo
fora do meu radar"*. E o argumento mais forte do que o conservadorismo que eu
havia dado: evidencia MEDIDA num alvo que o pesquisador nao tinha listado e
exatamente uma descoberta, e rebaixa-la sob a previsao do alvo ja conhecido
esconderia o unico resultado que o enriquecimento tem para dar. Por isso a
classe de evidencia e a chave EXTERNA e a procedencia a interna, e nao o
contrario -- a alternativa (procedencia por fora) foi considerada e descartada
aqui.

A extensao permeavel, se vier, e medida com conjunto independente.

`max` e nao soma para agregar varios alvos. **Correcao de 2026-11-26, medida ao
implementar a tarefa 5:** a versao anterior deste paragrafo dizia que soma
premiaria o composto promiscuo e que `max` era "o que foi medido". Nao se
sustenta -- no controle de 12, soma sobre os 9 alvos da a MESMA AUC de 0,750
que `max`. O conjunto nao discrimina as duas, entao `max` fica como escolha
DECLARADA, com o motivo que vale sem medicao: somar probabilidades de alvos
diferentes soma grandezas que nao se somam, e um composto que acerta cinco
alvos fracos ultrapassaria um que acerta o alvo certo com forca.

O que a medicao realmente rejeita e a CONTAGEM de promiscuidade (AUC 0,583, e a
amitriptilina zerada), nao a soma. `tests/test_ranque_auc.py` fixa os tres
numeros, inclusive este, para que a distincao nao se perca.

**Afirmativo** = `origem_alvo = usuario` (sub-faixa a) ou uma das
`ORIGENS_DE_ENRIQUECIMENTO` de `entrada.py:204`, hoje `open_targets` e
`reactome` (sub-faixa b). Nao ha terceira sub-faixa: as duas origens de
enriquecimento entram juntas, porque nada medido as distingue.

### Linha medida sem `pchembl`

Caso presente no unico dado de controle que existe: **1 das 49 linhas medidas
do `evidencias_controle.tsv` tem `pchembl` vazio** (ensaio cujo resultado nao
converte, por exemplo percentual de inibicao). Ela nao e comparavel, logo nao
pode ordenar. Decisao:

* ela **nao** entra no valor da faixa 1;
* o composto fica na faixa 1 se alguma OUTRA linha medida sua, em alvo
  afirmativo, tiver `pchembl`;
* se nenhuma tiver, o composto e ordenado pela previsao (faixa 2) e a medida
  nao comparavel sai como coluna (`medida_sem_pchembl`), visivel.

As duas alternativas foram rejeitadas: converter a medida numa escala
comparavel seria inventar grandeza, e rebaixar o composto para a faixa 3 diria
"sem evidencia" sobre um composto que FOI medido.

### O mesmo checkpoint nos dois lados da multiplicacao

`p_afirmativo` vem de `evidencias.tsv`, publicado pela etapa 1 com a
procedencia dela. `p_negativo` vem do forward desta fase. Se os dois sairem de
checkpoints diferentes, a multiplicacao `p x (1 - p_negativo)` combina dois
modelos diferentes num numero -- e a ordem passa a depender de qual par de
arquivos estava em disco.

Logo: a execucao **recusa rodar** (codigo 2) se o MD5 do checkpoint nao for o
que o manifesto da etapa 1 registra, nomeando os dois. E o mesmo guarda que a
etapa 1 ja aplica ao par checkpoint/rotulos, aplicado ao par etapa 1/fase 3.

### Sem os pesos do PP4CKD

A etapa 1 degrada para so evidencia medida; esta fase degrada tambem: sem
modelo nao ha canal negativo, a faixa 2 passa a ordenar por `max p` sozinho
(o `p` que a etapa 1 ja publicou), e o manifesto registra
`canal_negativo: false` com o motivo. O que nao pode acontecer e o canal
negativo valer zero em silencio -- zero e "nenhuma atividade indesejada", que
e afirmacao, nao ausencia.

## O canal negativo

    negativo = alvo NAO afirmativo
             & organismo != "Homo sapiens"
             & n_moleculas_treino >= MINIMO_TREINO (50)
             & pref_name != pref_name de qualquer alvo afirmativo

Agregado: `max` sobre esse conjunto. Calculado **no mesmo forward** que
pontua os alvos afirmativos -- o PP4CKD devolve os 7.676 de uma vez, entao a
varredura nao custa previsao a mais. Nada alem do agregado por composto e
guardado.

Medido no controle de 12 (adendo): 1.169 alvos no canal, AUC 0,806 contra
0,750 da linha de base. A exclusao de ortologo custa 9 alvos, nao muda a AUC e
cobre o CHEMBL304 (NET de rato) do requisito 2.

**Desligavel por flag**, e isto nao e luxo: n = 12 e os seis causadores do
controle sao todos antibioticos, logo a atividade nao humana os separa por
construcao do conjunto. `--sem-canal-negativo` tem de reproduzir 0,750.

## Entradas e saidas

Entradas: `evidencias.tsv` e `alvos_enriquecidos.tsv` (etapa 1), `portao.tsv` e
`propriedades.tsv` (etapa 2), checkpoint e `target_labels.tsv` do PP4CKD.

Saidas, num diretorio, **publicadas juntas**:

| arquivo | uma linha por | o que traz |
|---|---|---|
| `ranking.tsv` | composto ranqueado | posicao, faixa, valor, os componentes (`p_afirmativo`, `alvo_afirmativo`, `p_negativo`, `alvo_negativo`), o veredito do portao carregado, procedencia |
| `nao_ranqueados.tsv` | composto fora da ordem | o motivo: faixa 3, ou reprovado no portao, com o criterio |
| `manifesto.json` | — | versoes, MD5, parametros, dispositivo |

## O portao e o ranqueamento

`portao.tsv` traz `aprovado`. Decisao: **reprovado no portao nao ranqueia, mas
sai em `nao_ranqueados.tsv` com o criterio e o valor que decidiram** -- nada
desaparece em silencio. Composto apenas MARCADO ranqueia normalmente, com as
marcas como colunas: a etapa 2 marca justamente porque a medicao mostrou que
reprovar apaga os positivos de referencia, e o ranqueamento nao pode desfazer
isso por conta propria.

## Verificacao

Alem dos unitarios de faixa e de canal:

* o controle de 12 reproduz **0,750** com `--sem-canal-negativo` e **0,806**
  com o canal -- os dois numeros do adendo, como teste;
* um composto sem evidencia em alvo afirmativo sai em `nao_ranqueados.tsv`,
  nunca com posicao;
* ortologo de alvo pedido nao entra no canal negativo (CHEMBL304);
* repeticao por manifesto byte a byte;
* nenhum teste usa rede.

Codigos de saida: os das etapas 1 e 2 -- `0` completa, `3` parcial, `1` nenhum
composto utilizavel, `2` erro (nada publicado).

## Deixado ao implementador

* Nome e forma exata das flags, seguindo `aplicar_portao.py`.
* Como paginar/lotear o forward do PP4CKD sobre 108.935 compostos (a etapa 2 ja
  tem o padrao de lote em `aplicar_portao.py`).
* Formato da coluna que nomeia o alvo responsavel por `p_afirmativo` e por
  `p_negativo` (ID + nome + organismo numa celula, ou colunas separadas).
* Empate dentro de faixa: ordem estavel por `ligante_inchikey_parental`, para
  que duas execucoes deem o mesmo arquivo.
