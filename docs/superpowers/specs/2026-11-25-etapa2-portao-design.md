# Etapa 2 — O portão

Desenho da segunda fase da reestruturação do pipeline `druglikeness`.
Data: 2026-11-25. Etapa anterior: `2026-11-25-etapa1-coletar-evidencias-design.md`.

## Intenção e porquê

A etapa 1 produz `evidencias.tsv`: uma linha por evidência de par
ligante/alvo, medida ou prevista. A etapa 2 decide **o que segue para o
ranqueamento** e, principalmente, **o que fica marcado**.

O nome "portão" descreve a posição, não a severidade. A medição abaixo mostra
que quase todo critério plausível de filtragem, aplicado como exclusão,
apagaria os próprios positivos de referência do pipeline. Então o portão
**marca, anota e explica; por padrão quase nada reprova**. Quem quiser
estreitar o portão promove um critério a reprovação por linha de comando, e o
manifesto registra a escolha.

Isso também resolve três dívidas concretas:

1. **Deep-PK sai.** `mods/farmaco_completo.py` consulta uma API externa em
   lotes de 5, sequencialmente: ~21.800 requisições para a base inteira, com
   um bug que rotula falha de rede como "rejeitado pelo Deep-PK" (linhas
   583-585), avança o checkpoint de todo jeito (607) e o apaga no fim (613).
   Entra **ADMET-AI** (MIT, local, pesos no próprio pacote, sem rede).
2. **GEO**, adiado da etapa 1, é respondido aqui — com medição, abaixo.
3. **openFDA** volta como evidência negativa, com a regra corrigida: a atual
   erra nos dois sentidos (medido abaixo).

## Objetivos

1. Calcular ADMET localmente, sem rede, para a base inteira, em tempo viável.
2. Anotar cada ligante com alerta estrutural (PAINS/BRENK/NIH) **nomeado**.
3. Buscar na bula (openFDA) a condição-alvo como efeito adverso, separando
   seção forte de seção fraca e descontando a indicação.
4. Anotar cada alvo com sua expressão por tecido (GTEx/Tabula Sapiens).
5. Produzir um veredito por ligante (`aprovado`), com o motivo de cada
   critério e o valor que o decidiu — reversível por coluna.
6. Dizer, no resumo, quantos ligantes **cada critério reprovaria** se fosse
   promovido: o botão é medido antes de ser apertado.

## Não-objetivos

- Ranquear. É a fase 3.
- Decidir a biologia. O portão não escolhe alvo nem tecido; anota.
- Reescrever `ranking_multibase.py` ou `main.py`. A ligação é da fase 3.
- Treinar ou calibrar modelo ADMET. ADMET-AI entra como publicado.

## Anti-objetivos

- **Apagar em silêncio.** Um ligante reprovado continua na saída, com o
  critério e o valor que o reprovaram.
- **Reprovar por critério não medido.** Cada reprovação padrão tem um número
  no conjunto-controle abaixo.
- **Misturar predição com medida** (regra herdada da etapa 1). ADMET-AI é
  predição; bula é documento regulatório. Colunas e procedências distintas.
- **Tratar ausência como ausência de risco.** Ligante sem bula no openFDA é
  `indeterminado`, nunca "sem alerta".

## O que foi medido (e por que o portão quase não reprova)

Tudo abaixo foi medido nesta máquina, nos arquivos do repositório.

### Alerta estrutural

| catálogo | `input.exemplo.tsv` (200 produtos naturais) | controle (12 fármacos) |
|---|---|---|
| PAINS (A+B+C) | 30 (15%) | 2 (17%) |
| BRENK | 130 (65%) | 6 (50%) |
| NIH | 60 (30%) | 1 (8%) |

Quem o PAINS pega no controle: **amitriptilina** (trata) e linezolida (causa).
Quem o BRENK pega: **carbamazepina** e **capsaicina** (tratam), mais quatro.
Reprovar por PAINS apagaria a amitriptilina — o positivo que a etapa 1 usa
como critério de aceite. Reprovar por BRENK apagaria 65% da base.
**Conclusão: alerta estrutural marca, não reprova.**

### QED (drogabilidade clássica)

No controle, os dois QED mais altos são **ciprofloxacino (0,89)** e
**linezolida (0,89)** — ambos do grupo `causa`. Os que tratam ficam entre 0,54
e 0,81. QED não separa tratar de causar; é exatamente a premissa da
reestruturação, agora com número. **QED é coluna, nunca critério.**

### Expressão por tecido (a pergunta do GEO)

GTEx (bulk RNA-seq, TPM mediano) para os 9 alvos curados de dor neuropática:

| gene | maior TPM (tecido) | maior TPM em tecido nervoso |
|---|---|---|
| GRIN1 | 255,3 (córtex frontal) | 255,3 |
| CACNA2D2 | 39,2 (cerebelo) | 39,2 |
| CACNA2D1 | 37,6 (gastrocnêmio) | 16,7 |
| TRPV1 | 22,7 (cerebelo) | 22,7 |
| SCN1A | 12,0 (córtex pré-frontal) | 12,0 |
| OPRM1 | 5,9 (cerebelo) | 5,9 |
| **SLC6A4** (SERT) | 35,0 (pulmão) | **0,45** |
| **SLC6A2** (NET) | 5,5 (testículo) | **0,08** |

Um portão de "≥1 TPM em tecido nervoso" **reprovaria SERT e NET** — os dois
alvos da duloxetina e da amitriptilina, os mais bem validados da lista. A
causa é conhecida: RNA-seq de tecido inteiro dilui núcleos pequenos (locus
coeruleus, núcleos da rafe) no córtex amostrado. **Expressão marca, não
reprova.**

**GEO, decisão.** O caminho pelo GEO cru não entra. Dois motivos, o segundo
medido: (a) séries do GEO são experimentos individuais, sem vocabulário de
tecido harmonizado nem normalização comum — o trabalho de integração é maior
que a etapa inteira; (b) a via mantida e harmonizada para esse mesmo dado é o
Expression Atlas, que o Open Targets indexa como fonte `expression_atlas`, e
**para `EFO_0005762` (dor neuropática) ela devolve ZERO evidências** nos alvos
testados (SLC6A2, TRPV1, SCN9A). Não há transcriptômica humana de
caso-controle para essa doença — não se biopsia gânglio da raiz dorsal. Então
o que o portão anota é expressão **basal** (GTEx + Tabula Sapiens, CC0, via
`target.baselineExpression` do Open Targets), e a coluna de expressão
diferencial existe, vem do escore de tipo de dado `rna_expression` da
associação, e sai zerada nesta doença — registrado, não escondido.

### Alerta de bula (openFDA)

A regra atual (`mods/checar_alertas_adversos.py`) erra nos dois sentidos.
Medida no controle, buscando `neuropathy|neuropathic|neuritis|paresthesia` em
`boxed_warning`/`warnings`/`warnings_and_cautions`/`adverse_reactions` com
`limit=1`:

- **4 de 6 falsos positivos** no grupo `trata`: pregabalina, duloxetina,
  amitriptilina e carbamazepina são sinalizadas, porque a bula de um fármaco
  para dor neuropática lista parestesia e neuropatia entre as reações
  adversas;
- **falso negativo no caso que motivou o script**: ciprofloxacino **não** é
  sinalizado, porque `limit=1` devolve uma bula arbitrária entre centenas e a
  que traz a tarja preta de neuropatia periférica não é a primeira.

A regra nova, medida no mesmo controle:

1. `limit=50` bulas por princípio ativo, agregadas (qualquer bula com o termo
   conta);
2. **seção forte** (`boxed_warning`, `warnings`, `warnings_and_cautions`)
   conta como alerta; `adverse_reactions` é registrada e **não** conta;
3. **indicação desconta**: termo presente em `indications_and_usage` significa
   que o fármaco *trata* a condição, e cancela o alerta.

| grupo | alerta pela regra nova |
|---|---|
| `trata` (6) | **0** — pregabalina e duloxetina batem em seção forte, mas são indicadas; as outras quatro não batem |
| `causa` (6) | **4** — ciprofloxacino (tarja preta), linezolida, etambutol, metronidazol |

Os dois `causa` não pegos são honestos: cloranfenicol tem **1** bula no
openFDA (só seção fraca) e etionamida **nenhuma** — ambos saem
`indeterminado`, não "sem alerta". 4/6 de sensibilidade com 0/6 de falso
positivo, contra 4/6 de falso positivo e o caso-chave perdido na regra atual.

**Este é o único critério que reprova por padrão**, e por isso: tem
especificidade medida, é a queixa que originou o projeto, e o rastro inteiro
(seção, termo, trecho, número de bulas) fica em `alertas_bula.tsv`.

### Desempenho do ADMET-AI

ADMET-AI 2.0.1, MIT, pesos dentro do pacote (10 modelos Chemprop: 5 de
classificação + 5 de regressão), 104 colunas por molécula (52 valores + 52
percentis contra os aprovados do DrugBank). Nenhuma chamada de rede.

| medição | tempo |
|---|---|
| carga do modelo | 0,1 s (6,3 s de import) |
| 1.600 moléculas, `num_workers=8` (padrão com GPU) | 384 s |
| 1.600 moléculas, `num_workers=0` | **8,6 s** |

**45× de diferença**: com GPU, o `ADMETModel` escolhe `num_workers=8` e paga o
nascimento de 8 processos por modelo, dez vezes. O pipeline passa
`num_workers=0` por padrão, com `--workers` para quem quiser outro valor. Nos
números de `num_workers=0`: physchem 7 ms/mol, chemprop 5,4 ms/mol, percentis
3 ms/mol → a base inteira (108.935) sai em ~28 min, contra as ~21.800
requisições sequenciais do Deep-PK.

## Restrições

**Licença.** ADMET-AI: MIT (código e pesos no pacote). Conjunto de referência
do DrugBank embutido no pacote: usado apenas em memória, para calcular
percentil; **não é redistribuído por este repositório**. openFDA: domínio
público. Open Targets (incluindo `baselineExpression`): CC0. GTEx e Tabula
Sapiens chegam pelo Open Targets, não baixados daqui. KEGG e Reactome
Functional Interaction continuam proibidos.

**Escala.** ADMET-AI aguenta a base inteira (medido). A bula **não**: é uma
requisição por ligante a um serviço público com limite de taxa, e 108.935
consultas não são aceitáveis. Portanto a checagem de bula tem teto explícito
(`--limite-bula`, padrão 500) e o que passar do teto sai `indeterminado` com
esse motivo escrito. Cache em disco, como na etapa 1.

**Reprodutibilidade.** Manifesto próprio da etapa, no formato da etapa 1:
versões (ADMET-AI, RDKit, Open Targets, openFDA `last_updated`), MD5 das
entradas, parâmetros, política de critérios, contagens. Modo de repetição
(`--manifesto`) pelas mesmas regras.

**Degradação.** Sem `admet_ai` instalado a etapa roda e diz que rodou sem
ADMET: as colunas de ADMET-AI saem vazias e o critério `toxicidade` sai
`indeterminado`. Sem rede, `--sem-rede` pula bula e expressão. Nenhuma dessas
ausências vira aprovação ou reprovação silenciosa.

## Abordagens consideradas

**A — Penalizar na pontuação, como hoje.** `checar_alertas_adversos.py`
multiplica a pontuação por 0,5. *Descartada*: a pontuação é da fase 3, e
embutir a penalidade aqui faz o portão opinar sobre o ranking antes de ele
existir; além disso uma penalidade não é auditável por coluna.

**B — Portão que marca, com reprovação opcional por critério.** Veredito por
critério, `aprovado` derivado da política, motivo e valor em colunas, e o
resumo dizendo quanto cada critério custaria se promovido. *Escolhida*: é a
única forma compatível com as medições acima, onde quase todo critério
plausível apaga um positivo de referência.

**C — Portão duro com limiares de literatura.** Reprovar por PAINS, BRENK,
hERG > 0,5, etc. *Descartada por medição*: apagaria amitriptilina,
carbamazepina, capsaicina e 65% da base.

## Entradas

- `--ligantes` (obrigatório): o mesmo TSV da etapa 1 (`nome`, `smiles`, …),
  lido pelo mesmo `entrada.ler_ligantes`.
- `--evidencias` (opcional): `evidencias.tsv` da etapa 1. Quando dado,
  restringe o portão aos ligantes que têm evidência, e o resumo diz quantos
  ficaram de fora.
- `--alvos` (opcional): `alvos_enriquecidos.tsv` da etapa 1; liga a anotação
  de expressão.
- `--tecidos` (opcional): termos UBERON ou nomes de tecido de interesse, um
  por linha. Sem isso a expressão é anotada sem recorte.
- `--termos-alerta` (opcional): termos da condição a buscar na bula; sem eles
  a checagem de bula não roda.

## Saídas

Todas TSV, cabeçalho `#`, escape reversível, ponto decimal, escrita atômica,
publicadas juntas — as regras de formato da etapa 1, sem exceção.

### `propriedades.tsv` — uma linha por ligante

Identidade (`ligante`, `ligante_inchikey_parental`, `ligante_chembl_id`,
`smiles`), descritores do RDKit (peso molecular, doadores/aceptores de HB,
TPSA, ligações rotacionáveis, refratividade molar, átomos pesados, QED),
alertas estruturais **nomeados** (`alerta_pains`, `alerta_brenk`,
`alerta_nih`, cada um com os nomes dos filtros casados), e as **104 colunas do
ADMET-AI com os nomes originais** — é o que torna a coluna citável. Mais
`procedencia` (`ADMET-AI <versão>` / `RDKit <versão>`).

### `portao.tsv` — uma linha por ligante

`ligante`, `ligante_inchikey_parental`, `aprovado` (`sim`/`nao`), `motivos`
(critérios que reprovaram, `;`), e por critério: `<criterio>_veredito`
(`passa`/`reprova`/`marca`/`indeterminado`) e `<criterio>_valor` (o número ou
o texto que decidiu). `politica` registra, por linha, a política em vigor.

### `alertas_bula.tsv` — uma linha por achado

`ligante`, `termo`, `secao`, `forca` (`forte`/`fraca`/`indicacao`), `trecho`,
`n_bulas_consultadas`, `campo_do_nome`, `procedencia` (`openFDA
<last_updated>`). Uma linha por (termo, seção): o rastro é o produto.

### `alvos_anotados.tsv` — uma linha por alvo

As colunas de `alvos_enriquecidos.tsv` mais: `expressao_maxima`,
`tecido_maximo`, `expressao_nos_tecidos_pedidos`, `tecidos_pedidos_casados`,
`unidade_expressao`, `fonte_expressao` (`gtex` / `tabula_sapiens`),
`escore_rna_expression` (expressão diferencial da associação; zero quando a
doença não tem dado), `expressao_veredito`.

### `manifesto.json`

Como na etapa 1, mais: versão do ADMET-AI, `num_workers`, política de
critérios, teto de bula, contagens por veredito, e **quantos ligantes cada
critério reprovaria se promovido**.

## Critérios e política padrão

| critério | o que olha | padrão | por quê |
|---|---|---|---|
| `smiles` | RDKit lê o SMILES | **reprova** | sem molécula não há o que avaliar |
| `alerta_bula` | bula (regra nova) | **reprova** | 0/6 falso positivo medido |
| `alerta_estrutural` | PAINS/BRENK/NIH | marca | reprovaria amitriptilina e 65% da base |
| `toxicidade` | ADMET-AI (hERG, AMES, DILI, ClinTox, Carcinogens) acima do percentil `--percentil-toxicidade` (padrão 90) dos aprovados do DrugBank | marca | percentil é relativo a fármaco aprovado, mas não foi validado neste conjunto |
| `expressao` | expressão basal do alvo nos tecidos pedidos | marca | reprovaria SERT e NET |

`--reprovar C[,C]` promove; `--marcar C[,C]` rebaixa. Critério sem dado é
`indeterminado` e **nunca** reprova, qualquer que seja a política.

## Códigos de saída

Os da etapa 1: `0` completa, `3` parcial, `1` nenhuma entrada utilizável,
`2` erro (nada publicado).

## O que fica a cargo de quem implementa

- Tamanho do lote do ADMET-AI (a carga é fixa; o lote só afeta memória).
- Formato do cache da bula e da expressão (o da etapa 1 serve).
- Como casar nome de tecido pedido com UBERON. **Não usar casamento por
  trecho de nome**: na medição acima `cortex of kidney` casou com um padrão de
  "nervoso". Usar o identificador UBERON ou a ancestralidade do Biosample.
- Ordenação determinística das saídas (necessária para repetição).

## Como se verifica que a etapa 2 funcionou

1. O controle de 12 fármacos produz `alerta_bula` em ciprofloxacino,
   linezolida, etambutol e metronidazol, e em **nenhum** dos seis que tratam.
2. `input.exemplo.tsv` (200 ligantes) sai com 104 colunas de ADMET-AI
   preenchidas, sem rede, e 130 marcados por BRENK — nenhum reprovado por
   isso.
3. Com `--reprovar alerta_estrutural`, o resumo diz que a amitriptilina seria
   reprovada; sem a flag, ela é aprovada e marcada.
4. Sem `admet_ai` instalado, a etapa termina com código 0, as colunas de ADMET
   vazias e `toxicidade = indeterminado`.
5. `--manifesto` de uma execução anterior reproduz as mesmas saídas e recusa
   rodar se uma entrada, a versão do ADMET-AI ou a política mudar.

---

## Adendo de 2026-11-25 (durante a implementação): a janela de 50 bulas

Medições feitas ao implementar a tarefa da bula, que **mudaram a regra** acima.
Ficam aqui porque a regra original foi decidida por medição e só outra medição
a corrige.

### O defeito

O desenho tirou o `limit=1` porque *"devolve uma bula arbitrária entre
centenas"*. Com `limit=50` o problema não desapareceu, só mudou de escala:

- a **capsaicina** tem **211 bulas** no openFDA, e a **única** que a indica
  para dor neuropática (Qutenza) ficava fora das 50 devolvidas;
- o desconto por indicação falhava em silêncio, a capsaicina era reprovada, e
  ela **trata** a condição: o falso positivo do grupo de controle ia de
  **0/6 para 1/6**.

### O que a API permite (medido)

| caminho | o que vê | custo (ciprofloxacino, 166 bulas) |
|---|---|---|
| `limit=50` sem ordenar | 50 arbitrárias | 8,6 MB |
| `limit=50&sort=effective_time:desc` | as 50 vigentes | 8,6 MB |
| `limit=1000` (máximo da API) | todas | **27 MB** |
| consulta dirigida (`search` por seção e termo, `limit=1`) | **o acervo inteiro** | **202 KB** |

**"As 50 mais acessadas" não existe:** a bula do openFDA tem 34 campos e
nenhum de acesso, popularidade, visualização ou contagem de uso. A data
(`effective_time`) é o único critério de relevância disponível, e significa
"rotulagem vigente".

### A regra nova

1. **A decisão vem de duas consultas dirigidas** — seção forte e indicação —
   que rodam no acervo inteiro (265.009 bulas no export de 2026-11-25). A
   seção fraca não tem consulta própria porque não decide nada.
2. **O rastro vem da janela**, agora ordenada por `effective_time:desc`: o
   trecho citado em `alertas_bula.tsv` sai da rotulagem vigente.
3. **O peso da evidência passa a ser coluna** (`n_bulas_com_o_termo`): 118 das
   166 bulas do ciprofloxacino trazem `neuropathy` em seção forte, contra
   **1 de 211** na indicação da capsaicina. Uma contagem feita dentro da
   janela não distingue um alerta robusto de um fio de cabelo.

Custo: **três requisições por ligante** (janela, forte, indicação).

A ordenação por data, só por si, também salvaria a capsaicina (a Qutenza é de
2026-07-07 e entra nas 50 mais recentes). As duas coisas ficam porque a
ordenação resolveu aquele caso **por sorte de data** — a indicação estava em 1
bula de 211 — e uma bula antiga com genéricos novos por cima escaparia.
Recência correlaciona com "ainda comercializado", não com "carrega o alerta".

Verificado depois da mudança, contra o openFDA ao vivo: **0/6 de falso
positivo e 4/6 pegos**, os mesmos números, agora com o peso visível.

### Requisito para a base inteira (não implementado)

O teto `--limite-bula` (padrão 500) existe porque **nenhum** dos caminhos acima
escala para os 108.935 ligantes da base:

- 3 requisições por ligante = ~327 mil requisições;
- o openFDA permite **240 por minuto** sem chave (1.000 com chave) →
  **~23 horas** no melhor caso, e o serviço é público.

O caminho que escala é o **download em massa**: o openFDA publica o acervo de
bulas como JSON particionado — **1,7 GB em 14 partes, 265.009 registros**,
export de 2026-11-25, em `https://api.fda.gov/download.json`. Baixar uma vez e
indexar localmente por nome troca 327 mil requisições por um download e
elimina o limite de taxa. Fica como requisito de quem for rodar a base inteira;
até então o teto é a proteção, e o que passa dele sai `indeterminado` com esse
motivo escrito.

### Adendo 2: as seções que faltavam

A regra original olhava `boxed_warning`, `warnings` e `warnings_and_cautions`
como seção forte. Medido depois, uma consulta dirigida por seção no
conjunto-controle:

| seção | falso positivo (grupo que trata) | ganho (grupo que causa) |
|---|---|---|
| `contraindications` | **0/6** | só etambutol (8/8) — já pego |
| `precautions` | **0/6** | cipro (1) e etambutol (8) — já pegos |
| `warnings_and_precautions` | 0 | 0 — campo não usado |
| `drug_interactions` | 0 | 0 |
| `use_in_specific_populations` | pregabalina 227, duloxetina 113, capsaicina 1 | **0** |

**Entram como seção forte: `contraindications` e `precautions`.** Nenhuma das
duas muda veredito neste controle — e é justamente por isso que entrar é
seguro. A razão de entrarem não é sensibilidade:

- uma **contraindicação** que nomeia a condição é a afirmação mais forte que
  uma bula pode fazer ("não use neste paciente"). Omiti-la deixava o critério
  dependente de o fabricante ter repetido o aviso em outra seção. Nestes 12 ele
  repetiu; 12 fármacos não autorizam supor que repita sempre;
- **`precautions`** fecha uma **lacuna de formato**: bula de formato antigo
  (não-PLR) usa `precautions` onde a nova usa `warnings_and_cautions`. Sem ela,
  o pipeline não olhava o campo equivalente nas bulas antigas.

**`use_in_specific_populations` fica de fora, com o número que a excluiu:** três
acertos, todos no grupo que **trata**, zero no que causa. Um teste prende a
exclusão para que ninguém a adicione sem medir de novo.

Custo: **zero requisição a mais** — as seções entram no mesmo `OR` da consulta
dirigida que já existia. Verificado ao vivo depois da mudança: **0/6 de falso
positivo e 4/6 pegos**, e o rastro passa a citar `contraindications` e
`precautions` no etambutol, que antes só apareciam como `warnings`.

### Adendo 3: a GPU não reproduz, e por isso a CPU é o padrão

O critério de verificação 5 pede que `--manifesto` reproduza as mesmas saídas.
Na primeira repetição real, `propriedades.tsv` **diferiu** — e as outras três
saídas não. Medido:

| dispositivo | colunas diferentes entre duas execuções da mesma molécula |
|---|---|
| GPU | **36 de 312** |
| CPU | **0 de 312** |

A diferença nos valores é de ~1e-8 relativo (`23.619516` contra `23.619518`),
ordem de redução nos kernels da GPU. Mas o **percentil amplifica**: ~0,035 de
diferença, porque um percentil é uma posição entre os ~2.850 aprovados do
DrugBank e um fio de diferença inverte uma comparação.

Nenhum veredito mudou — `portao.tsv` saiu idêntico nas duas. O portão decide
por percentil ≥ 90, e 0,035 não atravessa esse limiar. Mas uma coluna que muda
sozinha entre execuções não é citável, e o modo de repetição existe exatamente
para isso.

Preço de exigir CPU, medido nas 200 moléculas de `input.exemplo.tsv`:

| dispositivo | por molécula | base de 108.935 |
|---|---|---|
| GPU | 20,5 ms | ~37 min |
| CPU | 36,4 ms | ~66 min |

**1,8×, não 10×.** Por 29 minutos numa execução que precisa ser citada, a CPU
passa a ser o padrão. `--gpu` troca reprodutibilidade por tempo, avisa na hora
que faz isso, e:

- o manifesto registra o dispositivo **real**, não o pedido (`dispositivo()`
  responde pelo que vai valer; um manifesto que registrasse o pedido mentiria
  quando o pedido não pegasse);
- `dispositivo` é **chave da repetição**: repetir uma execução de GPU na CPU é
  recusado;
- repetir uma execução feita na GPU **avisa** que as colunas de ADMET podem
  diferir na última casa e que isso não é falha da repetição — em vez de
  prometer bit-exatidão que a GPU não entrega.

Verificado depois da mudança: `propriedades.tsv`, `portao.tsv` e
`alertas_bula.tsv` saem byte a byte idênticos na repetição.
