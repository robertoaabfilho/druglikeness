# Etapa 1 — Coletar evidências

Desenho da primeira fase da reestruturação do pipeline `druglikeness`.
Data: 2026-11-25.

## Intenção e porquê

O pipeline hoje responde à pergunta errada. Ele ranqueia compostos pela
**semelhança de propriedades físico-químicas** a um perfil médio de fármacos
conhecidos para a doença — mede "parece um fármaco", não "age neste alvo". Foi
por isso que ciprofloxacino apareceu no top-10 de uma busca por analgesia: ele
se parece com analgésicos em peso molecular, TPSA e ADMET, sem nenhuma relação
farmacológica.

A reestruturação troca o eixo: **o alvo molecular passa a comandar**. Esta
etapa 1 constrói a base factual sobre a qual as fases seguintes (filtrar,
ranquear, validar) vão operar.

A regra que dá coesão ao conjunto, e que esta etapa estabelece: **cada fase só
acrescenta linhas e colunas, nunca sobrescreve, e toda afirmação carrega sua
procedência**. Uma previsão de modelo e um ensaio de bancada nunca se misturam
na mesma célula.

O objetivo mais amplo é um repositório público no GitHub, sem impedimento de
licença, reprodutível o bastante para ser citado na seção de métodos de um
artigo.

## Objetivos

1. Receber três entradas: sintomas da doença, alvos moleculares (do usuário **e**
   enriquecidos pelo pipeline) e ligantes no padrão do prospector.
2. Enriquecer a lista de alvos a partir dos sintomas, via Open Targets.
3. Produzir `evidencias.tsv` com **uma linha por evidência** de par
   ligante/alvo — medida ou prevista — rastreável até o ensaio e o documento.
4. Oferecer um painel Streamlit onde o usuário sobe os três arquivos e vê o
   resultado, sem precisar de linha de comando.
5. Gravar um manifesto por execução (versões das bases, checksum do checkpoint,
   parâmetros, sementes) que torne o resultado reproduzível e citável.

## Não-objetivos

- Ranquear, filtrar ou pontuar. Isso é fase 2 e 3.
- Decidir a biologia da doença. A lista de alvos é entrada respeitada como
  dada; o enriquecimento **acrescenta**, nunca substitui nem remove.
- Reescrever `ranking_multibase.py`. Ele continua funcionando como está até a
  fase 3.
- Substituir o Deep-PK nesta etapa (fase 2 faz isso, com ADMET-AI).

## Anti-objetivos

Coisas que tornariam esta etapa um fracasso mesmo funcionando:

- **Misturar medido com predito.** Se uma coluna puder conter tanto
  `pchembl=7,2` quanto `p=0,93`, o desenho falhou. Uma previsão não vale o
  mesmo que um ensaio.
- **Perder o rastro.** Uma linha de evidência que não diga de qual ensaio,
  documento e base veio não serve para publicação.
- **Casar alvo pelo nome sem organismo.** `noradrenaline transporter` existe em
  humano (CHEMBL222), camundongo (CHEMBL2370) e rato (CHEMBL304); juntar os
  três em silêncio invalida o resultado.
- **Dar confiança igual a alvos desiguais.** Alvos não humanos têm mediana de
  26 moléculas ativas no treino contra 67 dos humanos, e 798 têm menos de 10.
  Uma probabilidade de 0,9 sobre um alvo de 8 exemplos precisa ser distinguível
  de 0,9 sobre 4.000.

## Restrições

**Licença.** Repositório público exige cadeia limpa:

| componente | licença | consequência |
|---|---|---|
| PP4CKD (código) | Apache-2.0 | ok, atribuir |
| PP4CKD (pesos) | CC-BY-4.0 | ok, atribuir |
| RDKit | BSD-3 | ok |
| Open Targets | CC0 | ok |
| ChEMBL (dados) | CC BY-SA 3.0 | **share-alike**: arquivo derivado versionado no repositório carrega a mesma licença |
| Reactome (anotação) | CC0 | ok |
| Reactome (Functional Interaction) | **restrito** | **proibido**: não pode ser redistribuído sem licença da Pathway Solutions, e dependia de conteúdo KEGG |
| openFDA, PubChem | domínio público | ok |

O repositório **não tem LICENSE hoje** — sem ele, no GitHub o código fica
tecnicamente sem permissão de reuso. Decisão do usuário: **MIT**, a mesma do
ADMET-AI que entra na fase 2. O `LICENSE` vai na raiz nesta etapa.

Dado derivado do ChEMBL que seja versionado vai para `dados/` com `LICENSE`
próprio apontando CC BY-SA 3.0 — o share-alike alcança os dados, não o código.

**Reprodutibilidade.** O `pp4ckd_models.tar.gz` do Zenodo **está truncado na
origem**: o MD5 confere com o publicado, mas o tar termina cedo e só 22 dos ~77
modelos saem íntegros — `ecfp4_full_model.pt` não é um deles. O repositório
precisa documentar isso e o contorno (`ecfp4_fold_4.pt`), ou um terceiro não
reproduz.

**Formato de arquivo.** 42,8% dos nomes de produto natural na base real contêm
vírgula. Entrada e saída são **TSV**; `.csv` é aceito na leitura, com aviso
quando a detecção indicar quebra por vírgula.

**Escala.** Decisão do usuário: a etapa 1 precisa aguentar a **base inteira**,
os 108.935 compostos do prospector. Isso é arquitetura, não parâmetro, e
impõe três coisas:

- a evidência medida vem do `evidencias_chembl.csv` já baixado, não da API do
  ChEMBL ao vivo (a API não aguenta 100 mil ligantes); a API fica como opção
  para listas curtas;
- a escrita das saídas é em fluxo, não em memória — 100 mil ligantes × ~44
  alvos chega a milhões de linhas previstas;
- o painel executa em segundo plano, com progresso, em vez de um botão que
  trava.

**Dois repositórios separados.** O `prospector` não é importado. Este pipeline
conhece o *formato* dos dados dele, não o código.

## Abordagens consideradas

**A — Acrescentar etapas ao desenho atual.** Manter `main.py` + `config.txt` +
um script por etapa; entram duas etapas novas e as quatro desligadas voltam.
Rápido e com pouco churn. *Descartada*: deixa dois rankings concorrentes sem
hierarquia, mantém o monólito de 1.245 linhas do `ranking_multibase.py`, e não
registra o que gerou cada número — exatamente as duas coisas pedidas (coesão e
publicável) ficam de fora.

**B — Reorganizar em torno de uma tabela de evidências.** Quatro fases com
contrato explícito (coletar → filtrar → ranquear → validar), cada uma só
acrescentando colunas, com manifesto por execução. O `perfil` antigo deixa de
comandar e vira coluna de baseline, que é o que a comparação de um artigo
precisa ter ao lado. *Escolhida pelo usuário*, a ser entregue em fatias
verificáveis.

**C — Motor de workflow (Snakemake/Nextflow).** Padrão em artigo de
bioinformática, reprodutibilidade forte. *Descartada*: desproporcional para um
pipeline de um usuário, acrescenta uma DSL e um motor inteiro, e a
reprodutibilidade que compra se obtém com o manifesto de execução da opção B.

## Decisões tomadas pelo usuário

1. Abordagem **B**, trabalhada etapa por etapa, começando pela 1.
2. `evidencias` tem **uma linha por evidência**, não por par. Um par com três
   ensaios medidos e uma previsão vira quatro linhas. Foi a escolha explícita,
   contra a alternativa de uma linha por par com colunas paralelas, porque
   preserva o rastro até o documento e não força fundir origens.
3. Painel interativo em **Streamlit**, com upload dos três arquivos.
4. Colunas mínimas de `evidencias`: `ligante`, `alvo`, fonte da evidência,
   atividade relatada, mais o que embase a evidência.
5. **Corte do enriquecimento: colapso de família explícito.** Corte de escore
   e depois um representante por família, com os demais registrados numa
   coluna — controle auditável, em vez do critério implícito que eu havia
   recomendado (exigir presença no modelo).
6. **Escala: base inteira**, ~100 mil ligantes.
7. **Licença: MIT.**
8. **Reactome entra como fonte obrigatória** de expansão por via, ao lado do
   Open Targets. KEGG foi descartado por licença (API acadêmica apenas, e
   "prestar serviço com KEGG" exige licença paga via Pathway Solutions —
   exatamente a limitação que o repositório público não pode ter), e GEO foi
   adiado para a fase 2, como portão de expressão e não como fonte de alvos.

## Entradas

### `ligantes.tsv`

| coluna | obrigatória | para quê |
|---|---|---|
| `nome` | sim | rótulo legível |
| `smiles` | sim | entrada do modelo |
| `inchikey_parental` | recomendada | chave de junção; faz "cloridrato de X" e "X" serem o mesmo |
| `inchikey` | opcional | identidade exata, preserva sal e carga |
| `chembl_id` | opcional | liga à evidência medida |
| `fonte_smiles`, `conferencia` | opcional | procedência, herdadas do prospector |

Mínimo: `nome` + `smiles`. O que faltar, o pipeline recalcula — a química já
existe em `mods/importar_prospector.py` (Cleanup → maior fragmento orgânico →
neutralização).

### `alvos.tsv`

| coluna | obrigatória | para quê |
|---|---|---|
| `alvo_chembl_id` | sim, **ou** `alvo`+`organismo` | chave sem ambiguidade |
| `alvo` | alternativa | nome preferencial, casamento textual |
| `organismo` | sim, se usar `alvo` | evita ortólogo silencioso |
| `gene` | opcional | símbolo que o biólogo escreve (SCN9A) |
| `observacao` | opcional | por que o alvo está na lista; vai inteiro para a saída |

Mínimo: uma coluna, `alvo_chembl_id`.

### `sintomas.txt`

Um sintoma por linha, em inglês (é o idioma das bases). `#` comenta.

## Saídas

### `alvos_enriquecidos.tsv`

A lista do usuário **mais** o que o Open Targets associou aos sintomas.

`alvo_chembl_id`, `alvo`, `organismo`, `tipo_alvo`, `gene`, `origem`
(`usuario` | `open_targets` | `reactome`), `efo_id`, `sintoma`,
`escore_associacao`, `familia_colapsados`, `via_origem`, `alvo_semente`,
`n_moleculas_treino`, `no_modelo` (se existe no espaço do PP4CKD),
`observacao`.

### `evidencias.tsv` — uma linha por evidência

Identidade:
`ligante`, `ligante_inchikey_parental`, `ligante_chembl_id`,
`alvo_chembl_id`, `alvo`, `organismo`, `tipo_alvo`, `origem_alvo`.

Evidência:
`tipo_evidencia` (`medido` | `predito`), `fonte` (`ChEMBL` | `PP4CKD`),
`atividade_relatada` (texto legível: `IC50 = 120 nM`, `p = 0,93 (ECFP4)`),
`medida`, `relacao`, `valor`, `unidade`, `pchembl`.

Lastro:
`tipo_ensaio` (B/F/A), `ensaio_descricao`, `triagem_hts` (booleano),
`confianca_ensaio`, `n_moleculas_treino`, `documento`, `pmid`, `ano`,
`procedencia` (carimbo completo numa string).

`triagem_hts` é coluna de primeira classe porque foi medida: em ensaio
dirigido o recall@10 é 0,272, em triagem HTS cai para 0,111. Somar os dois sem
distinguir estraga qualquer leitura.

### `manifesto.json`

Por execução: data, versão do pipeline (commit), versão de cada base
consultada (release do Open Targets, versão do ChEMBL), nome e MD5 do
checkpoint do PP4CKD, parâmetros, sementes, e contagem de linhas de cada saída.
É o que permite citar o repositório e reproduzir.

## Painel Streamlit

```
┌──────────────────────────────────────────────────────────┐
│  druglikeness · etapa 1 — coletar evidências             │
├──────────────────────────────────────────────────────────┤
│  1. Sintomas        [ sintomas.txt        ] [Procurar]   │
│     neuropathic pain, peripheral neuropathy              │
│                                                          │
│  2. Alvos           [ alvos.tsv           ] [Procurar]   │
│     9 alvos · 9 resolvidos no modelo · 0 sem organismo   │
│     ⚠ 1 alvo com menos de 50 moléculas de treino         │
│                                                          │
│  3. Ligantes        [ ligantes.tsv        ] [Procurar]   │
│     200 ligantes · 200 SMILES válidos · 100 com ChEMBL   │
│                                                          │
│  [ Enriquecer alvos ]   → +37 alvos do Open Targets      │
│  [ Coletar evidências ]                                  │
├──────────────────────────────────────────────────────────┤
│  evidencias.tsv · 4.812 linhas                           │
│   medido 1.204 · predito 3.608 · triagem HTS 311         │
│  [tabela com filtro por tipo, organismo, fonte]          │
│  [⤓ evidencias.tsv] [⤓ alvos_enriquecidos.tsv]           │
│  [⤓ manifesto.json]                                      │
└──────────────────────────────────────────────────────────┘
```

O painel valida na hora do upload e **mostra o que vai dar errado antes de
rodar**: alvo sem organismo, SMILES ilegível, alvo ausente do modelo, alvo com
pouco treino. O botão de enriquecer é separado do de coletar para o usuário
poder ver e editar a lista de alvos antes de gastar a coleta.

O painel não é a interface única: tudo continua acessível por linha de comando,
e é a linha de comando que o artigo cita.



## Expansão por via (Reactome)

Fonte obrigatória, ao lado do Open Targets. Serve para achar o alvo que o
pesquisador não pensou em listar e que a associação doença→alvo não pega,
porque está a um passo na mesma via.

**Caminho:** `UniProt` → `GET /ContentService/data/mapping/UniProt/{acc}/pathways?species=9606`
→ para cada via, `GET /data/participants/{stId}/referenceEntities` → UniProt
dos participantes → mesmo mapeamento para ChEMBL da decisão 1.

**O problema, medido.** Expansão sem regra explode pior que o espalhamento
por família. Com quatro alvos semente (SLC6A2, SLC6A4, OPRM1, TRPV1):

| regra | proteínas únicas |
|---|---|
| sem teto | 573 |
| teto de 150 proteínas por via | 263 |
| teto de 80 proteínas por via | 90 |

Um único alvo, OPRM1, puxa 759 proteínas em 6 vias — 317 só de "G alpha (i)
signalling events". Com ~30 alvos semente, mesmo o teto de 150 passaria de
mil alvos.

**Regra adotada**, dois filtros em série porque o teto sozinho não basta:

1. **Teto de 150 proteínas por via.** Derruba os concentradores genéricos
   ("G alpha (i) signalling events" com 317, "Peptide ligand-binding
   receptors" com 200) e preserva os específicos ("Opioid Signalling" 90,
   "TRP channels" 28, "Serotonin clearance from the synaptic cleft" 3).
2. **O alvo expandido precisa existir no espaço do PP4CKD.** Sem isso o
   preditor não pontua e a entrada só acrescenta ruído. Quanto esse segundo
   filtro corta ainda **não foi medido** — medir é a primeira coisa a fazer
   na tarefa que implementa isto, e o número vai para o README.

**Auditabilidade.** Todo alvo vindo daqui sai marcado com `origem = reactome`,
`via_origem` (o `R-HSA-…` que o trouxe) e `alvo_semente` (qual alvo da lista o
puxou). Nenhuma adição é anônima, e dá para remover a expansão inteira
filtrando uma coluna.

**Restrição de licença.** Só dados de **anotação**, que são CC0. Os arquivos
de *Functional Interaction* do Reactome são proibidos aqui: não podem ser
redistribuídos sem licença da Pathway Solutions e dependiam de conteúdo KEGG
licenciado. Sem essa regra, o problema de licença do KEGG volta pela porta
dos fundos.

## Decisões que eu tomei (confira-as)

Um revisor independente leu este desenho e levantou dez pontos. Três foram
para o usuário (corte, escala, licença — acima). Os outros sete resolvi, e
ficam registrados aqui para serem conferidos, não enterrados.

**1. Como um gene vira `alvo_chembl_id`.** Verificado na prática, não suposto:

```
associatedTargets → target.proteinIds (UniProt)
                  → ChEMBL /target?target_components__accession=<up>&organism=...
                  → target_chembl_id
```

Testado em quatro alvos: P23975→CHEMBL222, P31645→CHEMBL228, P35372→CHEMBL233,
Q8NER1→CHEMBL4794 — exatamente os IDs curados à mão para dor neuropática. O
ChEMBL devolve mais de um alvo por UniProt (vem também `SELECTIVITY GROUP`);
a regra é **preferir `SINGLE PROTEIN`**.

**2. O que é "família", para o colapso.** A classificação de proteína do
próprio ChEMBL (`protein_class`), não uma definição nossa. Todo o resto do
pipeline já está ancorado em ChEMBL, então a regra fica citável. Os paralogos
colapsados vão para a coluna `familia_colapsados`, nunca sumem.

**3. Alvo enriquecido que não mapeia** fica na saída com `alvo_chembl_id`
vazio e `no_modelo = false`. Nunca é descartado em silêncio: o enriquecimento
só acrescenta.

**4. Ligante sem `chembl_id`.** Resolução por InChIKey exato e, falhando,
pelo bloco de conectividade (os 14 primeiros caracteres), preferindo a forma
parental. O ChEMBL não indexa InChIKey parental — a chave dele é por forma
salina —, então o bloco de conectividade é o que reconhece sal e base livre
como o mesmo composto. O revisor notou que o critério de aceite 2 não passaria
sem isso: a fixture `exemplos/controle.dor_neuropatica.tsv` não tem
`chembl_id`. Ela ganha a coluna.

**5. Quais registros do ChEMBL viram linha.** Todos, rotulados — é o
contrato "acrescenta e marca, não filtra". Com uma correção que o revisor
apontou com razão: `valor` + `unidade` **nunca são comparáveis entre
unidades** (120 nM e 78% não se ordenam), e por isso `pchembl` é a única
coluna comparável entre linhas, ficando vazia quando não existe. Registro com
relação `>` entra, com a relação marcada. Sem piso de confiança na coleta;
`confianca_ensaio` é carregada para a fase 2 decidir.

**6. Quais pares ganham linha `predito`.** Todo par (ligante × alvo que
exista no modelo) acima de `--corte`, cujo **padrão é 0,1** e vai no
manifesto, porque na base inteira o corte 0 produziria milhões de linhas de
ruído. Alvo fora do modelo não ganha linha prevista — nem uma linha com zero,
que seria afirmar algo que o modelo não disse.

**7. Quem resolve alvo fora do modelo.** A API de alvo do ChEMBL é o
resolvedor canônico; o `target_labels.tsv` só responde "está no modelo?".
Alvo fora do modelo é cidadão pleno para evidência **medida**.

**8. Ponto de entrada.** CLI novo e autônomo, `mods/coletar_evidencias.py`. A
ligação ao `main.py` acontece quando a fase 3 chegar — antes disso não há o
que orquestrar. Saídas em `results/` (ignorado pelo git); só um arquivo
pequeno de evidência do conjunto-controle é versionado, em `dados/`, para os
testes rodarem sem rede.

**9. Manifesto é saída e, opcionalmente, entrada.** Sempre gravado. Com
`--manifesto <arquivo>`, a execução entra em modo de repetição: recusa rodar
se a versão do Open Targets, a do ChEMBL ou o MD5 do checkpoint divergirem, e
serve as respostas de rede do cache. Sem isso, o critério de aceite 4 não se
sustenta — o Open Targets não garante a mesma resposta entre releases.

**10. Sintoma → EFO não é escolhido em silêncio.** O `ranking_multibase.py`
hoje pega `hits[0]` sem avisar, e um EFO errado produz uma lista de alvos
errada sem nenhum sinal. O painel mostra o EFO casado com as alternativas e
deixa trocar; a CLI aceita `--efo` para fixar.

**11. Pesos e rótulos no painel.** Caminho na barra lateral, com variável de
ambiente como padrão. O `target_labels.tsv` é obrigatório (o painel precisa
dele já no upload, para dizer "9 resolvidos no modelo"). Os pesos são
opcionais: sem eles o painel roda em **modo degradado**, só evidência medida,
e diz isso na tela em vez de falhar.

## O que fica a cargo de quem implementa

- Biblioteca HTTP e política de repetição nas chamadas a Open Targets e ChEMBL
  (o repositório já usa `requests`; manter).
- Formato do cache das respostas de rede. Que ele exista é requisito (ver
  decisão 9); o formato é livre.
- Como paginar as consultas ao Open Targets e à API de alvo do ChEMBL.
- Layout interno dos módulos: um por fonte de evidência ou um só com funções
  separadas.
- Separador decimal nas colunas de máquina: ponto. A vírgula fica confinada ao
  texto de exibição em `atividade_relatada`.
- Ordenação determinística das linhas de saída (necessária para o critério de
  aceite 4).
- Reaproveitar a heurística de `PUBCHEM_BIOASSAY`/`qHTS` já existente em
  `mods/validar_alvos_medidos.py` para `triagem_hts`, e o `confidence_score`
  do ensaio do ChEMBL para `confianca_ensaio`.

## Como se verifica que a etapa 1 funcionou

1. `sintomas.txt` com "neuropathic pain" devolve alvos do Open Targets que
   incluam os canais de cálcio α2δ e os transportadores de monoamina.
2. O conjunto-controle `exemplos/controle.dor_neuropatica.tsv` produz
   evidências em que duloxetina e amitriptilina têm linhas `medido` para o
   transportador de noradrenalina.
3. Nenhuma linha mistura `pchembl` com probabilidade de modelo.
4. Rodar duas vezes com o mesmo manifesto produz o mesmo arquivo.
