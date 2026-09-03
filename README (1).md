# Drogabilidade — Pipeline de Análise e Ranking de Compostos

Pipeline em duas etapas para (1) calcular propriedades de drogabilidade/ADMET de uma lista de compostos e (2) ranqueá-los por semelhança ao perfil químico de drogas já usadas para tratar uma doença específica.

## Visão geral

```
input.txt ──► farmaco_completo.py ──► farmaco.csv ──► ranking_multibase.py ──► ranking.csv
```

1. **`farmaco_completo.py`** — recebe uma lista de SMILES, calcula descritores de drogabilidade (RDKit) e consulta a API ADMET do Deep-PK, salvando tudo combinado em `farmaco.csv`.
2. **`ranking_multibase.py`** — recebe o nome de uma doença, busca drogas já usadas para tratá-la em 4 bases públicas, deriva pesos por característica a partir delas, e usa esses pesos para ranquear os compostos do `farmaco.csv`, por um de dois métodos (perfil ou TOPSIS).

Os dois scripts precisam estar na **mesma pasta** — o `ranking_multibase.py` importa funções diretamente do `farmaco_completo.py` (reaproveita o cálculo de descritores em vez de duplicá-lo).

---

## Requisitos

```bash
pip install rdkit requests numpy
pip install psycopg2-binary   # opcional, só para a fonte DrugCentral do ranking_multibase.py
```

---

## 1. `farmaco_completo.py`

Une o cálculo de descritores de drogabilidade com a consulta ADMET em um único fluxo.

### Formato do `input.txt`

Um composto por linha, aceitando qualquer um destes formatos (pode misturar entre linhas):

```
CC(=O)Oc1ccccc1C(=O)O
Aspirina,CC(=O)Oc1ccccc1C(=O)O
Aspirina;CC(=O)Oc1ccccc1C(=O)O
```

Linhas vazias ou iniciadas com `#` são ignoradas.

### Uso

```bash
python farmaco_completo.py --input input.txt --output farmaco.csv
```

| Argumento | Descrição | Padrão |
|---|---|---|
| `--input` | Arquivo de entrada com os SMILES | `input.txt` |
| `--output` | CSV final de saída | `farmaco.csv` |
| `--pred-type` | Categoria de predição do Deep-PK (`admet`, `absorption`, `distribution`, `metabolism`, `excretion`) | `admet` |
| `--batch-size` | Quantos SMILES enviar por job ao Deep-PK | `5` |
| `--checkpoint` | Arquivo de progresso, para retomar após interrupção | `deeppk_progress.json` |

### O que ele calcula

- **Descritores locais (RDKit):** peso molecular, doadores/aceptores de ligação H, TPSA, ligações rotacionáveis, refratividade molar, átomos pesados.
- **Regras de drogabilidade:** Regra dos 5 de Lipinski (parcial, sem LogP) e Regra de Veber.
- **ADMET (Deep-PK):** propriedades gerais, absorção, distribuição, metabolismo, excreção e toxicidade.

### Resiliência

- Salva progresso incrementalmente (`farmaco.deeppk_raw.csv` + `deeppk_progress.json`) — pode ser interrompido e retomado.
- Se uma molécula específica for rejeitada pela API do Deep-PK, o script **isola e pula só ela** (fica sem dados ADMET, mas não trava o resto do processamento nem fica preso em loop).

---

## 2. `ranking_multibase.py`

Ranqueia os compostos de um `farmaco.csv` por semelhança ao perfil químico de drogas já usadas para tratar uma doença.

### Uso

```bash
python ranking_multibase.py --doenca "Guillain-Barre syndrome" --farmaco farmaco.csv --output ranking.csv
```

Vários termos podem ser passados juntos (cada um é buscado separadamente e depois tudo é mesclado):

```bash
python ranking_multibase.py --doenca "Guillain-Barre syndrome" "Miller Fisher syndrome" --farmaco farmaco.csv --output ranking.csv
```

> **As buscas nas 4 bases são em inglês** e por correspondência literal de texto — siglas em português (ex: "SGB") ou termos que não aparecem no texto/nome em inglês da condição não retornam nada. Use sempre o nome da doença em inglês.

| Argumento | Descrição | Padrão |
|---|---|---|
| `--doenca` | Palavra(s)-chave da doença | *obrigatório* |
| `--farmaco` | CSV de entrada (gerado pelo `farmaco_completo.py`) | `farmaco.csv` |
| `--output` | CSV de ranking de saída | `ranking.csv` |
| `--metodo` | `perfil` ou `topsis` (ver seção abaixo) | `perfil` |
| `--top-por-fonte` | Quantas drogas manter de cada fonte, priorizando maior relevância | `10` |
| `--no-opentargets` / `--no-chembl` / `--no-drugcentral` / `--no-openfda` | Desativa uma fonte específica | — |
| `--no-admet` | Não consulta ADMET das drogas de referência (mais rápido) | — |

### Como funciona

**Etapa 1 — Buscar drogas de referência em 4 bases (gratuitas, sem chave):**

| Fonte | O que traz | Critério de corte (top 10) |
|---|---|---|
| **Open Targets** | Drogas conhecidas + fase clínica | Maior fase clínica (Aprovado > Fase IV > ... > Pré-clínico) |
| **ChEMBL** | Indicações por ID EFO/MONDO ou texto | Maior fase clínica da indicação |
| **DrugCentral** | Indicações já aprovadas (via Postgres público) | Correspondência exata ao nome da doença, depois nº de produtos aprovados |
| **openFDA** | Menções em bula (`indications_and_usage`) | Nº de rótulos de bula diferentes que mencionam a doença |

Cada fonte é consultada **por inteiro primeiro** (sem cortar durante a busca) e só depois reduzida às 10 mais relevantes.

**Etapa 2 — Mesclar as 4 listas:** deduplicação por **ID ChEMBL** (mais confiável que nome em texto), com fallback para nome normalizado nas fontes que não fornecem ID (DrugCentral, openFDA).

**Etapa 3 — Garantir o SMILES** de cada droga de referência: usa o já trazido (DrugCentral) → senão ChEMBL pelo ID → senão PubChem pelo nome (pedindo `IsomericSMILES`, `CanonicalSMILES` e `ConnectivitySMILES` juntos, já que o PubChem às vezes devolve só uma dessas chaves).

**Etapa 4 — Calcular os descritores** dessas drogas via `processar_lista_druglikeness()`, importada do `farmaco_completo.py` (mesma lógica, mesmas colunas do `farmaco.csv`).

**Etapa 5 — (opcional) ADMET das drogas de referência** via Deep-PK, preenchendo as mesmas colunas ADMET do `farmaco.csv`.

**Etapa 6 — Calcular os pesos** por característica, de acordo com a ocorrência entre as drogas de referência:

- **Numéricas** (PesoMolecular, TPSA...): peso maior quanto mais **consistente** o valor entre as drogas — `peso = 1 / (1 + desvio/média)`.
- **Booleanas** (Passa_Ro5, Passa_Veber): peso = frequência do valor mais comum (moda).
- **Categóricas** (ADMET): extrai a categoria do texto (ex: `"Yes (High confidence)"` → `"Yes"`) e usa a frequência da categoria mais comum.

**Etapa 7 — Pontuar e ranquear**, por um dos dois métodos (`--metodo`):

#### Método `perfil` (padrão)

Para cada composto, soma `peso × similaridade` por característica (similaridade numérica via kernel gaussiano em torno da média de referência; booleana/categórica via match sim/não), normaliza para 0–100.

#### Método `topsis`

Implementação do **TOPSIS** (Technique for Order of Preference by Similarity to Ideal Solution) híbrida:

1. Monta uma matriz de decisão real (compostos × critérios).
2. Injeta o **perfil das drogas de referência** (média para numéricas, moda para booleanas/categóricas) como uma linha-alvo adicional.
3. Normaliza tudo pela **norma vetorial** de cada coluna (`r_ij = x_ij / √Σx_ij²`) — matriz de compostos e linha-alvo juntas, na mesma escala.
4. Aplica os pesos: `v_ij = peso_j × r_ij`.
5. A **Solução Ideal Positiva (A+)** é a própria linha-alvo (o perfil de referência) já normalizada e ponderada — não o máximo/mínimo do próprio `farmaco.csv`, por isso "híbrido".
6. A **Solução Ideal Negativa (A-)** é, para cada critério, o valor observado nos compostos mais distante do alvo.
7. Calcula a distância euclidiana até A+ (`D_mais`) e A- (`D_menos`), e o coeficiente de proximidade relativa `C = D_menos / (D_mais + D_menos)` — a `Pontuacao_Final` do ranking é `C × 100`.

O `ranking.csv` gerado por este método traz duas colunas extras: `D_mais` e `D_menos` (as distâncias euclidianas), além da `Pontuacao_Final`.

> `perfil` e `topsis` respondem perguntas ligeiramente diferentes: `perfil` é "o quanto esse composto se parece com a média das drogas de referência, característica por característica"; `topsis` é "o quanto esse composto está próximo do perfil de referência **em relação à variação observada entre os seus próprios candidatos**" — mais sensível à dispersão do seu próprio conjunto de compostos.

### Arquivos de saída

- **`ranking.csv`** — ranking final (`Ranking`, `Nome`, `SMILES`, `Pontuacao_Final`, e `D_mais`/`D_menos` se `--metodo topsis`).
- **`ranking.pesos.json`** — pesos calculados por característica (transparência de por que cada composto pontuou como pontuou).
- **`ranking.referencia.csv`** — drogas de referência usadas, com a coluna `Fontes` mostrando de onde cada uma veio (ex: `chembl+opentargets`).

---

## Dicas de busca por doença

- Use sempre o **nome em inglês** da condição (as 4 bases são todas em inglês).
- Para doenças autoimunes pós-infecciosas (ex: Síndrome de Guillain-Barré ligada a *Campylobacter jejuni*), busque pelo **nome da síndrome resultante**, não pelo patógeno — o tratamento indexado nas bases é o da condição autoimune (imunomoduladores), não antimicrobianos contra o agente causador.
- Se o interesse é nos **sintomas** (ex: dor neuropática) e não na doença de base em si, pode fazer mais sentido buscar pelo sintoma/condição geral (`"neuropathic pain"`, `"peripheral neuropathy"`) em vez do nome específico da doença rara — bases com mais dados geram pesos mais confiáveis.
- Pode combinar vários termos numa mesma chamada (`--doenca termo1 termo2 ...`), mas termos muito distintos entre si podem diluir a especificidade do perfil de referência.

---

## Limitações conhecidas

- **DrugCentral** depende de uma conexão direta ao Postgres público (não há API REST oficial estável); requer `psycopg2`. Se a tabela `product` não existir nesse dump, o script cai automaticamente para um critério de corte mais simples.
- **openFDA** não fornece SMILES nem ID ChEMBL — sempre usa o fallback do PubChem, e os nomes costumam vir como sal (ex: "DONEPEZIL HYDROCHLORIDE"), então pode não deduplicar com as outras fontes.
- Compostos biológicos (anticorpos), extratos vegetais/homeopáticos e sais inorgânicos simples encontrados como referência frequentemente não têm um SMILES único válido e são descartados com aviso — isso é esperado, não é um bug.
- A API do Deep-PK pode rejeitar moléculas específicas; o script pula essas moléculas automaticamente (ficam sem dados ADMET) em vez de travar.
- O método `topsis` não é "TOPSIS de manual" — é uma adaptação (ideal externo ao invés de derivado só do pool de alternativas). Se o objetivo for rigor metodológico formal para publicação, vale documentar essa escolha explicitamente.
