# Drogabilidade — Pipeline de Análise e Ranking de Compostos

Pipeline em duas etapas para (1) calcular propriedades de drogabilidade/ADMET de uma lista de compostos e (2) ranqueá-los por semelhança ao perfil químico de drogas já usadas para tratar uma doença específica.

## Visão geral

```
input.txt ──► farmaco_completo.py ──► farmaco.csv ──► ranking_multibase.py ──► ranking.csv
```

1. **`farmaco_completo.py`** — recebe uma lista de SMILES, calcula descritores de drogabilidade (RDKit) e consulta a API ADMET do Deep-PK, salvando tudo combinado em `farmaco.csv`.
2. **`ranking_multibase.py`** — recebe o nome de uma doença, busca drogas já usadas para tratá-la em 4 bases públicas, deriva pesos por característica a partir delas, e usa esses pesos para ranquear os compostos do `farmaco.csv`.

---

## Requisitos

```bash
pip install rdkit requests
pip install psycopg2-binary   # opcional, só para a fonte DrugCentral do ranking_multibase.py
```

Os dois scripts precisam estar na **mesma pasta** — o `ranking_multibase.py` importa funções diretamente do `farmaco_completo.py` (reaproveita o cálculo de descritores em vez de duplicá-lo).

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
python ranking_multibase.py --doenca "Alzheimer's disease" --farmaco farmaco.csv --output ranking.csv
```

| Argumento | Descrição | Padrão |
|---|---|---|
| `--doenca` | Palavra(s)-chave da doença | *obrigatório* |
| `--farmaco` | CSV de entrada (gerado pelo `farmaco_completo.py`) | `farmaco.csv` |
| `--output` | CSV de ranking de saída | `ranking.csv` |
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

**Etapa 3 — Garantir o SMILES** de cada droga de referência: usa o já trazido (DrugCentral) → senão ChEMBL pelo ID → senão PubChem pelo nome.

**Etapa 4 — Calcular os descritores** dessas drogas via `processar_lista_druglikeness()`, importada do `farmaco_completo.py` (mesma lógica, mesmas colunas do `farmaco.csv`).

**Etapa 5 — (opcional) ADMET das drogas de referência** via Deep-PK, preenchendo as mesmas colunas ADMET do `farmaco.csv`.

**Etapa 6 — Calcular os pesos** por característica, de acordo com a ocorrência entre as drogas de referência:

- **Numéricas** (PesoMolecular, TPSA...): peso maior quanto mais **consistente** o valor entre as drogas — `peso = 1 / (1 + desvio/média)`.
- **Booleanas** (Passa_Ro5, Passa_Veber): peso = frequência do valor mais comum (moda).
- **Categóricas** (ADMET): extrai a categoria do texto (ex: `"Yes (High confidence)"` → `"Yes"`) e usa a frequência da categoria mais comum.

**Etapa 7 — Pontuar e ranquear:** para cada composto do `farmaco.csv`, soma `peso × similaridade` por característica, normaliza para 0–100, ordena → salva em `ranking.csv`.

### Arquivos de saída

- **`ranking.csv`** — ranking final (`Ranking`, `Nome`, `SMILES`, `Pontuacao_Final`).
- **`ranking.pesos.json`** — pesos calculados por característica (transparência de por que cada composto pontuou como pontuou).
- **`ranking.referencia.csv`** — drogas de referência usadas, com a coluna `Fontes` mostrando de onde cada uma veio (ex: `chembl+opentargets`).

---

## Limitações conhecidas

- **DrugCentral** depende de uma conexão direta ao Postgres público (não há API REST oficial estável); requer `psycopg2`. Se a tabela `product` não existir nesse dump, o script cai automaticamente para um critério de corte mais simples.
- **openFDA** não fornece SMILES nem ID ChEMBL — sempre usa o fallback do PubChem, e os nomes costumam vir como sal (ex: "DONEPEZIL HYDROCHLORIDE"), então pode não deduplicar com as outras fontes.
- Compostos biológicos (anticorpos), extratos vegetais/homeopáticos e sais inorgânicos simples encontrados como referência frequentemente não têm um SMILES único válido e são descartados com aviso — isso é esperado, não é um bug.
- A API do Deep-PK pode rejeitar moléculas específicas; o script pula essas moléculas automaticamente (ficam sem dados ADMET) em vez de travar.
