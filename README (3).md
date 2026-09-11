# Drogabilidade — Pipeline de Análise e Ranking de Compostos

Pipeline para (1) calcular propriedades de drogabilidade/ADMET de uma lista de compostos e (2) ranqueá-los por semelhança ao perfil químico de drogas já usadas para tratar um sintoma/condição específica — com ferramentas de validação e um orquestrador que roda tudo a partir de um único arquivo de configuração.

## Estrutura de pastas

```
./main.py                  <- orquestrador (único arquivo solto na raiz)
./config.txt                <- você cria/edita (main.py gera um modelo na 1ª vez)
./input.txt                 <- você edita (lista de SMILES)
./mods/                     <- todos os scripts do pipeline
    configuracao.py
    farmaco_completo.py
    ranking_multibase.py
    comparar_metodos_ranking.py
    validar_enriquecimento.py
    ranking_embeddings.py
    calibrar_deeppk.py
./results/                  <- criada automaticamente; todas as saídas vão aqui
    farmaco.csv
    farmaco.deeppk_raw.csv
    deeppk_progress.json
    ranking.csv
    ranking.referencia.csv
    ranking.pesos.json
    comparacao_metodos.csv          (se rodar_comparacao_metodos = true)
    comparacao_metodos.spearman.csv
    validacao_enriquecimento.csv    (se rodar_validacao_enriquecimento = true)
    validacao_enriquecimento.cortes.csv
    ranking_embeddings.csv          (se rodar_embeddings = true)
```

## Visão geral do fluxo

```
input.txt ──► mods/farmaco_completo.py ──► results/farmaco.csv
                                                  │
config.txt (sintomas) ──► mods/ranking_multibase.py ──► results/ranking.csv
                                                  │        results/ranking.referencia.csv
                                                  │        results/ranking.pesos.json
                    ┌─────────────────┬───────────┘
                    ▼                 ▼                 ▼
     comparar_metodos_ranking.py  validar_enriquecimento.py  ranking_embeddings.py
        (validação: 4 métodos)   (validação: leave-one-out)  (ranking alternativo)
```

`mods/calibrar_deeppk.py` existe como ferramenta complementar, fora do fluxo do `main.py` (decisão do projeto: não seguir com essa abordagem — ver seção 6).

---

## Requisitos

Só `rdkit` e `requests` são estritamente necessários para começar — **o `main.py` verifica e instala automaticamente** (via pip) qualquer dependência que faltar, etapa por etapa, conforme você for ativando-as no `config.txt`. Se preferir instalar tudo de uma vez:

```bash
pip install rdkit requests numpy scipy
pip install psycopg2-binary      # opcional -- só para a fonte DrugCentral (degrada sozinha sem ela)
pip install transformers torch    # opcional -- só para ranking_embeddings.py (pesado, alguns GB)
```

---

## Uso rápido

```bash
python main.py
```

Na primeira vez, se `config.txt` não existir, o `main.py` cria um modelo comentado e para — edite-o e rode de novo.

### O que o `main.py` faz em cada execução

1. Lê `config.txt`.
2. Cria a pasta `results/` se não existir.
3. Para cada etapa ativada: **verifica se as dependências dela estão instaladas e instala automaticamente as que faltarem** (via `pip install`), depois roda o script correspondente em `mods/` como subprocesso.
4. Se uma etapa obrigatória (1 ou 2) falhar — incluindo falha ao instalar uma dependência —, o pipeline para ali. Falhas nas etapas de validação (3, 4, 5) só emitem aviso, já que o ranking principal já foi gerado.

### `config.txt`

```ini
# Sintoma(s)/condição a pesquisar nas bases de referência (em inglês).
# Pode ser mais de um, separado por vírgula.
sintomas = neuropathic pain, peripheral neuropathy

# Pasta onde todos os resultados são salvos
results_dir = results

# --- farmaco_completo.py (input.txt -> farmaco.csv) ---
farmaco_input = input.txt
farmaco_output = farmaco.csv
pred_type = admet
batch_size = 5

# --- ranking_multibase.py (sintoma -> ranking.csv) ---
ranking_output = ranking.csv
top_por_fonte = 10
no_admet = false
pular_revisao_referencias = false

# --- validar_enriquecimento.py ---
n_decoys = 50

# --- ranking_embeddings.py (opcional -- requer transformers+torch, pesado) ---
modelo_embeddings = seyonec/ChemBERTa-zinc-base-v1
similaridade_embeddings = centroide
peso_embedding = 0.5
peso_admet_embeddings = 0.5

# --- Quais etapas rodar (main.py) ---
rodar_farmaco = true
rodar_ranking = true
rodar_comparacao_metodos = false
rodar_validacao_enriquecimento = false
rodar_embeddings = false
```

Linhas começando com `#` são comentários. Argumentos passados por linha de comando (quando você roda um script de `mods/` diretamente) sempre têm prioridade sobre o `config.txt` — o arquivo só preenche o que não foi passado explicitamente.

### As 5 etapas orquestradas

| # | Script | O que faz | Controlado por | Dependências verificadas |
|---|---|---|---|---|
| 1 | `mods/farmaco_completo.py` | `input.txt` → `results/farmaco.csv` | `rodar_farmaco` | rdkit, requests |
| 2 | `mods/ranking_multibase.py` | sintoma(s) → `results/ranking.csv` (+ revisão interativa) | `rodar_ranking` | rdkit, requests (+ psycopg2 opcional) |
| 3 | `mods/comparar_metodos_ranking.py` | validação vs. TOPSIS clássico/QED/Tanimoto | `rodar_comparacao_metodos` | rdkit, numpy, scipy |
| 4 | `mods/validar_enriquecimento.py` | validação leave-one-out + decoys | `rodar_validacao_enriquecimento` | numpy, scipy, requests |
| 5 | `mods/ranking_embeddings.py` | ranking alternativo via ChemBERTa | `rodar_embeddings` | numpy, **torch, transformers** |

---

## 1. `mods/farmaco_completo.py`

Une o cálculo de descritores de drogabilidade com a consulta ADMET em um único fluxo.

### Formato do `input.txt`

Um composto por linha, aceitando qualquer um destes formatos (pode misturar entre linhas):

```
CC(=O)Oc1ccccc1C(=O)O
Aspirina,CC(=O)Oc1ccccc1C(=O)O
Aspirina;CC(=O)Oc1ccccc1C(=O)O
```

Linhas vazias ou iniciadas com `#` são ignoradas.

### Uso direto (fora do `main.py`)

```bash
python mods/farmaco_completo.py --input input.txt --output results/farmaco.csv
```

### O que ele calcula

- **Descritores locais (RDKit):** peso molecular, doadores/aceptores de ligação H, TPSA, ligações rotacionáveis, refratividade molar, átomos pesados.
- **Regras de drogabilidade:** Regra dos 5 de Lipinski (parcial, sem LogP) e Regra de Veber.
- **ADMET (Deep-PK):** propriedades gerais, absorção, distribuição, metabolismo, excreção e toxicidade.

### Resiliência

- Salva progresso incrementalmente (`<output>.deeppk_raw.csv` + checkpoint) — pode ser interrompido e retomado.
- Se uma molécula específica for rejeitada pela API do Deep-PK, o script **isola e pula só ela** (fica sem dados ADMET, mas não trava o resto nem fica preso em loop).

---

## 2. `mods/ranking_multibase.py`

Ranqueia os compostos de um `farmaco.csv` por semelhança ao perfil químico de drogas já usadas para tratar um sintoma/condição.

### Uso direto

```bash
python mods/ranking_multibase.py --sintoma "neuropathic pain" --farmaco results/farmaco.csv --output results/ranking.csv
```

> **As buscas nas 4 bases são em inglês** e por correspondência literal de texto — siglas em português ou termos que não aparecem no texto/nome em inglês da condição não retornam nada.

### Como funciona

**Etapa 1 — Buscar drogas de referência em 4 bases (gratuitas, sem chave):**

| Fonte | O que traz | Critério de corte (top N) |
|---|---|---|
| **Open Targets** | Drogas conhecidas + fase clínica | Maior fase clínica (Aprovado > Fase IV > ... > Pré-clínico) |
| **ChEMBL** | Indicações por ID EFO/MONDO ou texto | Maior fase clínica da indicação |
| **DrugCentral** | Indicações já aprovadas (via Postgres público) | Correspondência exata ao termo, depois nº de produtos aprovados |
| **openFDA** | Menções em bula (`indications_and_usage`) | Nº de rótulos de bula diferentes que mencionam a condição |

Cada fonte é consultada **por inteiro primeiro** (sem cortar durante a busca) e só depois reduzida às `--top-por-fonte` mais relevantes.

**Etapa 2 — Mesclar as 4 listas:** deduplicação por **ID ChEMBL**, com fallback para nome normalizado nas fontes que não fornecem ID (DrugCentral, openFDA).

**Etapa 3 — Garantir o SMILES** de cada droga de referência: DrugCentral → ChEMBL pelo ID → PubChem pelo nome (pedindo `IsomericSMILES`, `CanonicalSMILES` e `ConnectivitySMILES` juntos).

**Revisão interativa** (a menos que `--pular-revisao-referencias` ou `pular_revisao_referencias = true`): mostra a lista numerada de drogas de referência e pergunta quais remover antes de prosseguir:

```
[ 1] Donepezila       (fontes: opentargets+chembl)
[ 2] Memantina        (fontes: chembl)
[ 3] Rivastigmina     (fontes: drugcentral)

Digite os números das drogas a REMOVER, separados por vírgula (ou Enter para manter todas):
```

Como o `main.py` chama esse script via subprocesso **sem** redirecionar entrada/saída, essa pergunta aparece normalmente no seu terminal mesmo rodando pelo orquestrador. Em ambiente não-interativo, é pulada automaticamente (mantém todas as referências).

**Etapa 4 — Calcular os descritores** dessas drogas via `processar_lista_druglikeness()`, importada do `farmaco_completo.py`.

**Etapa 5 — (opcional) ADMET das drogas de referência** via Deep-PK (`--no-admet` pula).

**Etapa 6 — Calcular os pesos** por característica, de acordo com a ocorrência entre as drogas de referência:

- **Numéricas:** peso maior quanto mais **consistente** o valor — `peso = 1 / (1 + desvio/média)`.
- **Booleanas** (Passa_Ro5, Passa_Veber): peso = frequência do valor mais comum.
- **Categóricas** (ADMET): extrai a categoria do texto e usa a frequência da mais comum.

**Etapa 7 — Pontuar e ranquear:** soma `peso × similaridade` por característica (gaussiana pra numéricas, match sim/não pra booleanas/categóricas), normaliza para 0–100.

> Este é o único método de ranking do script (`perfil`). Um método "TOPSIS híbrido" foi explorado e removido do projeto — mostrou-se estruturalmente instável mesmo após correções; a comparação com TOPSIS de verdade é coberta pelo `comparar_metodos_ranking.py` (abaixo).

### Arquivos de saída

- **`ranking.csv`** — ranking final (`Ranking`, `Nome`, `SMILES`, `Pontuacao_Final`).
- **`ranking.pesos.json`** — pesos calculados por característica.
- **`ranking.referencia.csv`** — drogas de referência usadas (já refletindo a revisão interativa), com a coluna `Fontes`.

---

## 3. `mods/comparar_metodos_ranking.py` — validação contra âncoras da literatura

```bash
python mods/comparar_metodos_ranking.py --farmaco results/farmaco.csv --referencia results/ranking.referencia.csv --pesos results/ranking.pesos.json --output results/comparacao_metodos.csv
```

Calcula 4 pontuações por composto:

| Método | Papel |
|---|---|
| `perfil` | O método usado no ranking principal |
| `topsis_classico` | **Comparador direto** — TOPSIS de manual, ideal/anti-ideal do próprio `farmaco.csv`, direção benefício/custo por Lipinski/Veber |
| `qed` | **Checagem de coerência** — drug-likeness geral (Bickerton et al., 2012) |
| `tanimoto_max` | **Checagem de coerência** — maior similaridade estrutural (Morgan/ECFP4) contra qualquer referência |

Gera uma matriz de correlação de Spearman entre todos os pares (`*.spearman.csv`). Correlação alta com `topsis_classico` é evidência de robustez; correlação baixa com `qed`/`tanimoto_max` não invalida o método por si só — medem construtos diferentes.

---

## 4. `mods/validar_enriquecimento.py` — leave-one-out com decoys

```bash
python mods/validar_enriquecimento.py --referencia results/ranking.referencia.csv --n-decoys 50 --output results/validacao_enriquecimento.csv
```

Responde: "uma nota alta significa algo, ou é só um número qualquer?" Metodologia padrão de triagem virtual retrospectiva (Mysinger et al., 2012, *J. Med. Chem.*):

1. Para cada droga de referência: remove ela do conjunto, **recalcula os pesos sem ela**.
2. Busca decoys no ChEMBL — drogas aprovadas para **outras** indicações (controle negativo).
3. Mistura a droga removida com os decoys, roda o ranking, registra o percentil de topo.
4. Repete para todas as referências.

Reporta: percentil médio/mediano, recuperação por corte (top 1%/5%/10%/20%, contagem bruta + Enrichment Factor), e ROC-AUC (droga conhecida vs. decoys).

---

## 5. `mods/ranking_embeddings.py` — ranking alternativo via embeddings aprendidos

```bash
python mods/ranking_embeddings.py --farmaco results/farmaco.csv --referencia results/ranking.referencia.csv --pesos results/ranking.pesos.json --output results/ranking_embeddings.csv
```

Desativado por padrão no `config.txt` (`rodar_embeddings = false`) porque `transformers`+`torch` são pesados. Ative com `rodar_embeddings = true` — o `main.py` instala as dependências automaticamente na primeira vez.

Usa ChemBERTa (`seyonec/ChemBERTa-zinc-base-v1` por padrão) para gerar embeddings, mede proximidade (centroide do grupo ou vizinho mais próximo) e reordena por ADMET:

```
embeddings das referências → embeddings dos candidatos → similaridade → ranking → reordenação por ADMET
```

| Chave do config.txt | Descrição | Padrão |
|---|---|---|
| `modelo_embeddings` | Checkpoint HuggingFace | `seyonec/ChemBERTa-zinc-base-v1` |
| `similaridade_embeddings` | `centroide` ou `max_individual` | `centroide` |
| `peso_embedding` / `peso_admet_embeddings` | Pesos do blend final | `0.5` / `0.5` |

**Limitação conhecida:** não foi validado ponta-a-ponta com o modelo real durante o desenvolvimento (ambiente sem acesso ao Hugging Face). A lógica de similaridade/ranking/reordenação foi testada com embeddings simulados; a primeira execução real (que baixa o modelo, alguns GB) precisa ser conferida.

---

## 6. `mods/calibrar_deeppk.py` — parada, não recomendada

Compara predições Deep-PK contra texto farmacocinético experimental do PubChem. **Decisão do projeto: não seguir com essa abordagem** (extração de números de texto livre foi julgada pouco confiável). Continua no repositório por referência, mas não é chamada pelo `main.py`.

---

## Dicas de busca por sintoma/condição

- Use sempre o **nome em inglês** (as 4 bases são todas em inglês).
- Para condições autoimunes pós-infecciosas (ex: Síndrome de Guillain-Barré ligada a *Campylobacter jejuni*), busque pelo **nome da síndrome resultante**, não pelo patógeno.
- Se o interesse é em **sintomas** específicos mais do que na doença de base, buscar pelo sintoma geral (`"neuropathic pain"`, `"peripheral neuropathy"`) costuma trazer mais dados (pesos mais confiáveis) do que o nome de uma condição rara.
- Vários termos podem ser combinados numa mesma chamada, mas termos muito distintos entre si podem diluir a especificidade do perfil de referência.

---

## Limitações conhecidas

- **DrugCentral** depende de conexão direta ao Postgres público; requer `psycopg2`. Se a tabela `product` não existir nesse dump, cai automaticamente para um critério de corte mais simples.
- **openFDA** não fornece SMILES nem ID ChEMBL — sempre usa o fallback do PubChem, e os nomes costumam vir como sal (ex: "DONEPEZIL HYDROCHLORIDE"), então pode não deduplicar com as outras fontes.
- Compostos biológicos, extratos vegetais/homeopáticos e sais inorgânicos simples encontrados como referência frequentemente não têm SMILES único válido e são descartados com aviso — esperado, não é bug.
- A API do Deep-PK pode rejeitar moléculas específicas; o script pula automaticamente (sem travar).
- **O método de ranking não é QSAR** no sentido formal (OECD, 2004) — não há endpoint de atividade biológica mensurado nem modelo de regressão treinado sobre ele. É um ranking por similaridade a um perfil de referência combinado com decisão multicritério, que incorpora predições de modelos QSAR/QSPR de terceiros (Deep-PK) como parte de suas variáveis de entrada.
- A fórmula de pesos (inverso do coeficiente de variação, invertido em relação ao uso convencional em MCDM) e a similaridade gaussiana não têm precedente direto na literatura — são adaptações do projeto, documentadas como tal. A gaussiana centrada na média é conceitualmente próxima de uma função de desejabilidade "alvo é o melhor" (Derringer & Suich, 1980), com agregação compensatória em vez da média geométrica não-compensatória clássica.
- **Calibração do Deep-PK contra dados experimentais:** não realizada (ver seção 6) — em aberto.
- **Nenhuma das validações (3, 4, 5) foi executada com dados reais do projeto** até o momento — os scripts estão testados com cenários sintéticos controlados. Rodar com `results/ranking.referencia.csv` real é o próximo passo para obter evidência de fato.
