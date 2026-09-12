# Drogabilidade — Pipeline de Análise e Ranking de Compostos

Pipeline para (1) calcular propriedades de drogabilidade/ADMET de uma lista de compostos, (2) ranqueá-los por semelhança ao perfil químico de drogas já usadas para tratar um sintoma/condição específica, e (3) sinalizar compostos cuja bula indica que eles **causam** (em vez de tratar) esse mesmo efeito — com ferramentas de validação e um orquestrador que roda tudo a partir de um único arquivo de configuração.

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
    checar_alertas_adversos.py
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
    ranking.alertas.csv             (se rodar_alertas_adversos = true)
```

## Visão geral do fluxo

```
input.txt ──► mods/farmaco_completo.py ──► results/farmaco.csv
                                                  │
config.txt (sintomas) ──► mods/ranking_multibase.py ──► results/ranking.csv
                                                  │        results/ranking.referencia.csv
                                                  │        results/ranking.pesos.json
                    ┌─────────────┬───────────────┼───────────────┐
                    ▼             ▼               ▼               ▼
     comparar_metodos_ranking.py  │  validar_enriquecimento.py  ranking_embeddings.py
        (validação: 4 métodos)    │  (validação: leave-one-out)  (ranking alternativo)
                                   ▼
                     checar_alertas_adversos.py
                (sinaliza compostos que CAUSAM o sintoma-alvo)
```

`mods/calibrar_deeppk.py` existe como ferramenta complementar, fora do fluxo do `main.py` (decisão do projeto: não seguir com essa abordagem).

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
3. Para cada etapa ativada: **verifica se as dependências dela estão instaladas e instala automaticamente as que faltarem**, depois roda o script correspondente em `mods/` como subprocesso.
4. Se uma etapa obrigatória (1 ou 2) falhar, o pipeline para ali. Falhas nas etapas de validação/checagem (3 a 6) só emitem aviso, já que o ranking principal já foi gerado.

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

# --- checar_alertas_adversos.py (opcional -- checa bula por efeitos adversos) ---
# Um ou mais termos, separados por vírgula. Qualquer um que aparecer nas
# seções de advertência da bula já sinaliza o composto.
termos_alerta = neuropathy, neuritis, nerve damage
penalidade_alerta = 0.5

# --- Quais etapas rodar (main.py) ---
rodar_farmaco = true
rodar_ranking = true
rodar_comparacao_metodos = false
rodar_validacao_enriquecimento = false
rodar_embeddings = false
rodar_alertas_adversos = false
```

Linhas começando com `#` são comentários. Argumentos passados por linha de comando (quando você roda um script de `mods/` diretamente) sempre têm prioridade sobre o `config.txt`.

### As 6 etapas orquestradas

| # | Script | O que faz | Controlado por | Dependências verificadas |
|---|---|---|---|---|
| 1 | `mods/farmaco_completo.py` | `input.txt` → `results/farmaco.csv` | `rodar_farmaco` | rdkit, requests |
| 2 | `mods/ranking_multibase.py` | sintoma(s) → `results/ranking.csv` (+ revisão interativa) | `rodar_ranking` | rdkit, requests (+ psycopg2 opcional) |
| 3 | `mods/comparar_metodos_ranking.py` | validação vs. TOPSIS clássico/QED/Tanimoto | `rodar_comparacao_metodos` | rdkit, numpy, scipy |
| 4 | `mods/validar_enriquecimento.py` | validação leave-one-out + decoys | `rodar_validacao_enriquecimento` | numpy, scipy, requests |
| 5 | `mods/ranking_embeddings.py` | ranking alternativo via ChemBERTa | `rodar_embeddings` | numpy, **torch, transformers** |
| 6 | `mods/checar_alertas_adversos.py` | sinaliza compostos que causam (em vez de tratar) o sintoma-alvo | `rodar_alertas_adversos` | requests |

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

### O que ele calcula

- **Descritores locais (RDKit):** peso molecular, doadores/aceptores de ligação H, TPSA, ligações rotacionáveis, refratividade molar, átomos pesados.
- **Regras de drogabilidade:** Regra dos 5 de Lipinski (parcial, sem LogP) e Regra de Veber.
- **ADMET (Deep-PK):** propriedades gerais, absorção, distribuição, metabolismo, excreção e toxicidade.

### Resiliência

- Salva progresso incrementalmente (`<output>.deeppk_raw.csv` + checkpoint) — pode ser interrompido e retomado.
- Se uma molécula específica for rejeitada pela API do Deep-PK, o script **isola e pula só ela**.

> **Observação de execuções reais:** cerca de 30-35% dos compostos submetidos podem ficar sem nenhum dado ADMET (rejeitados pela API).

---

## 2. `mods/ranking_multibase.py`

Ranqueia os compostos de um `farmaco.csv` por semelhança ao perfil químico de drogas já usadas para tratar um sintoma/condição.

> **As buscas nas 4 bases são em inglês** e por correspondência literal de texto.

### Como funciona

**Etapa 1 — Buscar drogas de referência em 4 bases (gratuitas, sem chave):**

| Fonte | O que traz | Critério de corte (top N) |
|---|---|---|
| **Open Targets** | Drogas conhecidas + fase clínica | Maior fase clínica |
| **ChEMBL** | Indicações por ID EFO/MONDO ou texto | Maior fase clínica da indicação |
| **DrugCentral** | Indicações já aprovadas (via Postgres público) | Correspondência exata ao termo, depois nº de produtos aprovados |
| **openFDA** | Menções em bula (`indications_and_usage`) | Nº de rótulos de bula diferentes que mencionam a condição |

**Etapa 2 — Mesclar as 4 listas:** deduplicação por **ID ChEMBL**, com fallback para nome normalizado nas fontes que não fornecem ID.

> **Correção aplicada:** a mesma droga podia aparecer duas vezes quando encontrada **com** ID ChEMBL numa fonte e **sem** ID em outra. Um índice nome→chave agora funde as duas entradas independente da ordem em que as fontes chegam — confirmado corrigindo um caso real (Pregabalina, Duloxetina, Tapentadol, Capsaicina apareciam duplicadas).

**Etapa 3 — Garantir o SMILES** de cada droga de referência: DrugCentral → ChEMBL pelo ID → PubChem pelo nome.

**Revisão interativa** (a menos que `--pular-revisao-referencias`): mostra a lista numerada de drogas de referência e pergunta quais remover antes de prosseguir. Funciona normalmente mesmo via `main.py`.

**Etapa 4 — Calcular os descritores** dessas drogas via `processar_lista_druglikeness()`.

**Etapa 5 — (opcional) ADMET das drogas de referência** via Deep-PK.

**Etapa 6 — Calcular os pesos** por característica (numéricas: consistência; booleanas/categóricas: frequência da moda).

**Etapa 7 — Pontuar e ranquear:** soma `peso × similaridade` por característica, normaliza para 0–100.

> Este é o único método de ranking do script (`perfil`). Um método "TOPSIS híbrido" foi explorado e removido do projeto por instabilidade; a comparação com TOPSIS de verdade é coberta pelo `comparar_metodos_ranking.py`.

### Arquivos de saída

- **`ranking.csv`**, **`ranking.pesos.json`**, **`ranking.referencia.csv`** (com a coluna `Fontes`).

---

## 3. `mods/comparar_metodos_ranking.py` — validação contra âncoras da literatura

| Método | Papel |
|---|---|
| `perfil` | O método usado no ranking principal |
| `topsis_classico` | **Comparador direto** — TOPSIS de manual |
| `qed` | **Checagem de coerência** — drug-likeness geral |
| `tanimoto_max` | **Checagem de coerência** — proximidade estrutural às referências |

**Resultado de uma execução real:** correlações moderadas e significativas com os três (ρ 0,28–0,48, p<0,05).

---

## 4. `mods/validar_enriquecimento.py` — leave-one-out com decoys

Remove cada droga de referência, recalcula pesos sem ela, mistura com decoys do ChEMBL, roda o ranking, mede o percentil de topo.

**Resultado de uma execução real (N=25):** percentil médio 77%, mediano 90%. Enrichment Factor de 8,0× (top 1%) a 4,8× (top 5%) acima do acaso.

---

## 5. `mods/ranking_embeddings.py` — ranking alternativo via embeddings aprendidos

Desativado por padrão (`rodar_embeddings = false`, requer `transformers`+`torch`). Usa ChemBERTa para gerar embeddings, mede proximidade e reordena por ADMET.

**Resultado de uma execução real:** executado com sucesso; o top-14 produzido diverge do método `perfil`.

---

## 6. `mods/checar_alertas_adversos.py` — compostos que causam (em vez de tratar) o sintoma-alvo

### O problema que motivou esse script

O método `perfil` mede semelhança físico-química/ADMET, não mecanismo de ação. Numa execução real, o ranking pra "dor neuropática" trouxe no top 10 **Ciprofloxacino** (tarja preta da FDA para neuropatia periférica), **Linezolida**, **Etionamida**, **Cloranfenicol** e **Etambutol** — drogas famosas por **causar** neuropatia como efeito adverso, não por tratá-la. Drogas que tratam dor neuropática e drogas que a causam como toxicidade podem convergir no mesmo espaço físico-químico (ambas precisam penetrar tecido neural), sem nenhuma relação farmacológica real entre elas.

### Como funciona

```
python mods/checar_alertas_adversos.py --ranking results/ranking.csv --termos "neuropathy,neuritis" --output results/ranking.alertas.csv
```

1. Pra cada composto do ranking, busca a bula no openFDA (por nome genérico, princípio ativo ou comercial).
2. Procura o(s) termo(s) de alerta **somente** em `boxed_warning`, `warnings`, `warnings_and_cautions` e `adverse_reactions` — **nunca** em `indications_and_usage` (testado explicitamente: uma droga que *trata* o sintoma não é confundida com uma que o *causa*).
3. Aplica uma penalidade multiplicativa (`penalidade_alerta`, padrão 0.5) na pontuação dos sinalizados, e gera um ranking ajustado.

### Múltiplos termos de alerta

`termos_alerta` (ou `--termos`) aceita vários termos separados por vírgula — qualquer um que aparecer já sinaliza o composto:

```
termos_alerta = neuropathy, neuritis, nerve damage, peripheral nerve
```

O CSV de saída (`Trecho_Alerta`) sempre mostra qual termo específico bateu e o trecho exato da bula, mesmo usando vários termos.

### Por que a penalidade é parcial, não uma remoção

Busca em texto livre não entende negação — uma bula que diz "não há relatos de neuropatia" ainda bateria no termo. Por isso os compostos sinalizados continuam visíveis no CSV (não somem silenciosamente), com o trecho da bula exposto para revisão humana antes de descartar qualquer composto.

### Colunas de saída (`ranking.alertas.csv`)

`Ranking_Ajustado`, `Ranking` (original), `Nome`, `SMILES`, `Pontuacao_Final` (original), `Pontuacao_Ajustada`, `Alerta_Adverso` (SIM/nao), `Campo_Alerta`, `Trecho_Alerta`.

---

## 7. `mods/calibrar_deeppk.py` — parada, não recomendada

Compara predições Deep-PK contra texto farmacocinético experimental do PubChem. **Decisão do projeto: não seguir com essa abordagem.** Continua no repositório por referência, mas não é chamada pelo `main.py`.

---

## Dicas de busca por sintoma/condição

- Use sempre o **nome em inglês**.
- Para condições autoimunes pós-infecciosas, busque pelo **nome da síndrome resultante**, não pelo patógeno.
- Sintomas gerais costumam trazer mais dados que condições raras específicas.
- **Ative `rodar_alertas_adversos` sempre que o sintoma-alvo também for um efeito adverso conhecido de outras drogas** (dor, neuropatia, convulsão, arritmia etc.) — é exatamente esse tipo de sintoma que corre risco de trazer "causadores" no topo do ranking em vez de "tratadores".

---

## Limitações conhecidas

- **DrugCentral** depende de conexão direta ao Postgres público; requer `psycopg2`. Se a tabela `product` não existir nesse dump, cai automaticamente para um critério de corte mais simples.
- **openFDA** não fornece SMILES nem ID ChEMBL — usa o fallback do PubChem, e os nomes costumam vir como sal.
- Compostos biológicos, extratos vegetais/homeopáticos e sais inorgânicos simples frequentemente não têm SMILES único válido e são descartados com aviso.
- A API do Deep-PK pode rejeitar moléculas específicas; o script pula automaticamente. Em execuções reais, isso afetou ~30-35% dos candidatos.
- **O método de ranking não é QSAR** no sentido formal (OECD, 2004) — é um ranking por similaridade a um perfil de referência combinado com decisão multicritério, que incorpora predições de modelos QSAR/QSPR de terceiros (Deep-PK) como parte de suas variáveis de entrada.
- A fórmula de pesos e a similaridade gaussiana não têm precedente direto na literatura — são adaptações do projeto. A gaussiana centrada na média é conceitualmente próxima de uma função de desejabilidade "alvo é o melhor" (Derringer & Suich, 1980).
- **A checagem de alertas adversos (`checar_alertas_adversos.py`) não entende negação em texto livre** — trechos sinalizados precisam de revisão humana antes de descartar um composto; por isso a penalidade é parcial por padrão, não uma exclusão.
- **A validação por decoys não testa especificidade entre classes farmacológicas** — mostra que drogas de referência conhecidas superam controles negativos aleatórios, mas não testa diretamente se compostos de uma classe não relacionada estariam sendo indevidamente favorecidos. A checagem de alertas adversos (etapa 6) é um primeiro passo nessa direção, mas específico a efeitos já documentados em bula, não uma solução geral.
- **Calibração do Deep-PK e comparação formal com outras ferramentas ADMET continuam sem execução.**
