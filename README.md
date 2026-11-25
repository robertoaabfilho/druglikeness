# Drogabilidade — Pipeline de Análise e Ranking de Compostos

Pipeline para (1) calcular propriedades de drogabilidade/ADMET de uma lista de compostos, (2) ranqueá-los por semelhança ao perfil químico de drogas já usadas para tratar um sintoma/condição específica, e (3) sinalizar compostos cuja bula indica que eles **causam** (em vez de tratar) esse mesmo efeito — com ferramentas de validação e um orquestrador que roda tudo a partir de um único arquivo de configuração.

## Estrutura de pastas

```
./main.py                  <- orquestrador (único arquivo solto na raiz)
./config.txt                <- você cria/edita (main.py gera um modelo na 1ª vez)
./input.txt                 <- você edita (lista de compostos)
./input.exemplo.tsv         <- 200 compostos reais, para testar sem montar entrada
./mods/                     <- todos os scripts do pipeline
    configuracao.py
    farmaco_completo.py
    ranking_multibase.py
    comparar_metodos_ranking.py
    validar_enriquecimento.py
    ranking_embeddings.py
    checar_alertas_adversos.py
    calibrar_deeppk.py
    importar_prospector.py
    prever_alvos.py
    validar_alvos_medidos.py
./exemplos/                 <- fixtures de validação (ver exemplos/LEIA-ME.md)
    controle.dor_neuropatica.tsv
    alvos.dor_neuropatica.txt
    alvos.patogeno.txt
    validacao.alvos_medidos.md
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

**Formato recomendado — TSV com cabeçalho iniciado por `#`:**

```
#nome	smiles	inchikey	inchikey_parental	fonte_smiles	conferencia	chembl_id
Cardamomin	COc1cc(O)cc(O)c1C(=O)/C=C/c1ccccc1	NYSZJNUIVUBQMM-BQYQJAHWSA-N	NYSZJNUIVUBQMM-BQYQJAHWSA-N	NPASS	mesma molécula	CHEMBL378104
```

Só `nome` e `smiles` são obrigatórios; as outras colunas são preservadas e
seguem para o `farmaco.csv`. Veja `input.exemplo.tsv` (200 compostos reais).

**Formatos simples** (uma molécula por linha, sem cabeçalho) continuam valendo:

```
CC(=O)Oc1ccccc1C(=O)O
Aspirina	CC(=O)Oc1ccccc1C(=O)O
Aspirina;CC(=O)Oc1ccccc1C(=O)O
Aspirina,CC(=O)Oc1ccccc1C(=O)O
```

> **Use TAB, não vírgula, para nome de produto natural.** O separador é
> procurado nesta ordem: TAB, `;`, `,` — e, para `;` e `,`, a divisão é feita
> a partir da **direita**, já que o SMILES é sempre o último campo. Isso é
> obrigatório porque nome de produto natural quase sempre tem vírgula
> (`1,2,3,4,6-pentagaloilglucose`, `3,3'-diindolilmetano`): numa base real,
> 42,8% dos nomes. Dividir pela *primeira* vírgula truncava o nome e
> corrompia o SMILES, que o RDKit então rejeitava — invalidando o composto
> em silêncio.

### Trazendo compostos de um pipeline de prospecção

`mods/importar_prospector.py` converte a tabela de compostos de um pipeline
de prospecção de plantas (`resultados/atividades_compostos_parcial.csv`) para
o formato acima:

```bash
python mods/importar_prospector.py \
    --origem ~/prospector/pipeline_prospeccao_v4/resultados/atividades_compostos_parcial.csv \
    --saida input.txt

# amostra de teste, metade com evidência ChEMBL (útil para validação)
python mods/importar_prospector.py --origem ... --saida input.exemplo.tsv \
    --limite 200 --preferir-com-chembl
```

Ele replica o tratamento químico da prospecção, de forma independente (os dois
projetos são repositórios separados — nada do código de lá é importado):

- **SMILES:** escolhe o da primeira fonte cujo SMILES o RDKit consegue ler, na
  ordem **NPASS > COCONUT > LOTUS > PubChem**, e grava a forma canônica.
- **Conferência:** preserva o veredito da comparação entre as bases
  (`mesma molécula`, `estereoquímica diferente: …`, `estruturas diferentes: …`).
  Compostos divergentes entram **marcados, não descartados** — a decisão é do
  pesquisador.
- **Identidade:** recalcula dois InChIKeys — `inchikey` da molécula como está e
  `inchikey_parental` depois de *Cleanup → maior fragmento orgânico →
  neutralização*. **A chave parental é a correta para cruzar as duas bases**,
  porque reconhece que sal e base livre são o mesmo composto, e a deduplicação
  é feita por ela. Exemplo real: a neferina tem dois `[NH+]` no SMILES, então
  seu `inchikey` sai `MIBATSHDJRIUJK-ROJLCIKYSA-P`, enquanto o `-…-N` que as
  bases declaram só aparece na forma parental.
- Divergência entre o InChIKey declarado e o recalculado também é reportada em
  `conferencia`, sem descartar a linha.

> **Procedência dos dados:** `input.exemplo.tsv` traz estruturas públicas de
> NPASS, COCONUT, LOTUS e PubChem, trazidas pela prospecção. Ao publicar
> resultados, cite as bases de origem.

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

## 8. `mods/prever_alvos.py` — compostos × alvos (qualquer organismo)

Responde à pergunta central do repositório: **dada uma lista de compostos e uma
lista de alvos, quais compostos têm chance de atuar sobre esses alvos?**

A lista de alvos é **entrada**, não dedução. Quem a fornece (ou a enriquece a
partir da doença) é o pesquisador.

Usa os pesos pré-treinados do **PP4CKD** ([Zenodo
10.5281/zenodo.23147411](https://zenodo.org/records/23147411), pesos CC-BY-4.0,
código Apache-2.0): uma MLP que recebe o fingerprint de uma molécula e devolve
uma probabilidade para cada um dos 7.676 alvos do ChEMBL 36.

```
SMILES -> ECFP4 (Morgan r=2, 4096 bits) -> Linear(4096,1000) -> ReLU
       -> Linear(1000,500) -> ReLU -> Linear(500,7676) -> sigmoide
```

### Alvos não humanos são de primeira classe

Dos 7.676 alvos, **4.140 são de *Homo sapiens* e 3.536 não são** — 666 de
camundongo, 646 de rato, e centenas de patógeno, distribuídos por 1.084
organismos. O modelo traz até alvos de **organismo inteiro** (ensaio
fenotípico, não proteína): `CHEMBL360` = *M. tuberculosis* com 11.049 moléculas
ativas no treino.

Isso serve a dois usos diferentes:

| uso | o alvo é | exemplo |
|---|---|---|
| **afirmativo** | do patógeno | doença infecciosa/parasitária: `CHEMBL360`, `CHEMBL364` (*P. falciparum*), `CHEMBL368` (*T. cruzi*) |
| **negativo** | atividade indesejada | excluir do ranking o composto que é antibacteriano, não analgésico |

O uso negativo foi **medido** no conjunto-controle deste repositório (6 fármacos
que tratam dor neuropática contra 6 que a causam):

| escore | AUC |
|---|---|
| só alvos humanos de dor | 0,750 |
| `p_humano × (1 − p_bacteriano)` | **0,806** |
| só `CHEMBL360` separando antibacteriano de analgésico | **1,000** |

Com o alvo bacteriano no cálculo, o ciprofloxacino cai do 4º para o 7º lugar.

### Como pedir os alvos

Um por linha. `#` comenta a linha inteira ou o fim dela. `| organismo` fixa o
organismo **daquela linha**:

```
CHEMBL360                      # M. tuberculosis, célula inteira
CHEMBL1849                     # InhA de Mtb
DNA gyrase | Escherichia coli  # nome + organismo, para não pegar o ortólogo errado
dihydrofolate reductase | tuberculose
```

Ou restrinja tudo com `--organismo` (repetível, aceita vírgula). Apelidos
(`humano`, `rato`, `camundongo`, `tuberculose`, `malaria`, `chagas`…) e trechos
do nome científico funcionam — `--organismo Plasmodium` alcança
`Plasmodium falciparum (isolate 3D7)`, porque as bases guardam cepa e isolado
dentro do nome do organismo.

**Pedir por ID ChEMBL vence o filtro de organismo.** Se você escreveu
`CHEMBL304`, você quer o transportador de noradrenalina **de rato** — o script
avisa e mantém.

### Descobrindo a grafia (sem rodar previsão)

Sem isto não há como escrever uma lista de alvos não humanos, porque a grafia
tem de ser a do ChEMBL:

```bash
python mods/prever_alvos.py --rotulos target_labels.tsv --listar-organismos
python mods/prever_alvos.py --rotulos target_labels.tsv --listar-organismos Plasmodium
python mods/prever_alvos.py --rotulos target_labels.tsv \
    --buscar-alvo gyrase --organismo Staphylococcus
```

### Uso

```bash
python mods/prever_alvos.py \
    --entrada input.exemplo.tsv \
    --modelos ~/modelos_baixados/pp4ckd/models \
    --rotulos ~/modelos_baixados/pp4ckd/target_labels.tsv \
    --alvos alvos.txt \
    --organismo humano \
    --saida results/alvos_preditos.tsv
```

Saída (TSV): `nome`, `inchikey_parental`, `chembl_id_composto`,
`alvo_chembl_id`, `alvo`, `organismo`, `tipo_alvo`, `moleculas_treino`,
`probabilidade`, `procedencia`.

### Duas colunas que existem por um motivo

**`moleculas_treino`** — quantas moléculas ativas treinaram aquele alvo. Alvo
não humano é bem mais esparso: mediana de **26** moléculas ativas contra 67 dos
humanos, 67% abaixo de 50, e 798 alvos abaixo de 10. Sem essa coluna você não
distingue p=0,89 sustentado por 3.493 exemplos de p=0,89 sustentado por 10. O
caso real: os três alvos mais prováveis do ciprofloxacino no modelo inteiro são
organismos com 10, 11 e 12 moléculas de treino, enquanto a DNA girase (596
ativos, o alvo de verdade) fica mais abaixo. `--minimo-treino N` (padrão 50)
avisa — **sem remover** — quando um alvo pedido está abaixo do limiar.

**`tipo_alvo`** — `SINGLE PROTEIN` é ligação a uma proteína; `ORGANISM` e
`CELL-LINE` são ensaio fenotípico de célula inteira. São afirmações diferentes e
não devem ser lidas como a mesma coisa.

### Procedência: predito nunca se mistura com medido

```
predito (PP4CKD ECFP4 ecfp4_fold_4.pt)   previsão de modelo, hipótese a testar
medido (ChEMBL pchembl=7.2)              bioatividade experimental
```

Uma previsão de 0,93 não tem o mesmo peso que um ensaio, e as duas nunca saem
na mesma coluna sem distinção. A coluna `procedencia` nomeia o checkpoint usado,
para que o resultado seja reproduzível.

### Limites conhecidos

- **Pró-fármaco.** O modelo é baseado em ligante: julga a estrutura que recebe.
  A etionamida pontua 0,80 em *M. tuberculosis* (fenótipo) e **0,0017 em InhA**
  (a proteína), porque só inibe a InhA depois de ativada pela EthA — a
  molécula-mãe não se liga.
- **Molécula pequena e polar.** Gabapentina e pregabalina pontuam baixo nos
  alvos de dor (0,0001 e 0,02): são análogos de aminoácido, e o alvo delas
  (α2δ-1) tem só 562 ativos no treino — a sub-unidade α2δ-2 tem 30. A AUC de
  0,750 do conjunto humano vem quase toda dessas duas falhas.
- **Ortólogo silencioso** é tratado, não eliminado: só o filtro de organismo e a
  coluna `organismo` evitam que `noradrenaline transporter` traga humano
  (`CHEMBL222`), camundongo (`CHEMBL2370`) e rato (`CHEMBL304`) juntos.

### Os pesos publicados estão truncados na origem

O `pp4ckd_models.tar.gz` do Zenodo tem **MD5 idêntico ao publicado**
(`f9380bb9cef3fe56d6d49d3e4a4d57ec`, 785.350.656 bytes) mas o tar termina antes
do fim: só **22 dos ~77 modelos** saem íntegros, e `ecfp4_full_model.pt` não é
um deles. Não é download corrompido — é o arquivo publicado.

Por isso:

- `--checkpoint arquivo.pt` aceita um `.pt` específico;
- a pasta `--modelos` cai automaticamente para `<fp>_fold_N.pt`, com aviso.

Um checkpoint de fold é um modelo treinado de verdade, em 9/10 dos dados —
válido para previsão e **preferível** na validação, porque permite medir fora do
que o modelo viu. Os resultados acima usam `ecfp4_fold_4.pt`.

## 9. `mods/validar_alvos_medidos.py` — recall@10 contra alvos medidos

Mede o `prever_alvos.py` contra bioatividade **medida** (não âncoras de
literatura, não decoys): quando um composto tem atividade experimental contra um
alvo, o modelo põe esse alvo no topo?

Resultados completos em [`exemplos/validacao.alvos_medidos.md`](exemplos/validacao.alvos_medidos.md).
Em resumo, sobre 6.000 compostos com alvo medido no ChEMBL:

- **O vazamento não foi eliminado.** A divisão de folds do PP4CKD
  (`cv_splits_10fold.npz`) não foi publicada, e o estrato abaixo do limiar de
  treino (`pchembl < 5`) tem só 21 pares depois de remover triagem HTS. Todo
  número é limite superior — está dito lá, não escondido.
- **O treino do alvo domina**: recall@10 vai de 0,127 a 0,477 conforme o alvo
  sai de <50 para ≥5.000 moléculas ativas. Fator de 3,8, e a razão de existir a
  coluna `moleculas_treino`.
- **Alvo não humano vai melhor**: 0,387 contra 0,209 do humano.
- **Triagem HTS é ruído**: 0,111 contra 0,272 de ensaio dirigido; entra como
  estrato próprio, nunca somada ao resto.
- **A medida que corresponde ao uso real**: o alvo medido fica na posição
  mediana **79 de 7.676** (top 1,0%); o melhor alvo de cada composto, na posição
  3. Numa lista fornecida de 10 alvos, o verdadeiro ficaria no top-3 em ~95% dos
  casos — teto otimista, porque alvos de uma mesma doença são correlacionados.

O recall@10 cru de 0,24 **não** é comparável ao 0,856 do artigo: o regime deles
tem ~1,3 rótulos por composto, aqui a mediana é 2 com cauda até 846 (painéis de
triagem), e o recall@10 é limitado por construção acima de 10 rótulos.

## 10. `mods/coletar_evidencias.py` — enriquecimento de alvos (Open Targets e Reactome)

Com `--alvos-enriquecidos alvos_enriquecidos.tsv --rotulos target_labels.tsv`, a etapa 1 grava uma segunda saída: a lista do usuário **mais** os alvos que o Open Targets associa aos sintomas e os que o Reactome acrescenta por via. O enriquecimento só acrescenta: nenhum alvo do usuário sai, muda de ordem ou é sobrescrito.

- **Sintoma -> EFO nunca é escolhido em silêncio.** Fixe com `--efo "neuropathic pain=EFO_0005762"` (repetível). Sintoma sem `--efo` não é enriquecido; a CLI lista os candidatos do Open Targets e pede a escolha.
- **Open Targets:** alvos com escore de associação >= `--corte-associacao` (padrão **0.4**, ajustável e gravado no manifesto), UniProt -> ChEMBL (`SINGLE PROTEIN` primeiro). Para a dor neuropática (934 alvos), 0.5 deixa 30 e 0.4 deixa 44; os 14 a mais são os canais de sódio dependentes de voltagem (SCN9A/Nav1.7 entre eles) e os adrenoceptores alfa-2, alvos reconhecidos de dor neuropática. O corte escolhe pelo escore do próprio alvo: promover o que o Reactome alcançou abaixo do corte faria a `origem` depender de acaso de via, e por isso não se faz; a `observacao` apenas avisa.
- **Colapso de família:** só dobra num representante os parálogos **ausentes do modelo**, na mesma classe de proteína do ChEMBL; o ID dobrado vai para `familia_colapsados` e a classe para `familia`. Parálogo presente no modelo continua sendo linha própria (NET e SERT, ambos em "SLC06 neurotransmitter transporter", são o mecanismo duplo da duloxetina).
- **Reactome:** só dados de anotação (CC0): `data/mapping/UniProt/{acc}/pathways` e `data/participants/{stId}/referenceEntities`. Os arquivos de Functional Interaction não são usados (licença da Pathway Solutions, conteúdo KEGG). Vias com mais de `--teto-via` proteínas (padrão 150) são ignoradas, e o alvo expandido só entra se existir no espaço do modelo. Cada adição traz `origem = reactome`, `via_origem` (todas as vias que a trouxeram, separadas por `;`) e `alvo_semente`. Alvo que o Reactome alcança e que o Open Targets também associa ao sintoma, mas abaixo do corte, continua `reactome` e ganha o escore e uma `observacao`.

**Cache de rede e repeticao.** Toda resposta de ChEMBL, Open Targets e Reactome vai para `--cache DIR` (padrao `cache_rede/` ao lado da saida), numa subpasta por versao dos servicos (`meta.dataVersion` do Open Targets, versao do ChEMBL e do Reactome); versao nova nunca le o cache velho. Erro, pagina truncada e timeout nunca sao gravados. `--manifesto manifesto.json` entra em **modo de repeticao**: sem rede, serve tudo do cache que o manifesto aponta (falta no cache e erro) e recusa rodar se mudou a versao dos servicos, o MD5 de checkpoint, rotulos, evidencias medidas, ligantes, alvos ou sintomas, ou algum parametro que muda a saida (`--corte`, `--limite`, `--corte-associacao`, `--teto-via`, `--efo`, `--sem-rede`). O novo manifesto sai ao lado da saida (`manifesto_repeticao.json` se colidir) e nunca sobrescreve o lido.

**Mudanca de `--manifesto`:** antes era o caminho do manifesto a ESCREVER; agora e o manifesto a REPETIR, e o manifesto de saida fica sempre em `manifesto.json` ao lado da saida. Quem passava `--manifesto saida.json` agora recebe erro de leitura, nao uma execucao silenciosamente diferente.

**Fluxo em dois passos.** A coleta de evidência **não** usa os alvos enriquecidos: a lista precisa ser revisada antes de gastar a previsão sobre todos os ligantes (210 alvos x ~109 mil ligantes são ~23 milhões de pares). A execução avisa isso ao final e o manifesto registra `alvos_enriquecidos_usados_na_evidencia: false`.

```bash
# 1. enriquecer (escreve alvos_enriquecidos.tsv; evidencias.tsv usa so os alvos informados)
python mods/coletar_evidencias.py --ligantes ligantes.tsv --alvos alvos.tsv --sintomas sintomas.txt \
    --saida passo1/evidencias.tsv --rotulos target_labels.tsv \
    --alvos-enriquecidos alvos_enriquecidos.tsv --efo "neuropathic pain=EFO_0005762"
# 2. revise/edite alvos_enriquecidos.tsv e colete a evidencia com ele
python mods/coletar_evidencias.py --ligantes ligantes.tsv --alvos alvos_enriquecidos.tsv --sintomas sintomas.txt \
    --saida evidencias.tsv --checkpoint ecfp4.pt --rotulos target_labels.tsv
```

**Quanto cada filtro corta, medido** (sementes CHEMBL222/NET, CHEMBL228/SERT, CHEMBL233/MOR, CHEMBL4794/TRPV1; modelo PP4CKD; Reactome e ChEMBL consultados em 2026-11-25):

| etapa | proteínas |
|---|---|
| nas vias das 4 sementes, sem filtro | 573 (+4 sementes) |
| após o teto de 150 por via (descarta "G alpha (i) signalling events", 317, e "Peptide ligand-binding receptors", 200) | 259 |
| sem alvo humano no ChEMBL | -53 |
| com alvo no ChEMBL, mas fora do espaço do modelo | -61 |
| **entram no enriquecido** | **145** |

O filtro de presença no modelo corta 114 das 259 proteínas (44%); os dois filtros juntos levam 573 a 145 (75% a menos). Os números dependem das sementes e do `target_labels.tsv`.

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
