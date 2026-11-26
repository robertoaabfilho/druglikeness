# Drogabilidade — Pipeline de Análise e Ranking de Compostos

Pipeline para (1) calcular propriedades de drogabilidade/ADMET de uma lista de compostos, (2) ranqueá-los por semelhança ao perfil químico de drogas já usadas para tratar um sintoma/condição específica, e (3) sinalizar compostos cuja bula indica que eles **causam** (em vez de tratar) esse mesmo efeito — com ferramentas de validação e um orquestrador que roda tudo a partir de um único arquivo de configuração.

## Estrutura de pastas

```
./main.py                  <- orquestrador (único arquivo solto na raiz)
./config.txt                <- você cria/edita (main.py gera um modelo na 1ª vez)
./input.txt                 <- você edita (lista de compostos)
./input.exemplo.tsv         <- 200 compostos reais, para testar sem montar entrada
./requirements.txt          <- instalavel num clone limpo
./requirements-opcionais.txt <- pesadas: torch, transformers, psycopg2-binary
./painel.py                 <- painel Streamlit da etapa 1
./dados/                    <- dado derivado do ChEMBL (CC BY-SA 3.0, dados/LICENSE)
    evidencias_controle.tsv
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
    coletar_evidencias.py    <- CLI da etapa 1 (seção 11)
    aplicar_portao.py        <- CLI da etapa 2 (seção 12)
    evidencias/              <- módulos da etapa 1
    portao/                  <- módulos da etapa 2
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

**O fluxo acima é o original, e está sendo substituído.** O motivo é medido: ele
ranqueia por semelhança físico-química com as drogas de referência, sem olhar
alvo nenhum, e no conjunto-controle de 12 fármacos os dois **QED mais altos
(0,89)** são ciprofloxacino e linezolida — ambos do grupo que **causa** dor
neuropática, enquanto os que a tratam ficam entre 0,54 e 0,81. O fluxo novo
separa o que estava misturado:

```
etapa 1  COLETAR     o que se sabe        -> nunca decide    (seção 11)
etapa 2  O PORTÃO    o que fica marcado   -> quase nunca reprova (seção 12)
fase 3   RANQUEAR    a ordem usando os alvos -> ainda não existe
```

Cada etapa tem CLI e manifesto próprios e não passa pelo `main.py`.

---

## Requisitos

Só `rdkit` e `requests` são estritamente necessários para começar — **o `main.py` verifica e instala automaticamente** (via pip) qualquer dependência que faltar, etapa por etapa, conforme você for ativando-as no `config.txt`. Se preferir instalar tudo de uma vez:

```bash
pip install -r requirements.txt              # instala num clone limpo, sem compilador
pip install -r requirements-opcionais.txt    # pesadas: torch, transformers, psycopg2-binary
```

Os dois arquivos são separados porque o pip **instala toda linha que encontrar**,
comentada como "opcional" ou não. Com `psycopg2` (sem `-binary`) no
`requirements.txt`, o primeiro comando do README falhava num clone limpo por
falta do `pg_config`; e `torch` + `transformers` somam alguns GB que nem o
painel nem os testes precisam. O que cada opcional liga está dito no próprio
`requirements-opcionais.txt`, e o código avisa na hora qual instalar:

| opcional | liga |
|---|---|
| `torch` | evidência **prevista** (PP4CKD). Sem ele a etapa 1 roda em modo degradado, só com evidência medida |
| `transformers` | `mods/ranking_embeddings.py` (ChemBERTa) |
| `psycopg2-binary` | fonte DrugCentral do `ranking_multibase.py`, que degrada sozinha sem ela |

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
| 6 | `mods/checar_alertas_adversos.py` | sinaliza compostos que causam (em vez de tratar) o sintoma-alvo — **substituído pela etapa 2** (seção 12) | `rodar_alertas_adversos` | requests |

As etapas 1 e 2 (seções 11 e 12) **não** são orquestradas pelo `main.py`: elas
têm CLI própria e manifesto próprio, e são o fluxo para o qual o pipeline está
sendo reestruturado.

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

> **O caminho Deep-PK está SUBSTITUÍDO pelo ADMET-AI da etapa 2** (seção 12).
> O Deep-PK é uma API externa que rejeita 30-35% dos compostos, não tem versão
> citável e não repete: a mesma molécula noutro dia pode voltar outro número.
> O ADMET-AI 2.0.1 roda **local** (MIT), devolve **104 colunas** para toda
> molécula que o RDKit lê, e a versão entra no manifesto como chave de
> repetição. Os descritores do RDKit desta seção continuam válidos — é só a
> consulta ADMET remota que foi substituída.

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

> **SUBSTITUÍDO pelo critério `alerta_bula` da etapa 2** (seção 12). Continua no
> repositório porque o `main.py` antigo o chama, mas a checagem de bula da
> etapa 2 faz o mesmo e mais: decide sobre o **acervo inteiro** do openFDA
> (265.009 bulas) em vez das 50 mais recentes, inclui as seções
> `contraindications` e `precautions`, registra as duas contagens que dão peso à
> evidência (118 de 166 bulas contra 1 de 211), e sai com veredito
> `indeterminado` — que nunca remove — quando não há bula, em vez de tratar
> ausência de dado como ausência de risco. Para trabalho novo, use a etapa 2.

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

**Fluxo em dois passos.** A coleta de evidência **não** usa os alvos enriquecidos: a lista precisa ser revisada antes de gastar a previsão sobre todos os ligantes (210 alvos x ~109 mil ligantes são ~23 milhões de pares). A execução avisa isso ao final, e o manifesto registra `alvos_enriquecidos_usados_na_evidencia` **a partir da procedência dos alvos que a evidência realmente usou** (`origem_alvo`), não a partir de qual flag foi passada: no passo 2, em que `--alvos-enriquecidos` nem aparece, o valor é `true`, e `origem_dos_alvos_da_evidencia` mostra a contagem por origem. Passar `--alvos` e `--alvos-enriquecidos` apontando para o **mesmo arquivo** é recusado — a execução leria a lista revisada e a reescreveria por cima da sua revisão.

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

## 11. Etapa 1 completa — painel, linha de comando e como reproduzir

A etapa 1 transforma três arquivos (ligantes, alvos, sintomas) numa tabela de evidências, `evidencias.tsv`, e num `manifesto.json`. Há duas portas de entrada para a mesma lógica: o **painel** e a **linha de comando**. A CLI é a interface que o artigo cita; o painel apenas a dispara, e nada de coleta mora nele.

### Pré-requisitos

```bash
pip install -r requirements.txt              # streamlit, rdkit, pandas, requests, scipy
pip install -r requirements-opcionais.txt    # torch -- só para evidência prevista
```

- `target_labels.tsv` do PP4CKD (**obrigatório**: diz quais alvos existem no modelo e quantas moléculas treinaram cada um);
- pesos `.pt` do PP4CKD (opcional: sem eles, e sem `torch`, não há evidência prevista);
- `evidencias_chembl.csv` (opcional: evidência medida).

> **O `.pt` e o `target_labels.tsv` têm de ser do mesmo PP4CKD.** A coleta
> **recusa** um par cuja largura de saída do checkpoint não seja o número de
> rótulos, nomeando as duas: a coluna *i* de um não é o alvo *i* do outro, e
> pontuar assim atribuiria a probabilidade de um alvo ao nome de outro. Como o
> tarball publicado está truncado e este README manda trocar o checkpoint por
> `ecfp4_fold_4.pt`, o par divergente é o caso esperado — troque os dois juntos.

### Como as entradas são lidas (regras que mudam o resultado)

- **TSV não tem citação.** Aspas fazem parte do nome: `"cis" isomer of X` é lido
  com as aspas. Com o dialeto padrão de CSV ele era lido como
  `cis isomer of X` — o rótulo de identidade reescrito em silêncio. Em `.csv` a
  citação continua valendo, porque é dela que a vírgula do nome depende.
- **`alvo_chembl_id` é normalizado para maiúsculas na leitura**, num lugar só.
  Antes, `chembl222` numa TSV dava `no_modelo = nao` no arquivo enriquecido e
  nenhuma semente do Reactome, enquanto a previsão pontuava o alvo — as duas
  saídas se contradiziam sobre o mesmo alvo.
- **Alvo repetido é dobrado e relatado.** A identidade é o ID, ou
  nome+organismo quando não há ID. Numa tabela de uma linha por evidência, o
  mesmo alvo pedido duas vezes (por ID numa linha e por nome noutra) dobrava
  cada ensaio medido dele, e a linha repetida era indistinguível de um segundo
  ensaio. Ligantes já eram deduplicados pela forma parental.
- **Todo alvo sai com `origem_alvo`** (`usuario` por padrão). A coluna `origem`
  do `alvos_enriquecidos.tsv` é aceita com esse nome, então realimentar aquele
  arquivo preserva a procedência de cada alvo em vez de perdê-la.
- **Alvo por nome casa com organismo EXATO.** `Homo` ou `sapiens` não casam
  mais `Homo sapiens` por trecho, e um nome que casa com mais de um alvo do
  modelo é **ambíguo**: nenhum é escolhido, e o caso é relatado (o painel
  mostra a lista). Escolher o primeiro era a última porta pela qual um
  ortólogo entrava em silêncio.

### Códigos de saída da CLI

| código | significado |
|---|---|
| 0 | execução completa |
| 3 | execução **parcial**: a saída existe, mas linhas de entrada foram descartadas ou houve falha de item. O manifesto registra `execucao_parcial: true` e `falhas_de_item` |
| 1 | nenhuma entrada utilizável (nada é escrito) |
| 2 | erro (nada é publicado; as saídas anteriores ficam intactas) |

O código 3 existe porque automação não lê o `stderr`: antes, uma coleta que
descartou linhas saía com 0 e era indistinguível de uma completa.

### Painel

```bash
export PP4CKD_ROTULOS=/caminho/target_labels.tsv     # padrão da barra lateral
export PP4CKD_PESOS=/caminho/ecfp4_fold_4.pt         # opcional
export EVIDENCIAS_CHEMBL=/caminho/evidencias_chembl.csv   # opcional
streamlit run painel.py
```

1. **Entradas.** Envie ligantes, alvos e sintomas. A validação roda no upload, antes de qualquer execução, e mostra o que vai dar errado: alvo por nome sem organismo, SMILES ilegível, alvo fora do modelo, alvo com menos de 50 moléculas de treino.
2. **Enriquecer alvos** (opcional, botão próprio). "Buscar EFO dos sintomas" lista as alternativas do Open Targets e **você escolhe**; sem escolha o sintoma não é enriquecido. A lista proposta aparece editável. **Os alvos enriquecidos não são usados na coleta até você clicar em "Usar esta lista na coleta".** Os botões são separados de propósito: uma execução real deu 235 alvos contra 108.935 ligantes, cerca de 23 milhões de pares, e o painel mostra essa conta antes da coleta.
3. **Coletar evidências.** Roda em segundo plano (um enriquecimento sem cache leva 23 a 25 min); o painel mostra o tempo decorrido e o relatório que a CLI escreve. A CLI não informa percentual, então a barra de tempo é estimativa, não progresso. O resultado traz tabela filtrável e botões para baixar `evidencias.tsv` e `manifesto.json`.

   **O resultado é lido da pasta nomeada pelas entradas que estão na tela** —
   os três arquivos enviados, a lista de alvos usada, o corte, o limite e as
   fontes de evidência —, do mesmo jeito que a lista enriquecida. Resultado de
   outra execução nunca aparece sob "Resultado" como se fosse o atual, e uma
   coleta cancelada ou que falhou não mostra a tabela da anterior (a execução
   nova limpa a saída daquela chave antes de começar). O painel ainda diz de
   que execução o resultado é (data do manifesto), avisa quando o manifesto
   está ausente e avisa quando a execução foi **parcial**. A tabela da tela
   mostra no máximo 5.000 linhas, lidas em fluxo; o arquivo inteiro vem pelo
   botão de download, que só abre o arquivo quando você clica.

**Modo degradado.** Sem os pesos, o painel roda só com evidência medida e diz isso na tela; sem pesos nem `evidencias_chembl.csv` ele só enriquece. O enriquecimento não usa pesos.

O painel não tem modo de repetição: para repetir uma execução use a CLI (abaixo).

### Linha de comando

Um passo (coleta sobre os alvos informados):

```bash
python mods/coletar_evidencias.py --ligantes ligantes.tsv --alvos alvos.tsv --sintomas sintomas.txt \
    --saida evidencias.tsv --checkpoint ecfp4_fold_4.pt --rotulos target_labels.tsv \
    --evidencias-chembl evidencias_chembl.csv
```

Dois passos, com enriquecimento (a revisão da lista fica entre eles):

```bash
# 1. enriquecer; evidencias.tsv deste passo usa so os alvos informados
python mods/coletar_evidencias.py --ligantes ligantes.tsv --alvos alvos.tsv --sintomas sintomas.txt \
    --saida passo1/evidencias.tsv --rotulos target_labels.tsv \
    --alvos-enriquecidos alvos_enriquecidos.tsv --efo "neuropathic pain=EFO_0005762"
# 2. revise alvos_enriquecidos.tsv e colete com ele
python mods/coletar_evidencias.py --ligantes ligantes.tsv --alvos alvos_enriquecidos.tsv \
    --sintomas sintomas.txt --saida evidencias.tsv --checkpoint ecfp4_fold_4.pt \
    --rotulos target_labels.tsv --evidencias-chembl evidencias_chembl.csv
```

Também se pode fazer tudo numa chamada só, passando `--alvos-enriquecidos` junto com `--checkpoint`: a execução grava a lista enriquecida e avisa que a coleta **não** a usou (o manifesto registra `alvos_enriquecidos_usados_na_evidencia: false`, porque todos os alvos lidos têm `origem_alvo = usuario`).

**As três saídas são publicadas juntas.** Tabela, lista enriquecida e manifesto
são escritas em temporários e só então trocadas pelos arquivos de destino. Uma
falha no meio não publica nada e deixa intactos o `evidencias.tsv`, o
`alvos_enriquecidos.tsv` revisado à mão e o `manifesto.json` anteriores — e
nunca sobra tabela nova sem manifesto, que é o que a torna citável.

### As 26 colunas de `evidencias.tsv`

Uma linha é **uma evidência** sobre um par ligante/alvo — medida ou prevista.
Um par com três ensaios medidos e uma previsão são quatro linhas. É esta tabela
que uma seção de métodos cita, então cada coluna está definida aqui.

**A regra que não pode ser ignorada:** `pchembl` é a **única** coluna
comparável entre linhas. `valor` + `unidade` **nunca** se comparam entre
unidades (120 nM e 78% não se ordenam), e `valor` de uma linha prevista é
probabilidade, não potência. Ordenar a tabela por `valor` mistura grandezas
diferentes; ordenar por `pchembl` (filtrando `tipo_evidencia = medido`) não.

*Identidade*

| coluna | o que é |
|---|---|
| `ligante` | nome como veio da entrada |
| `ligante_inchikey_parental` | InChIKey da forma parental (sal e base livre são o mesmo composto). É a chave de junção |
| `ligante_chembl_id` | ID ChEMBL do ligante, quando há. Vazio quando a correspondência não é de identidade |
| `correspondencia_ligante` | `exata`, `estereoisomero` ou `desconhecida`. **`estereoisomero` é outro composto**: a atividade medida nele não é a do ligante. A linha fica (a fase 2 decide), mas filtrável sem ler prosa. Linha prevista deixa vazio |
| `alvo_chembl_id`, `alvo`, `organismo`, `tipo_alvo` | o alvo. `tipo_alvo` `SINGLE PROTEIN` é ligação a uma proteína; `ORGANISM` e `CELL-LINE` são ensaio fenotípico — afirmações diferentes |
| `origem_alvo` | de onde veio o alvo: `usuario`, `open_targets` ou `reactome`. É o que permite remover a expansão inteira filtrando uma coluna |

*Evidência*

| coluna | o que é |
|---|---|
| `tipo_evidencia` | `medido` ou `predito`. **Nunca se misturam numa célula** |
| `fonte` | `ChEMBL` (medido) ou `PP4CKD` (predito) |
| `atividade_relatada` | texto de exibição (`IC50 = 120 nM`, `p = 0.93 (ECFP4)`). Para ler, não para calcular |
| `medida` | `IC50`, `Ki`, `Potency`… no medido; `probabilidade` no predito |
| `relacao` | `=`, `>`, `<`. Um `>` é **limite**, não medida: entra marcado |
| `valor`, `unidade` | número e unidade como a fonte os traz (ponto decimal). **Não comparáveis entre unidades** |
| `pchembl` | −log10 da potência molar, **a única coluna comparável entre linhas**. Vazia quando a fonte não a traz; nunca calculada nem preenchida com outra coisa. Sempre vazia no predito |

*Lastro (por que acreditar na linha)*

| coluna | o que é |
|---|---|
| `tipo_ensaio` | código do ChEMBL: `B` ligação, `F` funcional, `A` ADME |
| `ensaio_descricao` | descrição do ensaio, como o ChEMBL a publica |
| `triagem_hts` | `sim`/`nao`. Varredura HTS de painel em vez de ensaio dirigido (heurística `PUBCHEM_BIOASSAY`/`qHTS`). Coluna de primeira classe porque foi medida: recall@10 0,111 em HTS contra 0,272 em ensaio dirigido — **somar os dois estraga a leitura** |
| `confianca_ensaio` | `confidence_score` do ensaio. Vazia nesta etapa: o `evidencias_chembl.csv` não a traz |
| `n_moleculas_treino` | moléculas ativas que treinaram aquele alvo (só no predito). p = 0,9 com 8 exemplos não vale p = 0,9 com 4.000 |
| `documento`, `pmid`, `ano` | o documento de onde o ensaio veio. Linha sem rastro até o documento não serve para publicação |
| `procedencia` | carimbo completo numa string: `medido (ChEMBL, evidencias_chembl.csv)` ou `predito (PP4CKD ECFP4 ecfp4_fold_4.pt)`, com o checkpoint nomeado. Marca também `estereoisomero da molecula ensaiada <ID>` |

**Formato.** TSV, cabeçalho iniciado por `#`, ponto decimal em toda coluna de
máquina (a vírgula só aparece no texto de `atividade_relatada`), e escape
reversível estilo TSV do MySQL (`\\` `\t` `\n` `\r`) — então tabulação e quebra
de linha dentro de uma descrição de ensaio não quebram a tabela.

### O conjunto-controle versionado: `dados/evidencias_controle.tsv`

49 linhas de evidência **medida** sobre os 12 fármacos do conjunto-controle de
dor neuropática (`exemplos/controle.dor_neuropatica.tsv`) e os 9 alvos de
`exemplos/alvos.dor_neuropatica.txt`. Está no repositório para que os critérios
de aceite da etapa 1 sejam verificáveis **sem rede e sem o
`evidencias_chembl.csv` de 60 MB** — é o que `tests/test_controle.py` afirma
(critério 2: a amitriptilina tem linhas `medido` para o transportador de
noradrenalina, CHEMBL222, com documento e `pchembl`).

Foi gerado pela própria CLI documentada acima:

```bash
python mods/coletar_evidencias.py \
    --ligantes exemplos/controle.dor_neuropatica.tsv \
    --alvos exemplos/alvos.dor_neuropatica.txt \
    --sintomas exemplos/sintomas.dor_neuropatica.txt \
    --saida dados/evidencias_controle.tsv \
    --rotulos ~/modelos_baixados/pp4ckd/target_labels.tsv \
    --evidencias-chembl ~/prospector/pipeline_prospeccao_v4/resultados/evidencias_chembl.csv \
    --sem-rede
```

A fonte é o `evidencias_chembl.csv` do pipeline de prospecção (bioatividade
baixada do ChEMBL), e **não** a API ao vivo. Cobertura medida nessa fonte:
amitriptilina 175 linhas, ciprofloxacino 179, capsaicina 105, gabapentina 22,
pregabalina 15; duloxetina, carbamazepina e linezolida têm **zero** — não estão
no arquivo, que cobre os compostos garimpados pela prospecção. É por isso que o
critério de aceite 2 do desenho fala de amitriptilina, e não de duloxetina.
Só 49 das linhas da fonte casam com os 9 alvos pedidos.

> **Licença:** por derivar do ChEMBL, este arquivo é **CC BY-SA 3.0**
> (`dados/LICENSE`), não MIT. O mesmo vale para as fixtures derivadas do ChEMBL
> em `tests/dados` (`tests/dados/LICENSE`).

### Pesos do PP4CKD: o tarball do Zenodo está truncado na origem

O `pp4ckd_models.tar.gz` publicado no Zenodo tem MD5 igual ao checksum publicado (`f9380bb9cef3fe56d6d49d3e4a4d57ec`), mas o tar termina antes do fim: sobram íntegros 22 dos ~77 modelos, e `ecfp4_full_model.pt` não está entre eles. Não é download corrompido. **Contorno:** use `ecfp4_fold_4.pt`, um modelo real, treinado em 9/10 dos dados (detalhes na seção 8).

### O que o manifesto registra

Além dos parâmetros, das versões dos serviços e dos MD5 das entradas (checkpoint,
rótulos, evidência medida, ligantes, alvos, sintomas):

- `pipeline_commit` — o commit que rodou, com sufixo **`+sujo`** quando a árvore
  de trabalho tem modificação não commitada. Um hash limpo de código modificado
  seria uma afirmação falsa de reprodutibilidade. `+estado_desconhecido` quando
  o `git status` não pôde ser lido, `desconhecido` sem git.
- `pipeline_repositorio` — de qual repositório o commit foi lido, com um aviso
  explícito quando essa raiz não é a do pacote (pacote embutido em outro
  repositório: o hash não identificaria o código desta etapa).
- `python` e `rdkit` — o RDKit decide a química (InChIKey parental, fingerprint),
  então a versão dele muda o resultado. Caminhos não são resolvidos: os MD5
  prendem o conteúdo, que é mais forte.
- `execucao_parcial` e `falhas_de_item` (por fonte: `previsto`,
  `enriquecimento`, `medido`) — é o que impede uma coleta incompleta de parecer
  completa. Uma repetição compara as falhas de item com as da execução original
  e avisa quando diferem.
- `origem_dos_alvos_da_evidencia` e `alvos_enriquecidos_usados_na_evidencia` —
  derivados da procedência dos alvos lidos, não das flags passadas.

### Cache e repetição

Toda resposta de rede (ChEMBL, Open Targets, Reactome) é gravada em `--cache`, por versão dos serviços. `--manifesto manifesto.json` repete uma execução sem rede, servindo tudo do cache, e se recusa a rodar se uma versão, um MD5 de entrada ou um parâmetro mudou (seção 10). Duas execuções com o mesmo manifesto produzem o mesmo arquivo.

### Licenças

O **código** é MIT (`LICENSE`). Os **dados derivados do ChEMBL** são
**CC BY-SA 3.0**, e isso vale onde eles estiverem versionados: `dados/`
(`dados/LICENSE`) e as fixtures de teste derivadas do ChEMBL
(`tests/dados/LICENSE`, que nomeia os caminhos: `evidencias_amostra.csv` e as
respostas gravadas do ChEMBL). A cláusula share-alike vale só para esses dados
e não contamina o código. As fixtures do Reactome (só anotação) e do Open
Targets são CC0 e não exigem nada.

### Testes

```bash
python -m pytest tests/        # nenhum teste usa a rede
```

A suíte precisa de `requirements-opcionais.txt` (traz o `torch`):

```bash
pip install -r requirements-opcionais.txt
python -m pytest tests/
```

Sem o `torch` o `pytest` **se recusa a rodar** (código de saída 2, com a
mensagem dizendo o que instalar), em vez de pular os testes da evidência
prevista, incluindo o guarda de largura do checkpoint, e relatar sucesso.

Cobertura do painel, dita sem enfeite:

- A lógica que ele exibe (`mods/evidencias/resumo.py`) tem testes unitários (`tests/test_resumo.py`).
- A casca Streamlit **sobe in-process** pelo `AppTest` do Streamlit em `tests/test_painel.py` (sem navegador e sem rede): erro quando falta o `target_labels.tsv`, aviso de modo degradado sem pesos, avisos de validação no upload (com o `file_uploader` substituído pelos arquivos de `tests/dados`), a lista de alvos enriquecidos ligada aos arquivos e EFO que a geraram, o resultado ligado às entradas da tela (e o de outra execução *não* aparecendo), o portão de redesenho e o descarte de edições.
- **Não verificados:** navegador real, o widget de upload real, uma execução ao vivo (rede, pesos reais), o botão **Cancelar** e o `killpg` que ele usa (nenhum teste mata um processo de verdade, e não há fallback para `SIGKILL`), os **downloads diferidos** (a função só é chamada quando o runtime do Streamlit existe, o que o `AppTest` não exercita), e a tabela de resultado com dados reais de milhões de linhas.
  - O que **era** não verificado e agora é: o portão de redesenho (o aviso e o botão de descarte só aparecem quando um redesenho foi de fato suprimido) e o descarte das edições. O descarte deixou de depender de remover o estado do widget do `session_state` — que o Streamlit reenvia no rerun e que nunca foi verificado — e passou a criar uma **geração nova** do editor: outro widget, com estado vazio por construção.

---

## 12. Etapa 2 — o portão (`mods/aplicar_portao.py`)

A etapa 1 coleta e **nunca decide**. A etapa 2 é onde se decide — e o que ela
decide, medido, é muito menos do que o nome sugere.

### O nome descreve a posição, não a severidade

O portão quase não reprova, e isso não é timidez: é o que a medição manda. No
conjunto-controle de 12 fármacos (seis que tratam dor neuropática, seis que a
causam), cada critério plausível de exclusão apaga os próprios positivos de
referência do pipeline:

- reprovar por **PAINS** apagaria a **amitriptilina**;
- por **BRENK**, **65% da base** de 108.935 ligantes;
- por **expressão** abaixo de 1 TPM em tecido nervoso, **SERT e NET** — os
  alvos da duloxetina e da amitriptilina.

Com `--reprovar alerta_estrutural`, o controle de 12 perde 8 — entre eles
amitriptilina, carbamazepina e capsaicina, que **tratam** — e continua
aprovando o ciprofloxacino, que **causa**. Promover o critério inverte o
pipeline. Por isso o resumo diz o preço **antes**: para cada critério
rebaixado ele imprime quantos ligantes "reprovaria se promovido".

Então o portão **marca, anota e explica**. Por padrão só duas coisas reprovam:
SMILES que o RDKit não lê (não há molécula para avaliar) e alerta forte de
bula (0 de 6 falsos positivos medidos no controle).

### Os cinco critérios e a política padrão

| critério | o que olha | padrão | o número que decidiu o padrão |
|---|---|---|---|
| `smiles` | o RDKit lê a molécula | **reprovar** | sem molécula não há o que avaliar |
| `alerta_bula` | termos da condição na bula (openFDA) | **reprovar** | 0/6 de falso positivo e 4/6 pegos no controle |
| `alerta_estrutural` | PAINS, BRENK, NIH (RDKit `FilterCatalog`) | marcar | PAINS reprovaria a amitriptilina; BRENK, 65% da base |
| `toxicidade` | ADMET-AI, percentil ≥ `--percentil-toxicidade` | marcar | percentil não validado neste conjunto |
| `expressao` | expressão basal do alvo nos tecidos pedidos | marcar | reprovaria SERT (0,454 TPM) e NET |

`expressao` é critério de **alvo**, não de ligante: ele vive em
`alvos_anotados.tsv` e não entra em `portao.tsv`, porque um alvo de expressão
baixa não reprova um ligante — o ligante pode ter outros alvos.

`--reprovar C[,C]` promove critérios; `--marcar C[,C]` rebaixa. Nome que não
existe é **erro** (código 2), nunca flag ignorada em silêncio: uma política que
o usuário acredita estar valendo e não está é pior que uma mensagem de erro.
A coluna `politica` de cada linha diz sob qual regra aquela linha foi julgada.

### A semântica dos vereditos, uma só para todos os critérios

Cada critério devolve um veredito **intrínseco**, independente de política:

| veredito | significa |
|---|---|
| `passa` | o critério olhou e não achou problema |
| `marca` | achou algo que merece ser visto |
| `reprova` | achou algo que o próprio módulo considera impeditivo |
| `indeterminado` | **não havia dado para julgar** |

E três regras atravessam a etapa inteira:

1. **`indeterminado` NUNCA tira um ligante**, com política nenhuma. Ausência de
   dado não é evidência de problema — é o caso de quase todo produto natural na
   checagem de bula e de **todo alvo não humano** na expressão.
2. **Critério que a execução não pediu sai `indeterminado`, nunca `passa`.** Sem
   `--termos-alerta` a checagem de bula não roda, e "não perguntei" não é a
   mesma afirmação que "sem alerta". O mesmo vale para o que passou do
   `--limite-bula`: sai `indeterminado` com esse motivo escrito.
3. **Nada desaparece em silêncio.** Reprovado continua na tabela, com o
   critério e o valor que decidiu em `motivos`. Ligante cujo SMILES o RDKit não
   lê também: ele é **julgado**, não descartado na leitura.

### As quatro saídas (`--saida` é um diretório)

| arquivo | uma linha por | colunas |
|---|---|---|
| `propriedades.tsv` | ligante | 120: as 16 de identidade/descritores + as 104 do ADMET-AI |
| `portao.tsv` | ligante | 13: identidade, `aprovado`, `motivos`, `<critério>_veredito`/`<critério>_valor` × 4, `politica` |
| `alertas_bula.tsv` | (ligante, termo, seção) | 10: termo, seção, força, trecho literal, contagens, `campo_do_nome`, procedência |
| `alvos_anotados.tsv` | alvo | 30: as 18 da lista enriquecida da etapa 1 + 9 de expressão + `aprovado`, `motivos`, `politica` |

Mais o `manifesto.json`. **As cinco saídas são publicadas juntas**: cada uma
nasce num temporário no mesmo diretório e só toma o lugar do destino quando
todas estão prontas. Tabela sem manifesto não é citável.

Três detalhes das colunas que mudam leitura:

- **`alerta_pains`, `alerta_brenk` e `alerta_nih` trazem os nomes dos filtros**
  (`nitro_group`, não `1`). Um contador não diz o que casou, e é o nome que
  permite julgar se o alerta importa para aquele composto.
- **O ADMET-AI publica uma coluna `tpsa`, que colidia com a TPSA do RDKit.** As
  colunas do ADMET entram com os nomes originais dele, exceto as que colidem
  com o esquema deste repositório, que recebem o prefixo `admet_` —
  `admet_tpsa`. O manifesto registra cada renomeação. Antes disso o cabeçalho
  tinha `tpsa` duas vezes e o valor do RDKit era substituído em silêncio.
- **`toxicidade_valor` é percentil, não probabilidade crua.** 0,67 de DILI não
  quer dizer nada sozinho: a aspirina tem 0,67 e está aprovada. O percentil é
  contra os aprovados do DrugBank que o ADMET-AI traz.

### Por que a bula decide no acervo e cita na janela

A checagem de bula faz **três requisições por ligante**:

1. uma **janela de 50 bulas ordenada por data** (`sort=effective_time:desc`),
   que fornece o **trecho** citado em `alertas_bula.tsv` — da rotulagem
   vigente, não de uma bula sorteada;
2. duas **consultas dirigidas** (seção forte e indicação) que rodam no **acervo
   inteiro**, 265.009 bulas, e é delas que vem a **decisão**.

A janela sozinha não bastava, e isso foi medido: a capsaicina tem 211 bulas e a
única que a indica para dor neuropática (Qutenza) estava fora das 50 — o
desconto por indicação falhava em silêncio e o falso positivo ia de 0/6 para
1/6. Medido também o preço: 118 das 166 bulas do ciprofloxacino trazem
`neuropathy` em seção forte, e a resposta dirigida custa 202 KB contra 27 MB de
baixar as 166.

Daí as duas contagens em `alertas_bula.tsv`: `n_bulas_consultadas` é o tamanho
da janela que forneceu o trecho, e `n_bulas_com_o_termo` é quantas bulas do
princípio ativo casam o termo no acervo. A diferença é o peso da evidência —
**118 de 166** no ciprofloxacino contra **1 de 211** na indicação da capsaicina.
Uma contagem feita dentro da janela não distingue esses dois casos.

As seções fortes são `boxed_warning`, `warnings`, `warnings_and_cautions`,
`contraindications` e `precautions`. `adverse_reactions` é registrada mas não
decide, e `indications_and_usage` **desconta**: uma droga que trata o sintoma
não é confundida com uma que o causa.

### A GPU é mais rápida e não é reprodutível

O padrão é **CPU**. Medido entre duas execuções idênticas: na GPU, 36 de 312
colunas diferem (última casa, ~1e-8 relativo, que o ranqueamento por percentil
amplifica para ~0,035); na CPU, 0 de 312. Nenhum veredito mudou nessa medição,
mas um manifesto que promete repetição byte a byte e não entrega é pior que um
que avisa. O preço da CPU é ~66 min para a base inteira contra ~37 min na GPU.

`--gpu` avisa de saída, o manifesto registra **qual dispositivo realmente
rodou** (não o pedido), e repetir uma execução de GPU avisa em vez de prometer
igualdade bit a bit.

E `--workers 0` é o padrão porque é **45× mais rápido**: 1.600 moléculas saem
em 8,6 s com 0 e em 384 s com 8 — são dez modelos pagando o nascimento de N
processos cada.

### Comando completo

```bash
python mods/aplicar_portao.py \
    --ligantes ligantes.tsv \
    --evidencias evidencias.tsv \
    --alvos alvos_enriquecidos.tsv \
    --tecidos tecidos.txt \
    --termos-alerta "neuropathy,neuropathic,neuritis,paresthesia" \
    --saida portao/
```

- `--evidencias` restringe o portão aos ligantes **que têm** evidência da etapa
  1, e o resumo diz quantos ficaram de fora; sem ele, toda a entrada é julgada.
- `--tecidos` recebe ID UBERON ou nome **exato** do tecido, um por linha
  (`#` comenta). O casamento é por ID e por ancestralidade na UBERON, nunca por
  trecho de nome: `cortex of kidney` não casa com `cortex`. Sem tecidos a
  expressão é anotada mas o veredito fica `indeterminado`, porque não há com o
  que comparar.
- `--limite-bula` (padrão 500) é o teto de ligantes consultados no openFDA, que
  é serviço público com limite de taxa. Para a base inteira, 108.935 ligantes ×
  3 requisições a 240/min dariam ~23 h; o caminho é o download em massa
  (1,7 GB em 14 partições), **ainda não implementado**.
- `--sem-rede` pula bula e expressão; RDKit e ADMET-AI continuam.

Códigos de saída são os da etapa 1: `0` completa, `3` **parcial** (a saída
existe mas linhas foram descartadas; o manifesto diz isso), `1` nenhum ligante
utilizável, `2` erro (nada publicado).

### Cache e repetição

Toda resposta de rede (openFDA, Open Targets) é gravada em `--cache` (padrão:
`cache_rede/` dentro de `--saida`), em subpasta por versão dos serviços.
`--manifesto manifesto.json` repete a execução servindo tudo do cache, sem
rede — **falta no cache é erro**, não uma requisição nova.

A repetição **recusa rodar** se mudou um MD5 de entrada (ligantes, alvos,
evidências, tecidos), uma versão (`openfda`, `admet_ai`, `rdkit`,
`open_targets`) ou um parâmetro que muda a saída — inclusive a política e o
dispositivo. Caminhos de saída não entram nessas chaves: repetir para outra
pasta é repetir. O manifesto novo sai como `manifesto_repeticao.json` e nunca
sobrescreve o original.

### Testes

```bash
python -m pytest tests/test_portao.py tests/test_portao_cli.py
```

Nenhum usa rede: as respostas do openFDA e do Open Targets são fixtures
versionadas em `tests/dados/respostas/` (aparadas; o método e as consequências
estão em `tests/dados/LICENSE` — elas respondem fielmente só para os quatro
termos gravados). Os testes prendem as medições que decidiram os padrões: PAINS
pegando a amitriptilina, `indeterminado` não reprovando sob política nenhuma,
0/6 de falso positivo e exatamente `ciprofloxacin, linezolid, ethambutol,
metronidazole` pegos no controle, SERT e NET marcados e não reprovados, alvo
não humano saindo `indeterminado`, exatamente 3 requisições por ligante, e a
repetição byte a byte.

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
- **A checagem de alertas adversos não entende negação em texto livre** — em
  `checar_alertas_adversos.py` e também no critério `alerta_bula` da etapa 2:
  uma bula que diz "não há relatos de neuropatia" ainda casa o termo. Trechos
  sinalizados precisam de revisão humana, e é por isso que ambos expõem o
  trecho literal em vez de só um sinalizador. A etapa 2 mitiga um caso
  específico — a droga que *trata* o sintoma — descontando o que aparece em
  `indications_and_usage`, mas não resolve negação em geral.
- **A validação por decoys não testa especificidade entre classes farmacológicas** — mostra que drogas de referência conhecidas superam controles negativos aleatórios, mas não testa diretamente se compostos de uma classe não relacionada estariam sendo indevidamente favorecidos. A checagem de alertas adversos (etapa 6) é um primeiro passo nessa direção, mas específico a efeitos já documentados em bula, não uma solução geral.
- **Calibração do Deep-PK e comparação formal com outras ferramentas ADMET continuam sem execução.** Na etapa 2 a questão mudou de forma, não desapareceu: o ADMET-AI 2.0.1 substituiu o Deep-PK e é local, versionado e reprodutível, mas o **percentil de toxicidade não está validado neste conjunto** — é por isso que `toxicidade` só marca, nunca reprova, por padrão.
- **A checagem de bula da etapa 2 não cobre a base inteira numa execução.** O teto é `--limite-bula` (padrão 500) porque o openFDA é serviço público com limite de taxa: 108.935 ligantes × 3 requisições a 240/min dariam ~23 h. O caminho certo é o download em massa (1,7 GB em 14 partições), **registrado e não implementado**. O que passa do teto sai `indeterminado`, nunca "sem alerta".
- **A fase 3 (ranqueamento que usa os alvos) não existe ainda.** As etapas 1 e 2 entregam evidência e vereditos com procedência; quem ordena usando isso está por escrever. O requisito está em `docs/superpowers/specs/2026-11-25-fase3-alvos-no-ranqueamento-requisitos.md`.
