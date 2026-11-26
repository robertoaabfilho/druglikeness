# Etapa 2 — O portão (plano de execução)

Desenho: `../specs/2026-11-25-etapa2-portao-design.md`. Base já commitada:
`c9ebb78` (procedência compartilhada, openFDA no cache, chaves de repetição
parametrizáveis, `ensembl`/`escore_rna_expression` na lista enriquecida, e os
seis módulos de `mods/portao/`).

O que falta é a **montagem**: as peças existem e nada as costura. Este plano
entrega `mods/aplicar_portao.py`, os testes e a documentação.

## Estado de partida

| peça | situação |
|---|---|
| `mods/portao/esquema.py` | pronto (`EscritorTSV`, `colunas_portao`, `ler_tsv`) |
| `mods/portao/propriedades.py` | pronto (`linha_de_propriedades`) |
| `mods/portao/admet.py` | pronto (`disponivel`, `colunas`, `prever`, `veredito_de_toxicidade`) |
| `mods/portao/bula.py` | pronto (`buscar_bulas`, `achados`, `veredito`, `versao_openfda`) |
| `mods/portao/expressao.py` | pronto (`ler_tecidos`, `expressao_basal`, `anotar`) |
| `mods/portao/decisao.py` | pronto (`montar_politica`, `linha_de_portao`, `Contador`) |
| **CLI** | **não existe** |
| **testes** | **nenhum** para `mods/portao/` |
| README | não documenta a etapa 2 |

## Interface da linha de comando

`python mods/aplicar_portao.py`:

| flag | padrão | o que faz |
|---|---|---|
| `--ligantes` | obrigatório | TSV da etapa 1, lido por `entrada.ler_ligantes` |
| `--saida` | obrigatório | **diretório**: recebe as 4 tabelas e `manifesto.json` |
| `--evidencias` | — | `evidencias.tsv`; restringe o portão a quem tem evidência |
| `--alvos` | — | `alvos_enriquecidos.tsv`; liga a anotação de expressão |
| `--tecidos` | — | arquivo de tecidos (ID UBERON ou nome exato), um por linha |
| `--termos-alerta` | — | repetível; sem nenhum, a checagem de bula não roda |
| `--limite-bula` | 500 | teto de ligantes consultados no openFDA; acima dele, `indeterminado` com esse motivo |
| `--percentil-toxicidade` | 90 | percentil contra os aprovados do DrugBank |
| `--corte-expressao` | 1.0 | abaixo disso nos tecidos pedidos, `expressao` marca |
| `--workers` | 0 | `num_workers` do ADMET-AI (0 é 45× mais rápido; medido) |
| `--reprovar` / `--marcar` | — | promove/rebaixa critérios (`c[,c]`) |
| `--sem-rede` | — | pula bula e expressão (ADMET e RDKit continuam) |
| `--cache` | `cache_rede/` ao lado da saída | raiz do cache, subpasta por versão dos serviços |
| `--manifesto` | — | modo de repetição |
| `--limite` | — | só os N primeiros ligantes |

Códigos de saída: os da etapa 1 (`0` completa, `3` parcial, `1` nenhuma
entrada utilizável, `2` erro — nada publicado).

## Tarefas

### 1. Esqueleto do CLI e publicação atômica

`mods/aplicar_portao.py` com a mesma forma de `coletar_evidencias.py`: toda
validação antes de abrir saída; cada tabela nasce em temporário no mesmo
diretório (`atomico.temporario`) e **as cinco saídas são publicadas juntas** —
tabela sem manifesto não é citável. O cache abre com
`cache.abrir_ao_vivo(raiz, com_openfda=True)` só quando a execução usa rede
(bula ou expressão), e `abrir_para_repeticao` no modo manifesto.

Validações que param antes de qualquer escrita:
- `--saida` que aponte para um arquivo existente, ou para o diretório de uma
  entrada;
- `--tecidos` sem `--alvos` (expressão não tem o que anotar);
- `--termos-alerta` com `--sem-rede` (bula precisa de rede);
- critério desconhecido em `--reprovar`/`--marcar` (`CriterioDesconhecido`
  vira mensagem e código 2, não traceback).

### 2. `propriedades.tsv` e o critério `smiles`

Streaming: `EscritorTSV(COLUNAS_PROPRIEDADES + colunas_do_admet)`. Duas
passagens sobre os ligantes não cabem (108.935 × 120 colunas), então:

1. `admet.prever` por lote de `LOTE` SMILES, acumulando só o lote corrente;
2. para cada ligante do lote, `linha_de_propriedades` + a linha do ADMET.

Ligante com SMILES ilegível **fica na saída** com as colunas vazias e
`smiles_valor = "SMILES que o RDKit nao le"`; é a única reprovação estrutural
padrão. SMILES que o ADMET-AI descartou, ou ADMET-AI ausente, dá
`toxicidade = indeterminado`.

### 3. `alertas_bula.tsv` e o critério `alerta_bula`

Só com `--termos-alerta`. Por ligante, em ordem determinística, até
`--limite-bula`: `buscar_bulas` → `achados` → `veredito`. Uma linha por achado
(termo, seção), com `n_bulas_consultadas` e `campo_do_nome`. Acima do teto,
veredito `indeterminado` com `"acima do teto de bula (--limite-bula N)"`.
`procedencia = "openFDA <last_updated>"`, da versão que o cache já consultou.

### 4. `alvos_anotados.tsv` e o critério `expressao`

Só com `--alvos`. Por alvo: `ensembl` da própria coluna (etapa 1 já grava);
sem ela, `ensembl_do_gene`. Depois `expressao_basal` + `anotar(pedidos,
corte=--corte-expressao)`. Colunas originais de `alvos_enriquecidos.tsv` mais
`COLUNAS_ALVOS_EXTRA`. Alvo sem GTEx (não humano) é `indeterminado`, nunca
expressão zero.

**O critério `expressao` é de alvo, não de ligante** — `portao.tsv` é uma
linha por ligante. Então ele não entra em `portao.tsv`: entra em
`alvos_anotados.tsv` e no resumo, e o `Contador` o conta sobre os alvos. A
política vale igual (`--reprovar expressao` tira o alvo da lista anotada, com
o motivo na linha). Registrar isso no manifesto como `escopo: alvo`.

### 5. `portao.tsv`, resumo e manifesto

`linha_de_portao` por ligante, com os quatro critérios de ligante. `Contador`
alimenta o resumo em stderr, que diz **quantos cada critério reprovaria se
promovido** — o botão medido antes de apertado. Manifesto: parâmetros,
`versoes` (ADMET-AI, RDKit, serviços, MD5 das entradas), política, contagens
por veredito.

Chaves da repetição (`verificar_manifesto` já aceita por parâmetro):
- versões: `open_targets`, `openfda`, `admet_ai`, `rdkit`, `ligantes_md5`,
  `alvos_md5`, `evidencias_md5`, `tecidos_md5`;
- parâmetros: `politica`, `percentil_toxicidade`, `limite_bula`,
  `corte_expressao`, `termos_alerta`, `limite`, `sem_rede`.

`hashes_de_entrada_ausentes` com as entradas desta etapa (`ligantes_md5`
sempre; `alvos_md5`, `evidencias_md5`, `tecidos_md5` quando o parâmetro
correspondente foi usado).

### 6. Testes (`tests/test_portao.py`, `tests/test_portao_cli.py`)

Os cinco critérios de verificação do desenho, como teste:

1. **bula** — fixtures de resposta do openFDA para o controle: alerta em
   ciprofloxacino (tarja preta), linezolida, etambutol e metronidazol;
   **nenhum** nos seis que tratam (pregabalina e duloxetina batem em seção
   forte mas são indicadas); cloranfenicol e etionamida `indeterminado`.
2. **marca, não reprova** — ligante com BRENK sai `aprovado = sim` com
   `alerta_estrutural_veredito = marca`.
3. **política** — `--reprovar alerta_estrutural` reprova o mesmo ligante e o
   resumo diz quantos; `--marcar alerta_bula` para de reprovar.
4. **degradação** — sem `admet_ai`, código 0, colunas de ADMET vazias,
   `toxicidade = indeterminado` (monkeypatch em `admet.disponivel`).
5. **repetição** — `--manifesto` reproduz byte a byte e recusa quando entrada,
   versão ou política mudam.

Mais: rede bloqueada na suíte (o `conftest.py` já bloqueia socket), escrita
atômica (falha no meio não publica nada), e `indeterminado` nunca reprova com
política nenhuma.

Fixtures de rede em `tests/dados/respostas/`, no formato que o `Cache` grava,
com a licença já declarada em `tests/dados/LICENSE` (openFDA é domínio
público — acrescentar a nota).

### 7. README

Seção da etapa 2 depois da etapa 1: o que ela decide, a tabela de critérios
com a política padrão e o porquê medido de cada uma, as quatro saídas e suas
colunas, e o comando completo de exemplo. Marcar `checar_alertas_adversos.py`
e o caminho Deep-PK de `farmaco_completo.py` como substituídos.

## Ordem e verificação

Tarefas 1→5 em sequência (cada uma deixa a suíte verde), 6 junto de cada uma
(TDD: o teste do critério antes do critério), 7 no fim. Verificação final:
suíte completa verde, `input.exemplo.tsv` rodando sem rede, e o controle de 12
fármacos reproduzindo os números do desenho.

## Fora de escopo

Ranquear (fase 3), ligar `ranking_multibase.py`/`main.py` (fase 3), treinar ou
calibrar ADMET, e remover `farmaco_completo.py`/`checar_alertas_adversos.py`
(a remoção é da fase 3, quando nada mais os chama).
