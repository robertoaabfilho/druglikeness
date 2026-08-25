# farmaco_completo.py

Pipeline em Python que une, para uma lista de moléculas em SMILES:

1. **Descritores de drogabilidade locais** (RDKit) — peso molecular, doadores/aceptores de ligação de hidrogênio, TPSA, ligações rotacionáveis, refratividade molar, átomos pesados, Regra dos 5 (parcial) e Regra de Veber.
2. **Predições ADMET completas** (Absorção, Distribuição, Metabolismo, Excreção, Toxicidade) via API do [Deep-PK](https://biosig.lab.uq.edu.au/deeppk/).
3. **Padronização automática de SMILES com sais/pares iônicos** (`.`), já que tanto alguns cálculos quanto a API do Deep-PK rejeitam SMILES com múltiplos fragmentos desconectados.

O resultado final combinado é salvo em `farmaco.csv`.

## Requisitos

```bash
pip install rdkit requests
```

Python 3.10+ (usa sintaxe `list[tuple[str, str]]` e `X | None`).

## Uso básico

1. Crie um `input.txt` com um SMILES por linha (ver formatos abaixo).
2. Rode:

```bash
python farmaco_completo.py --input input.txt --output farmaco.csv --pred-type admet --batch-size 5
```

3. Acompanhe o progresso no terminal — `farmaco.csv` é atualizado incrementalmente a cada lote processado, não só no final.

Se o processo for interrompido (`Ctrl+C`, queda de conexão, etc.), rode o mesmo comando novamente: a consulta ao Deep-PK retoma do checkpoint salvo, sem reconsultar o que já foi obtido.

## Formato do `input.txt`

Um composto por linha. Formatos aceitos (pode misturar entre linhas):

```
CC(=O)Oc1ccccc1C(=O)O
Aspirina,CC(=O)Oc1ccccc1C(=O)O
Aspirina;CC(=O)Oc1ccccc1C(=O)O
Aspirina<TAB>CC(=O)Oc1ccccc1C(=O)O
```

- Linhas vazias ou iniciadas com `#` são ignoradas.
- Se nenhum nome for dado, a molécula recebe um nome automático (`Molecula_1`, `Molecula_2`, ...).
- SMILES com sais/pares iônicos (contendo `.`) são aceitos normalmente — ver seção seguinte.

## Padronização automática de SMILES com `.`

SMILES como `CCCC(CCC)C(=O)[O-].[Na+]` (valproato de sódio) representam mais de um fragmento desconectado. Isso é correto quimicamente, mas quebra tanto certos cálculos quanto a API do Deep-PK, que espera uma única molécula conectada.

Sempre que uma linha do `input.txt` contém `.`, o script:

1. Identifica todos os fragmentos.
2. Fica com o **maior fragmento** (heurística padrão de pipelines de padronização como ChEMBL/MolVS) — via `rdMolStandardize.LargestFragmentChooser`.
3. **Neutraliza** cargas remanescentes do fragmento escolhido (ex.: carboxilato → ácido carboxílico) — via `rdMolStandardize.Uncharger`.

Esse SMILES padronizado é o que é efetivamente usado nos descritores do RDKit **e** no que é enviado ao Deep-PK. O SMILES original digitado no `input.txt` é sempre preservado, para rastreabilidade.

### Ressalva: contra-íons metálicos

"Maior fragmento" é uma convenção mecânica, não uma avaliação farmacológica. Funciona bem quando o contra-íon é inerte — caso de Na⁺/K⁺ em sais de ácidos/bases orgânicos (ex.: valproato de sódio → ácido valproico, sem aviso). Mas quando o fragmento descartado contém outro metal (Zn, Li, Fe, Ca, Mg, Cu, Bi, Pt, ...), ele pode ser a própria espécie farmacologicamente ativa — caso do acetato de zinco, onde o Zn²⁺ é o princípio ativo e o acetato é só o contra-íon de solubilidade. Nesse cenário, prever PK sobre o ácido acético remanescente não faz sentido clínico sozinho.

Por isso, essas linhas **não são descartadas nem processadas em silêncio**: o script sinaliza um aviso na coluna `Aviso_Padronizacao` e imprime um alerta no terminal, para checagem manual.

## Saída: `farmaco.csv`

CSV com separador `;`. Colunas, na ordem:

| Coluna | Descrição |
|---|---|
| `Nome` | Nome informado no `input.txt` (ou gerado automaticamente) |
| `SMILES` | SMILES original, exatamente como veio do `input.txt` |
| `SMILES_Usado` | SMILES efetivamente analisado — igual ao original se não houver `.`; padronizado (maior fragmento + neutralizado) caso contrário |
| `Multi_Fragmento` | `True` se o SMILES original tinha `.` (múltiplos fragmentos) |
| `Contra_Ions_Removidos` | Fragmentos descartados na padronização, separados por `;` |
| `Aviso_Padronizacao` | Alerta de checagem manual quando um contra-íon metálico "não habitual" foi removido (vazio na maioria dos casos) |
| `Valido` | `True`/`False` — se o SMILES (padronizado) pôde ser interpretado pelo RDKit |
| `PesoMolecular` | Peso molecular (g/mol) |
| `DoadoresHB` | Nº de doadores de ligação de hidrogênio |
| `AceptoresHB` | Nº de aceptores de ligação de hidrogênio |
| `TPSA` | Área de superfície polar topológica (Ų) |
| `LigacoesRotacionaveis` | Nº de ligações rotacionáveis |
| `RefratividadeMolar` | Refratividade molar |
| `AtomosPesados` | Nº de átomos pesados (não-H) |
| `Violacoes_Ro5_parcial` | Nº de violações da Regra dos 5 (PM, HBD, HBA — sem LogP) |
| `Passa_Ro5_parcial` | `True` se ≤ 1 violação |
| `Passa_Veber` | `True` se TPSA ≤ 140 e ligações rotacionáveis ≤ 10 |
| `Geral_*`, `Absorcao_*`, `Distribuicao_*`, `Metabolismo_*`, `Excrecao_*`, `Toxicidade_*` | Parâmetros ADMET retornados pelo Deep-PK, uma coluna por parâmetro (descobertos dinamicamente na primeira resposta da API), no formato `predição (interpretação)` quando há interpretação disponível |
| `Outras_Propriedades` | Parâmetros retornados pelo Deep-PK que não caíram em nenhuma categoria reconhecida |

Moléculas com `Valido = False` (SMILES inválido, mesmo após padronização) não são enviadas ao Deep-PK e ficam com as colunas ADMET vazias.

## Opções de linha de comando

| Opção | Padrão | Descrição |
|---|---|---|
| `--input` | `input.txt` | Arquivo de entrada com os SMILES |
| `--output` | `farmaco.csv` | CSV final combinado |
| `--raw-deeppk` | `<output>.deeppk_raw.csv` | CSV intermediário com as respostas brutas do Deep-PK (permite retomada) |
| `--checkpoint` | `deeppk_progress.json` | Arquivo de checkpoint da consulta ao Deep-PK |
| `--pred-type` | `admet` | `absorption`, `distribution`, `metabolism`, `excretion` ou `admet` (todas as categorias) |
| `--batch-size` | `5` | Quantos SMILES por job enviado ao Deep-PK |
| `--email` | — | E-mail opcional para notificação da API ao final do job |

## Como funciona por baixo dos panos

1. **Leitura** (`ler_smiles_do_arquivo`): parseia o `input.txt` em pares `(nome, smiles)`.
2. **Descritores locais** (`calcular_descritores` + `processar_lista_druglikeness`): para cada molécula, padroniza o SMILES se necessário (`padronizar_smiles`), calcula os descritores RDKit sobre o SMILES padronizado e avalia Ro5 parcial / Veber. `farmaco.csv` já é salvo neste ponto, mesmo antes de consultar o Deep-PK.
3. **Consulta ao Deep-PK** (`consultar_deeppk`): envia os SMILES válidos (já padronizados) em lotes pequenos. Se um lote inteiro falhar, ele é dividido recursivamente ao meio (bissecção) até isolar exatamente qual SMILES está quebrando o job — o restante do lote continua sendo processado normalmente. Cada sucesso (lote ou sublote) é gravado incrementalmente no CSV intermediário e `farmaco.csv` é regravado com o progresso combinado.
4. **Combinação final** (`combinar_e_salvar`): junta descritores locais + resultados do Deep-PK pela posição original de cada molécula no `input.txt` e escreve `farmaco.csv`.

## Limitações conhecidas

- A padronização por "maior fragmento" é uma heurística mecânica: para sais onde o contra-íon é a espécie ativa (metais como Zn, Li, Fe...), a predição de PK resultante não é apropriada — sempre revise as linhas com `Aviso_Padronizacao` preenchido.
- A API do Deep-PK é externa e pode ficar indisponível ou lenta; o script espera até 30 minutos por job (`MAX_WAIT_SECONDS`) antes de desistir e isolar o lote por bissecção.
- Moléculas totalmente inorgânicas (sem carbono após a padronização, ex.: sais minerais puros) não são um caso de uso adequado para preditores de ADMET baseados em QSAR/deep learning treinados em química orgânica.
