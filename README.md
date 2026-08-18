# farmaco_completo.py

Script que une **descritores de drogabilidade locais (RDKit)** com **predições ADMET completas da API do Deep-PK** em um único fluxo, a partir de uma lista de SMILES, salvando tudo combinado em `farmaco.csv`.

---

## Índice

- [O que o script faz](#o-que-o-script-faz)
- [Requisitos](#requisitos)
- [Instalação](#instalação)
- [Formato do arquivo de entrada (`input.txt`)](#formato-do-arquivo-de-entrada-inputtxt)
- [Como usar](#como-usar)
- [Argumentos de linha de comando](#argumentos-de-linha-de-comando)
- [Arquivos gerados](#arquivos-gerados)
- [Colunas do `farmaco.csv`](#colunas-do-farmacocsv)
- [Como funciona por dentro](#como-funciona-por-dentro)
  - [Etapa 1 — Descritores de drogabilidade (RDKit)](#etapa-1--descritores-de-drogabilidade-rdkit)
  - [Etapa 2 — Consulta ao Deep-PK](#etapa-2--consulta-ao-deep-pk)
  - [Etapa 3 — Combinação e salvamento incremental](#etapa-3--combinação-e-salvamento-incremental)
- [Resiliência: checkpoint, retomada e isolamento de falhas](#resiliência-checkpoint-retomada-e-isolamento-de-falhas)
- [Barras de progresso e mensagens de log](#barras-de-progresso-e-mensagens-de-log)
- [Perguntas frequentes / solução de problemas](#perguntas-frequentes--solução-de-problemas)
- [Limitações conhecidas](#limitações-conhecidas)

---

## O que o script faz

Dado um arquivo `input.txt` com uma lista de moléculas (em SMILES), o script:

1. **Lê** o arquivo de entrada, aceitando SMILES puro ou acompanhado de um nome.
2. **Calcula, localmente e instantaneamente (via RDKit)**, os descritores clássicos de drogabilidade: peso molecular, doadores/aceptores de ligação de hidrogênio, TPSA, ligações rotacionáveis, refratividade molar, átomos pesados, além das avaliações da **Regra dos 5 de Lipinski** (parcial, sem LogP) e da **Regra de Veber**.
3. **Consulta a API pública do [Deep-PK](https://biosig.lab.uq.edu.au/deeppk/)** para obter as predições ADMET completas (propriedades Gerais, Absorção, Distribuição, Metabolismo, Excreção e Toxicidade) das mesmas moléculas — em lotes pequenos, com checkpoint/retomada e tratamento automático de falhas.
4. **Junta os dois conjuntos de resultados** pela posição original de cada molécula no arquivo de entrada e salva tudo combinado em `farmaco.csv`, **atualizado incrementalmente a cada lote processado**.

Em resumo: você fornece uma lista de SMILES e recebe uma única planilha com drogabilidade local + ADMET completo, sem precisar rodar dois scripts separados nem submeter molécula por molécula manualmente no site do Deep-PK.

---

## Requisitos

- **Python 3.10+** (o script usa sintaxe de tipos moderna, como `list[str]` e `str | None`)
- Pacotes Python:
  - [`rdkit`](https://www.rdkit.org/) — cálculo dos descritores de drogabilidade
  - [`requests`](https://requests.readthedocs.io/) — chamadas HTTP à API do Deep-PK
- **Conexão com a internet** liberada para `https://biosig.lab.uq.edu.au` (a etapa de RDKit funciona offline; só a consulta ADMET precisa de rede)

## Instalação

```bash
pip install rdkit requests
```

> Em ambientes gerenciados (Linux com Python "externally managed"), pode ser necessário:
> ```bash
> pip install rdkit requests --break-system-packages
> ```

---

## Formato do arquivo de entrada (`input.txt`)

Um SMILES por linha. Você pode misturar os formatos abaixo livremente entre linhas diferentes:

```text
CC(=O)Oc1ccccc1C(=O)O
Aspirina,CC(=O)Oc1ccccc1C(=O)O
Aspirina;CC(=O)Oc1ccccc1C(=O)O
Aspirina	CC(=O)Oc1ccccc1C(=O)O
```

Regras:

- **Separadores aceitos entre nome e SMILES:** vírgula (`,`), ponto e vírgula (`;`) ou tabulação (`\t`) — testados nessa ordem.
- **Se nenhum nome for informado**, a molécula recebe um nome automático: `Molecula_1`, `Molecula_2`, ...
- **Linhas vazias são ignoradas.**
- **Linhas iniciadas com `#` são tratadas como comentário** e ignoradas.
- A **ordem das linhas é preservada** e usada internamente como "índice original" para juntar corretamente os resultados do Deep-PK com os descritores locais — não reordene o arquivo entre execuções que usam o mesmo checkpoint.

---

## Como usar

Uso mais simples, com todos os valores padrão (lê `input.txt`, grava `farmaco.csv` no diretório atual):

```bash
python farmaco_completo.py
```

Especificando entrada e saída:

```bash
python farmaco_completo.py --input minhas_moleculas.txt --output resultado.csv
```

Rodando só a predição de Absorção (mais rápido que `admet` completo) e lotes maiores:

```bash
python farmaco_completo.py --pred-type absorption --batch-size 10
```

---

## Argumentos de linha de comando

| Argumento | Padrão | Descrição |
|---|---|---|
| `--input` | `input.txt` | Arquivo de entrada com os SMILES (um por linha). |
| `--output` | `farmaco.csv` | Arquivo CSV final, com os descritores locais + ADMET combinados. |
| `--raw-deeppk` | `<output_sem_extensao>.deeppk_raw.csv` | CSV intermediário onde as respostas brutas do Deep-PK são acumuladas (ver [Arquivos gerados](#arquivos-gerados)). |
| `--checkpoint` | `deeppk_progress.json` | Arquivo de checkpoint da consulta ao Deep-PK, usado para retomar uma execução interrompida. |
| `--pred-type` | `admet` | Tipo de predição no Deep-PK: `admet` (tudo), `absorption`, `distribution`, `metabolism` ou `excretion`. |
| `--batch-size` | `5` | Quantos SMILES enviar por job à API do Deep-PK (também é o tamanho inicial de lote antes de qualquer bissecção por falha — ver [Isolamento de falhas](#resiliência-checkpoint-retomada-e-isolamento-de-falhas)). |
| `--email` | *(nenhum)* | E-mail opcional passado à API do Deep-PK para notificação ao final do job. |

---

## Arquivos gerados

| Arquivo | O que é | Precisa manter depois? |
|---|---|---|
| **`farmaco.csv`** (ou o nome passado em `--output`) | **Resultado final** — descritores de drogabilidade + ADMET combinados, uma linha por molécula. | ✅ Sim, é o que você quer. |
| `*.deeppk_raw.csv` (ex: `farmaco.deeppk_raw.csv`) | **Cache intermediário** com as respostas brutas do Deep-PK por molécula, indexadas pela posição original no `input.txt`. Usado para (a) permitir retomada sem repetir chamadas à API já bem-sucedidas e (b) remontar o `farmaco.csv` a cada lote. | ⚠️ Só enquanto for continuar/repetir a execução no mesmo `input.txt`. Pode apagar com segurança depois que o `farmaco.csv` estiver completo e satisfatório. |
| `deeppk_progress.json` (ou o nome passado em `--checkpoint`) | Checkpoint: guarda o índice do próximo SMILES **válido** a consultar no Deep-PK. É apagado automaticamente quando a consulta termina com sucesso. | ⚠️ Mesma lógica do arquivo acima — só é útil enquanto a consulta ainda não terminou. |

> Se o script for interrompido (erro de rede, `Ctrl+C`, queda de energia), **não apague** o CSV intermediário nem o checkpoint antes de rodar de novo — é exatamente isso que permite retomar sem perder trabalho já feito.
>
> Se quiser recomeçar do zero (por exemplo, depois de editar o `input.txt`), apague os três arquivos antes de rodar novamente:
> ```bash
> rm -f farmaco.csv farmaco.deeppk_raw.csv deeppk_progress.json
> ```

---

## Colunas do `farmaco.csv`

O separador do CSV final é **ponto e vírgula (`;`)** — não vírgula — porque vários campos ADMET trazem texto com vírgulas embutidas (interpretações, faixas de valores, etc.).

### Descritores de drogabilidade (RDKit, sempre presentes)

| Coluna | Descrição |
|---|---|
| `Nome` | Nome da molécula (informado no `input.txt` ou gerado automaticamente). |
| `SMILES` | SMILES original da molécula. |
| `Valido` | `True`/`False` — se o RDKit conseguiu interpretar o SMILES. |
| `PesoMolecular` | Peso molecular (g/mol). |
| `DoadoresHB` | Nº de doadores de ligação de hidrogênio. |
| `AceptoresHB` | Nº de aceptores de ligação de hidrogênio. |
| `TPSA` | Área de superfície polar topológica (Å²). |
| `LigacoesRotacionaveis` | Nº de ligações rotacionáveis. |
| `RefratividadeMolar` | Refratividade molar. |
| `AtomosPesados` | Nº de átomos pesados (não-hidrogênio). |
| `Violacoes_Ro5_parcial` | Nº de violações da Regra dos 5 de Lipinski (**sem** o critério de LogP, que vem do Deep-PK). |
| `Passa_Ro5_parcial` | `True` se `Violacoes_Ro5_parcial <= 1` (padrão Lipinski). |
| `Passa_Veber` | `True` se `TPSA <= 140` **e** `LigacoesRotacionaveis <= 10`. |

Quando `Valido = False` (SMILES que o RDKit não conseguiu interpretar), todas as colunas numéricas acima ficam vazias, e a molécula **não é enviada ao Deep-PK**.

### Colunas ADMET (Deep-PK, dinâmicas)

O restante das colunas é descoberto **dinamicamente** a partir da própria resposta da API na primeira consulta bem-sucedida, então o conjunto exato de colunas pode variar conforme o `--pred-type` escolhido. Seguem o padrão `<Categoria>_<Parametro>`, agrupadas nesta ordem:

- **`Geral_*`** — propriedades físico-químicas gerais (ex: `Geral_Log_P`, `Geral_Log_S`, `Geral_pKa_Acid`, `Geral_Melting_Point`...)
- **`Absorcao_*`** — absorção (ex: `Absorcao_Caco_2_logPaap`, `Absorcao_Human_Intestinal_Absorption`, `Absorcao_P_Glycoprotein_Substrate`...)
- **`Distribuicao_*`** — distribuição (ex: `Distribuicao_Blood_Brain_Barrier`, `Distribuicao_Plasma_Protein_Binding`...)
- **`Metabolismo_*`** — metabolismo, principalmente interações com enzimas CYP (ex: `Metabolismo_CYP_3A4_Inhibitor`, `Metabolismo_CYP_2D6_Substrate`...)
- **`Excrecao_*`** — excreção (ex: `Excrecao_Clearance`, `Excrecao_Half_Life_of_Drug`)
- **`Toxicidade_*`** — toxicidade (ex: `Toxicidade_hERG_Blockers`, `Toxicidade_AMES_Mutagenesis`, `Toxicidade_Carcinogenesis`...)
- **`Outras_Propriedades`** — rede de segurança: qualquer parâmetro retornado pela API que não se encaixe nas categorias reconhecidas acima cai aqui, no formato `Categoria/Propriedade: valor`.

Cada célula ADMET traz o formato **`predição (interpretação)`** quando a API fornece uma interpretação textual (ex: `Negative (Non-blocker)`), ou só a predição quando não há interpretação disponível.

Moléculas que falharam permanentemente na consulta ao Deep-PK (ver seção seguinte) ficam com **todas as colunas ADMET vazias**, mas mantêm normalmente os descritores locais de drogabilidade.

---

## Como funciona por dentro

### Etapa 1 — Descritores de drogabilidade (RDKit)

Executa **localmente, sem rede**, e é praticamente instantânea mesmo para milhares de moléculas. Para cada linha do `input.txt`:

1. Tenta interpretar o SMILES com `Chem.MolFromSmiles`.
2. Se inválido, marca `Valido = False` e preenche o restante dos descritores com vazio.
3. Se válido, calcula os descritores via `rdkit.Chem.Descriptors` e avalia a Regra dos 5 (parcial) e a Regra de Veber.

O `farmaco.csv` já é gravado uma primeira vez logo após essa etapa (só com os descritores locais), antes mesmo de começar a falar com o Deep-PK.

### Etapa 2 — Consulta ao Deep-PK

Só as moléculas com `Valido = True` são enviadas à API — um único SMILES inválido dentro de um lote costuma derrubar o job inteiro no lado do Deep-PK.

O fluxo, por lote (tamanho definido por `--batch-size`):

1. Envia o lote de SMILES para `POST /api/predict` (upload de arquivo temporário) e recebe um `job_id`.
2. Consulta periodicamente (`GET /api/predict?job_id=...`, a cada 15s, por até 30 minutos) até o job terminar, mostrando uma barra de espera em tempo real.
3. Ao terminar, valida a resposta: se a API sinalizar erro de qualquer forma (status de erro, chave `error`/`message`, ou payload sem nenhuma molécula processada), a consulta desse lote é tratada como falha.
4. Em caso de falha do lote inteiro, o script **divide o lote ao meio e tenta cada metade recursivamente** (bissecção) — isso isola automaticamente qual(is) molécula(s) específica(s) está(ão) causando o problema, sem descartar o restante do lote, que costuma ser perfeitamente válido.
5. Se uma molécula falhar mesmo sozinha, o script contorna primeiro um detalhe conhecido da API — **o endpoint de lote do Deep-PK rejeita arquivos com um único SMILES** (erro `"No valid molecules were provided..."`) — enviando-a duplicada (2×) só para satisfazer esse mínimo, e usando apenas o primeiro resultado. Se **mesmo assim** falhar, a molécula é definitivamente marcada como falha e fica sem dados ADMET.

### Etapa 3 — Combinação e salvamento incremental

A cada lote **ou sublote** resolvido com sucesso (mesmo os menores, gerados pela bissecção), o script:

1. Anexa o resultado bruto ao CSV intermediário (`*.deeppk_raw.csv`).
2. Relê esse CSV intermediário por completo.
3. Junta com os descritores de drogabilidade (pela posição original da molécula no `input.txt`).
4. Regrava o `farmaco.csv` do zero com o progresso atualizado.

Ou seja: **você pode abrir o `farmaco.csv` a qualquer momento durante a execução** e ver o progresso mais recente, em vez de esperar o processo inteiro terminar.

---

## Resiliência: checkpoint, retomada e isolamento de falhas

- **Checkpoint automático:** a cada lote de nível superior concluído (com sucesso total, parcial, ou com falhas já registradas), o índice de progresso é salvo em `deeppk_progress.json`. Se o script for interrompido, rodar o mesmo comando de novo retoma a partir daí, sem repetir consultas já concluídas.
- **Isolamento de falhas por bissecção:** um lote que falha não derruba as outras moléculas do mesmo lote — o script encontra automaticamente o(s) culpado(s) e segue em frente com o resto.
- **Contorno do mínimo de 2 SMILES por job:** evita falsos positivos de "molécula inválida" quando na verdade é só uma limitação de formato do endpoint de lote da API.
- **Nenhum CSV "vazio silencioso":** versões anteriores deste script podiam gerar um `farmaco.csv` aparentemente completo, mas com todas as colunas ADMET em branco, sem avisar sobre o motivo. Isso foi corrigido — qualquer falha real da API agora é detectada e reportada.

Ao final da execução, se alguma molécula tiver falhado permanentemente, o script imprime um resumo:

```text
[AVISO] 1 SMILES sem dados ADMET (falharam sozinhos no Deep-PK):
  - índice 4: [C@@H]([C@@H](C(=O)O)S)(C(=O)O)S.[C@@H]... — nenhuma molécula retornada: {...}
```

Nesses casos, verifique se o SMILES é multi-componente (contém `.`, indicando fragmentos separados — sais, complexos, dímeros não covalentes) ou tente submetê-lo manualmente no site do Deep-PK para investigar melhor.

---

## Barras de progresso e mensagens de log

Durante a execução, o terminal mostra:

- **Progresso por lote**, após cada lote de nível superior:
  ```text
  [INFO] Deep-PK [##############################] 100% (5/5)
  ```
- **Barra de espera** enquanto um job específico está rodando na API (atualizada na mesma linha, sem poluir o log):
  ```text
  [INFO]   Aguardando Deep-PK [######--------------] 45s
  ```
- **Avisos de bissecção**, quando um lote precisa ser dividido por falha:
  ```text
  [AVISO]   lote de 5 falhou (nenhuma molécula retornada: {'status': 'ERROR while running!'}) — dividindo...
  ```
- **Confirmação de gravação** a cada atualização do arquivo final:
  ```text
  [INFO]     'farmaco.csv' atualizado (2 molécula(s) deste sublote).
  ```

Prefixos usados: `[INFO]` para progresso normal, `[AVISO]` para situações que não interrompem a execução mas merecem atenção, `[ERRO]` para falhas que encerram o script (impressas em `stderr`).

---

## Perguntas frequentes / solução de problemas

**O `farmaco.csv` foi gerado, mas todas as colunas ADMET estão vazias. O que houve?**
Verifique o log da execução — deve haver uma linha `[AVISO] Nenhuma coluna ADMET foi obtida do Deep-PK` ou um `[ERRO]` explicando o motivo (a API pode estar fora do ar, o `--pred-type` pode não bater com as colunas esperadas, etc.). Rode de novo — como o checkpoint provavelmente ficou "preso" achando que já terminou, apague `farmaco.deeppk_raw.csv` e `deeppk_progress.json` antes de tentar novamente.

**Uma molécula específica sempre falha no Deep-PK, mas funciona no site manualmente.**
O site costuma usar um endpoint diferente do endpoint de upload em lote usado por este script. SMILES com `.` (múltiplos fragmentos — sais, dímeros, complexos) são os candidatos mais prováveis a causar esse tipo de incompatibilidade. Considere separar os fragmentos em linhas distintas do `input.txt`, se fizer sentido quimicamente.

**Como faço para reprocessar só as moléculas que falharam?**
Crie um novo `input.txt` só com essas moléculas (usando os SMILES exibidos no resumo de falhas) e rode o script normalmente com um `--output` diferente; depois, se quiser, junte manualmente as linhas nos dois CSVs.

**Posso rodar em milhares de moléculas de uma vez?**
Sim — é exatamente para isso que existem o checkpoint, a gravação incremental e o `--batch-size`. Lotes menores (`--batch-size` menor) tendem a reduzir o "raio de explosão" de uma falha, mas fazem mais chamadas à API (mais lento). O padrão (`5`) é um meio-termo razoável.

**Erros de rede / `requests.exceptions.ConnectionError`.**
Verifique sua conexão com `https://biosig.lab.uq.edu.au`. Se seu ambiente usa proxy/firewall corporativo, pode ser necessário liberar esse domínio. Rodar o script de novo retoma de onde parou.

---

## Limitações conhecidas

- O script depende inteiramente da disponibilidade e do formato de resposta da **API pública do Deep-PK**, que pode mudar sem aviso prévio.
- A bissecção automática de lotes reduz, mas não elimina, o número de chamadas extras à API quando há falhas — em arquivos com muitas moléculas problemáticas, isso pode aumentar bastante o tempo total de execução.
- O CSV intermediário (`*.deeppk_raw.csv`) não é limpo automaticamente ao final — ver [Arquivos gerados](#arquivos-gerados) para saber quando é seguro apagá-lo manualmente.
- O script assume que a ordem das linhas do `input.txt` não muda entre uma execução interrompida e sua retomada (o checkpoint e o CSV intermediário são indexados pela posição original das moléculas).
