# `validacao_neuro_prospector.tsv` — conjunto de validação independente

244 compostos rotulados `trata` / `causa` / `ambiguo` no eixo dor-nervo, para
validar a fase 3 (seção 13 do README) **fora** dos 12 fármacos que derivaram as
decisões de desenho dela.

| grupo | compostos | com `chembl_id` |
|---|---|---|
| `trata` | 156 | 147 |
| `causa` | 66 | 53 |
| `ambiguo` | 22 | 20 |

## Por que ele é independente

As decisões de desenho da fase 3 — qual canal negativo usar (0,806 contra 0,694
e 0,583), o corte de 50 moléculas de treino, a exclusão de ortólogo, a exclusão
de organismo vazio — saíram todas de medições sobre
`exemplos/controle.dor_neuropatica.tsv`, 12 fármacos. Em termos metodológicos,
o 0,806 é um número de **treino**: a regra foi escolhida por maximizar a AUC
naquele conjunto, entre pelo menos cinco candidatas.

Os rótulos deste arquivo vêm de outra fonte, sobre outros compostos, e não
influenciaram nenhuma dessas escolhas.

**Uma exceção, marcada na coluna `no_controle_de_desenho`:** a **capsaicina** (`Capsaicin`,
como o nome sai no TSV) está nos dois conjuntos. Ela caiu como `ambiguo` aqui (tem termos dos dois
lados, o que é verdade farmacológica — trata dor neuropática e é irritante de
fibra C), mas a coluna existe para que a exclusão seja **verificável** em vez de
depender de inspeção visual. Qualquer validação deve filtrar
`no_controle_de_desenho == ""`.

## Procedência

**Origem:** `~/prospector/pipeline_prospeccao_v4/resultados/atividades_compostos_parcial.csv`
— projeto **separado**, lido em modo somente-leitura. MD5 em 2026-11-26:
`eac618972d4900e0bd684597b49f654c` (126.240 compostos, 137 MB).

**Rótulos:** coluna `atividades_duke`, da Dr. Duke's Phytochemical and
Ethnobotanical Database (USDA/ARS) — vocabulário **controlado**, com os dois
sentidos do mesmo eixo. Dos 126.240 compostos, **2.593** têm alguma atividade
Duke; destes, 264 têm algum termo das listas abaixo, e 20 foram descartados por
SMILES que o RDKit não lê.

**Química:** `nome`, `smiles`, `inchikey` e `inchikey_parental` vêm de
`mods/quimica.converter_linha`, o mesmo tratamento que a importação
normal aplica — então as chaves batem com as do pipeline. A deduplicação é por
`inchikey_parental`.

## A regra de seleção, completa

Um composto entra se `atividades_duke` tem ao menos um termo de uma das listas.
`trata` e `causa` exclusivos; termos dos dois lados → `ambiguo`.

**`trata`** (ação terapêutica no eixo dor/nervo):

```
Analgesic, Analgesic-Synergist, Antineuralgic, Antineuropathic, Antineuritic,
Antipolyneuritic, Anticonvulsant, Antiepileptic, Neuroprotective, Antimyalgic,
Antifibromyalgic, Antirheumatalgic, Antiodontalgic, Antimastalgic,
Antineurasthenic, Antineurotic
```

**`causa`** (dano neurológico):

```
Neurotoxic, Convulsant, Epileptigenic, Neuroparalytic, CNS-Paralytic,
Paralytic, Myoparalytic, Respiraparalytic, Cardioparalytic, Neuroinhibitor
```

**Deliberadamente fora das duas listas:** `Sedative`, `Narcotic`, `Anesthetic`,
`Myorelaxant`, `Neurotropic`, `Neurotransmitter`, `Neurogenic`. Eles descrevem
efeito no SNC **sem afirmar lesão**, e são exatamente o que as drogas que
*tratam* dor neuropática também fazem — a amitriptilina é sedativa. Incluí-los
em `causa` repetiria o erro do PAINS na etapa 2: apagar o positivo.

## Limitações — leia antes de usar

1. **O rótulo não é do mesmo nível de evidência do controle de 12.** Lá,
   `causa` vem de tarja preta e seção de advertência da FDA (openFDA). Aqui vem
   de alegações da literatura etnobotânica agregadas pelo Duke, **sem dose, sem
   via de administração e sem controle**. `(-)-Pinene` aparece como `causa`
   porque algum trabalho reportou neurotoxicidade — provavelmente em
   concentração que ninguém atinge por uso normal.
2. **Enantiômeros com rótulos opostos.** `(-)-Borneol` é `trata` (Analgesic) e
   `(+)-Borneol` é `causa` (Convulsant). Pode ser estereoquímica real ou ruído
   do Duke; a deduplicação por forma parental preserva a estereoquímica, então
   os dois ficam como compostos distintos. Trate pares assim como suspeitos.
3. **Desbalanceado:** 156 `trata` contra 66 `causa`. Use métrica que não
   dependa de prevalência (AUC serve; acurácia não).
4. **Cobertura rala:** só 2% dos compostos da prospecção têm atividade Duke, e
   quem tem tende a ser o produto natural mais estudado. Há viés de seleção para
   compostos famosos.
5. **Não há alvo aqui.** Este arquivo rotula o **desfecho** (trata/causa), não o
   alvo. Para ranquear, a fase 3 precisa de `evidencias.tsv` sobre esses
   compostos — ou seja, rodar a etapa 1 com eles e uma lista de alvos.
6. **Não é um conjunto de teste cego.** Eu escolhi as duas listas de termos
   olhando o vocabulário do Duke. Não ajustei nada da fase 3 depois de ver o
   resultado, e o compromisso está registrado neste arquivo — mas a escolha dos
   termos é um grau de liberdade que ficou comigo.

## Como regenerar

A regra acima é a especificação completa. O script que a aplicou vive fora do
repositório (foi um utilitário de sessão); se este arquivo precisar ser
regenerado com outra lista de termos ou outra carga da prospecção, vale
promovê-lo a `mods/`, para que o conjunto tenha a mesma garantia de repetição
que as tabelas das etapas 1 a 3.

## Licença

Ver `dados/LICENSE`. Resumo: rótulos do Duke são domínio público (USDA); os
SMILES vêm de NPASS, COCONUT, LOTUS e PubChem, com licenças próprias — **o
NPASS é acadêmico/não comercial**. `chembl_id` é identificador, não dado de
bioatividade, então o share-alike do ChEMBL não alcança este arquivo.
