# `referencias_alvo.tsv` — os ativos já conhecidos de cada alvo afirmativo

18.286 pares alvo-molécula: para cada um dos 9 alvos de
`input/alvos.txt`, as moléculas com atividade medida naquele alvo. Servem de
**referência de domínio**: é a distância estrutural a elas que decide se a
previsão do modelo tem habilidade medida ou é azar (`mods/ranque/dominio.py`).

| alvo | nome | referências |
|---|---|---|
| CHEMBL233 | receptor opioide mu | 5.790 |
| CHEMBL228 | transportador de serotonina (SERT) | 4.434 |
| CHEMBL4794 | TRPV1 | 3.085 |
| CHEMBL222 | transportador de noradrenalina (NET) | 2.856 |
| CHEMBL1904 | NMDA 2B | 1.549 |
| CHEMBL1919 | canal de cálcio α2δ-1 | 373 |
| CHEMBL2015 | NMDA 1 | 113 |
| CHEMBL1845 | canal de sódio tipo 1 α | 67 |
| CHEMBL3896 | canal de cálcio α2δ-2 | 19 |

## Procedência

API do ChEMBL (`/activity`), consultada **por ID de alvo** — não por nome, para
não casar com ortólogo de rato ou camundongo. Filtros:

* `pchembl_value >= 6.0` (≤ 1 µM)
* `assay_type` em `B` (ligação) ou `F` (funcional). Varredura HTS de painel
  (`A`, `T`, `P`, `U`) fica fora: o próprio projeto mediu que ela é ruído —
  recall@10 de 0,111 em HTS contra 0,272 em ensaio dirigido.
* uma linha por `molecule_chembl_id`, o primeiro SMILES canônico visto

Cliente genérico, sem identificar o usuário ao serviço externo.

## Por que os quatro alvos de baixo precisam de tratamento diferente

`sim_top5` é a média das 5 **maiores** similaridades. Com pouca referência, as
5 maiores são simplesmente as únicas que existem, e o valor cai por construção
— não porque o composto seja distante do conhecido, mas porque pouco é
conhecido. Medido contra os 107.123 ligantes de `input/ligantes.tsv`:

| referências | fração da biblioteca com `sim_top5 < 0,20` |
|---|---|
| 5.790 (CHEMBL233) | 28,5% |
| 1.549 (CHEMBL1904) | 75,7% |
| 373 (CHEMBL1919) | 96,4% |
| 19 (CHEMBL3896) | 99,9% |

O salto está entre 1.549 e 373. Por isso `dominio.MIN_REFERENCIAS = 500`:
abaixo disso a zona é `indeterminado`, e não `fora`. Chamar de `fora` rebaixaria
o alvo por falta de bibliografia em vez de por falta de sinal.

## O que este arquivo NÃO é

Não é conjunto de treino nem de validação. É o referencial de distância. As
moléculas daqui podem perfeitamente estar na biblioteca triada — e quando
estão, caem na zona `analogo`, que é onde o modelo acerta mais (precisão@5 de
52,6% contra 4,3% na banda de novidade) e onde há menos a descobrir.
