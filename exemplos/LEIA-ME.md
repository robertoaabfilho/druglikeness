# Fixtures de validação

Material usado para medir o `mods/modelo.py`. Os números citados na
seção 8 do README saem daqui, e podem ser reproduzidos.

| arquivo | o que é |
|---|---|
| `controle.dor_neuropatica.tsv` | 12 fármacos: 6 que **tratam** dor neuropática, 6 que a **causam** (neuropatia é efeito adverso documentado). SMILES e InChIKey vindos do PubChem. |
| `alvos.dor_neuropatica.txt` | 9 alvos **humanos** de dor neuropática, por ID ChEMBL. |
| `alvos.patogeno.txt` | alvos **não humanos** (*M. tuberculosis*, *E. coli*, *S. aureus*), exercitando a sintaxe `termo \| organismo`. |

O conjunto-controle é um experimento natural: contém etambutol e etionamida
(antituberculose) ao lado de ciprofloxacino e linezolida (antibacterianos de
outras classes) e de analgésicos que não têm ação antimicrobiana. Serve para
medir os dois caminhos — humano e não humano — com o mesmo conjunto.

## Reproduzindo

```bash
# alvos humanos
python mods/modelo.py \
    --entrada exemplos/controle.dor_neuropatica.tsv \
    --checkpoint modelo//models/ecfp4_fold_4.pt \
    --rotulos modelo//target_labels.tsv \
    --alvos exemplos/alvos.dor_neuropatica.txt --organismo humano \
    --saida results/pred_humano.tsv

# alvos de patógeno
python mods/modelo.py \
    --entrada exemplos/controle.dor_neuropatica.tsv \
    --checkpoint modelo//models/ecfp4_fold_4.pt \
    --rotulos modelo//target_labels.tsv \
    --alvos exemplos/alvos.patogeno.txt \
    --saida results/pred_patogeno.tsv
```

## O que foi medido (ECFP4, `ecfp4_fold_4.pt`)

Escore do composto = maior probabilidade entre os alvos pedidos.
AUC = P(composto que trata pontua acima de composto que causa).

| escore | AUC |
|---|---|
| só os 9 alvos humanos | 0,750 |
| `p_humano × (1 − p_patógeno)` | 0,806 |
| só `CHEMBL360` (*M. tuberculosis*), separando antibacteriano de analgésico | 1,000 |

Acertos: duloxetina 0,897 e amitriptilina 0,647 no transportador de
noradrenalina (as duas são inibidores de recaptação); ciprofloxacino em 1º na
DNA girase de *E. coli*, que é o alvo real dele.

Falhas, e por quê: gabapentina (0,0001) e pregabalina (0,02) são análogos de
aminoácido pequenos e polares, e a sub-unidade α2δ-2 tem só 30 moléculas
ativas no treino — quase toda a perda de AUC do conjunto humano vem daí.
Etionamida pontua 0,80 no fenótipo de *M. tuberculosis* e 0,0017 na InhA,
porque é pró-fármaco: a molécula-mãe não se liga à InhA.
