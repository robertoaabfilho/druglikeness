# -*- coding: utf-8 -*-
"""
esquema.py

As colunas das saidas da fase 3. O `EscritorTSV` de `mods/portao/esquema.py` e
o `ler_tsv` dele servem esta fase sem alteracao -- o `#` do cabecalho e do
escritor, nao do nome da coluna, por isso nenhum nome aqui comeca com ele.
"""

COLUNAS_RANKING = [
    "posicao", "ligante", "ligante_inchikey_parental", "faixa", "valor",
    "alvo", "p_afirmativo", "pchembl", "p_negativo", "alvo_negativo",
    "medida_sem_pchembl", "zona_dominio", "portao_aprovado",
    "portao_motivos", "observacao", "procedencia",
]
# `p_afirmativo` e `pchembl` saem SEPARADAS, e so uma delas e preenchida: a
# faixa 1 ordena por pchembl (medido) e a faixa 2 por probabilidade
# (previsto). Uma coluna "valor_bruto" unica poria as duas grandezas no mesmo
# lugar, que e exatamente o que a etapa 1 proibe.
#
# `zona_dominio` sai porque ela MUDA A ORDEM dentro da faixa prevista: sem a
# coluna, duas linhas vizinhas com `p_afirmativo` em ordem invertida
# pareceriam defeito. Vazia na faixa medida, que nao olha o dominio.

COLUNAS_NAO_RANQUEADOS = [
    "ligante", "ligante_inchikey_parental", "motivo", "faixa",
    "portao_aprovado", "portao_motivos", "procedencia",
]
# Sem `posicao`: quem nao ranqueia nao tem posicao, e uma coluna vazia
# convidaria a ordenar esta tabela como se fosse ranking.

AUSENTE_DO_PORTAO = "ausente"
PORTAO_NAO_CONSULTADO = "nao consultado"
MOTIVO_FAIXA_3 = "sem evidencia em alvo afirmativo"
