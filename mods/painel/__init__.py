# -*- coding: utf-8 -*-
"""
mods/painel/ - a biblioteca do painel.

`painel.py`, na raiz, e so a casca: as duas abas e os widgets. Tudo o que
decide alguma coisa mora aqui, em modulos que importam e sao testados SEM
servidor de Streamlit:

    leitura.py    as tabelas publicadas e os manifestos, com cache
    execucao.py   disparar o main.py, saber se esta de pe, medir o ritmo
    entradas.py   receber upload em input/ e conferir antes de rodar
    vistas.py     recortes das tabelas para a tela (funcoes puras)
    alvos.py      revisar a lista de alvos enriquecidos

O painel NAO reimplementa a sequencia das etapas: ele le a tabela `ETAPAS` do
`main.py`. Orquestracao morando no Streamlit e o erro que o desenho nomeia --
a sequencia viraria codigo de interface e sairia do alcance de qualquer
script.
"""
