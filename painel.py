# -*- coding: utf-8 -*-
"""
painel.py - a interface do pipeline: Operar e Analisar.

Duas abas, e a separacao entre elas e deliberada:

* **Operar** escreve. Entradas (`input/`), parametros (`json/config.json`), e
  disparar a execucao.
* **Analisar** so le. As tabelas de `tsv/`, com o manifesto de cada uma.

Evita o acidente de mexer num parametro enquanto se interpreta um resultado e
depois nao saber qual numero gerou o que.

**Esta casca nao decide nada.** Marcar, reprovar, ordenar e enriquecer ficam
nas etapas; a sequencia das etapas fica no `main.py`, de onde esta tela a le.
Toda a logica mora em `mods/painel/`, que se testa sem servidor.

**O painel le o sistema de arquivos.** Execucao disparada no terminal aparece
aqui; execucao disparada aqui e lida por um script. O sistema de arquivos e o
unico acoplamento entre as duas portas.
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import streamlit as st  # noqa: E402

import main  # noqa: E402
from mods import caminhos, config, dependencias  # noqa: E402
from mods.painel import alvos as palvos  # noqa: E402
from mods.painel import entradas, execucao, leitura  # noqa: E402
from mods.painel import estrutura, parametros, vistas  # noqa: E402
from mods.ranque.esquema import MOTIVO_FAIXA_3  # noqa: E402
from mods.ranque.faixa import FAIXAS  # noqa: E402

st.set_page_config(page_title="Pipeline druglikeness", page_icon="🧪",
                   layout="wide")

ROTULO_DO_PAPEL = {
    "ligantes": "Ligantes (TSV com `nome` e `smiles`)",
    "alvos": "Alvos (um ID do ChEMBL por linha, ou TSV)",
    "sintomas": "Sintomas (um por linha)",
}
DESCRICAO = {e["id"]: e["descricao"] for e in main.ETAPAS}


def _duracao(segundos: float) -> str:
    if segundos < 60:
        return f"{segundos:.0f}s"
    if segundos < 3600:
        return f"{segundos / 60:.0f}min"
    return f"{segundos / 3600:.1f}h"


def _cfg():
    """A configuracao, ou `None` com o erro na tela. Nao levanta: sem isso a
    aba inteira ficaria em branco por causa de uma chave mal escrita."""
    try:
        return config.ler()
    except config.ErroDeConfig as erro:
        st.error(f"`json/config.json` nao pode ser lido: {erro}")
        return None


# =========================================================================
# Operar
# =========================================================================
def bloco_entradas():
    st.subheader("1. Entradas")
    st.caption("As tres listas sao suas; o pipeline nao deduz nenhuma delas. "
               "Ha exemplos em `exemplos/`.")
    estado = entradas.estado()
    for papel, rotulo in ROTULO_DO_PAPEL.items():
        atual = estado[papel]
        colunas = st.columns([3, 2])
        with colunas[0]:
            envio = st.file_uploader(rotulo, key=f"envio_{papel}")
        with colunas[1]:
            if atual["existe"]:
                st.caption(
                    f"`{atual['caminho'].name}` · {atual['bytes']:,} bytes · "
                    f"md5 `{atual['md5'][:8]}` · "
                    f"{time.strftime('%d/%m %H:%M', time.localtime(atual['mtime']))}")
            else:
                st.caption("ainda nao fornecido")
        if envio is not None and st.button(f"Gravar em input/", key=f"gravar_{papel}"):
            relato = entradas.receber(envio, papel)
            if relato["substituiu"]:
                st.warning(
                    f"substituiu o arquivo anterior "
                    f"(md5 `{relato['md5_anterior'][:8]}` -> "
                    f"`{relato['md5'][:8]}`)")
            else:
                st.success(f"gravado em `{relato['destino']}`")
            _mostrar_conferencia(relato["conferencia"])
        if atual["existe"] and st.button("Conferir", key=f"conferir_{papel}"):
            _mostrar_conferencia(entradas.conferir(atual["caminho"], papel))


def _mostrar_conferencia(c: dict):
    st.write(f"**{c['aceitos']:,}** aceitos · {c['erros']} erros · "
             f"{c['avisos']} avisos")
    if c["mensagens"]:
        with st.expander(f"{c['problemas']} problema(s) — primeiros "
                         f"{len(c['mensagens'])}"):
            for mensagem in c["mensagens"]:
                st.text(mensagem)


def bloco_alvos_enriquecidos(cfg):
    """Revisar a lista que a etapa 1 enriqueceu. FERRAMENTA, nao portao.

    Fica em Operar porque e escrita. Tirar um alvo estreita o conjunto de
    alvos afirmativos da fase 3 -- nao reprova ligante nenhum: quem perder o
    unico alvo cai para a faixa 3 e aparece em `nao_ranqueados` com o motivo.
    """
    if cfg is None:
        return
    caminho = palvos.caminho_enriquecidos(cfg)
    if caminho is None:
        return
    df = palvos.ler_enriquecidos(caminho)
    if df is None:
        return
    st.markdown("**Lista de alvos enriquecidos** "
                f"— {len(df):,} alvos em `{caminho}`")
    st.caption("Editar aqui estreita ou amplia o que a fase 3 considera alvo "
               "afirmativo. Nao reprova ligante: quem ficar sem alvo cai para "
               "a faixa 3 e aparece em `nao_ranqueados` com o motivo.")
    editado = st.data_editor(df, num_rows="dynamic", key="editor_alvos",
                             height=300, width="stretch")
    mudou = palvos.diferenca(df, editado)
    if mudou["removidos"] or mudou["acrescentados"]:
        st.warning(
            f"a gravar: {len(mudou['removidos'])} a remover "
            f"({', '.join(mudou['removidos'][:5]) or '-'}), "
            f"{len(mudou['acrescentados'])} a acrescentar "
            f"({', '.join(mudou['acrescentados'][:5]) or '-'})")
    if st.button("Gravar a lista de alvos", key="gravar_alvos"):
        try:
            relato = palvos.gravar_enriquecidos(editado, caminho)
            st.success(
                f"{relato['total']} alvos gravados · "
                f"{len(relato['removidos'])} removidos · "
                f"{len(relato['acrescentados'])} acrescentados · "
                f"md5 `{(relato['md5_anterior'] or '')[:8]}` -> "
                f"`{relato['md5'][:8]}`")
        except ValueError as erro:
            st.error(str(erro))


def bloco_parametros(cfg):
    st.subheader("2. Parametros")
    st.caption("Cada numero vem com a medicao que o escolheu. Campo vazio "
               "devolve o parametro ao padrao da CLI.")
    if cfg is None:
        return
    mudancas: dict = {}
    for etapa in caminhos.ETAPAS:
        secao = cfg.get(etapa, {})
        with st.expander(f"{etapa} — {DESCRICAO[etapa][:70]}..."):
            formulario = _formulario_da_etapa(etapa, secao)
            mudancas[etapa] = parametros.somente_mudancas(secao, formulario)
    total = sum(len(m) for m in mudancas.values())
    if total:
        st.caption(f"{total} parametro(s) a gravar: " + ", ".join(
            f"{e}.{c}" for e, m in mudancas.items() for c in m))
    if st.button("Gravar configuracao", key="gravar_config", disabled=not total):
        try:
            config.gravar(mudancas)
            st.success(f"`json/config.json` atualizado ({total} parametro(s)); "
                       "as notas e as chaves intocadas ficaram como estavam.")
        except config.ErroDeConfig as erro:
            st.error(str(erro))


def _formulario_da_etapa(etapa: str, secao: dict) -> dict:
    """Um widget por chave conhecida, com o `_nota_*` como ajuda.

    A conversao de texto para valor fica em `mods/painel/parametros.py`: era o
    unico pedaco de logica que ficou dentro do Streamlit, e foi exatamente ali
    que nasceu o bug que separava `reprovar` caractere a caractere."""
    valores: dict = {}
    for chave, tipos in config.CHAVES[etapa].items():
        nota = secao.get(f"_nota_{chave}")
        atual = secao.get(chave)
        campo = f"p_{etapa}_{chave}"
        if bool in tipos:
            valores[chave] = st.checkbox(chave, value=bool(atual), key=campo,
                                         help=nota)
        elif int in tipos or float in tipos:
            inteiro = tipos == (int,)
            valores[chave] = st.number_input(
                chave, value=(int(atual) if inteiro else float(atual))
                if atual is not None else None,
                step=1 if inteiro else 0.05, key=campo, help=nota,
                # Formato explicito: o padrao de float e "%0.2f", e com ele
                # `corte: 0.125` voltaria 0.13 do navegador, em silencio.
                format=None if inteiro else "%g")
        elif list in tipos or str in tipos:
            rotulo = (f"{chave} (separados por vírgula)" if list in tipos
                      else chave)
            texto = st.text_input(rotulo, value=parametros.para_texto(atual),
                                  key=campo, help=nota)
            valores[chave] = parametros.de_texto(texto, anterior=atual)
        else:
            valores[chave] = st.text_input(
                chave, value=parametros.para_texto(atual), key=campo,
                help=nota) or None
    return valores


def bloco_o_que_rodar(cfg):
    st.subheader("3. O que rodar")
    emc = execucao.em_curso()
    escolhidas = st.multiselect(
        "Etapas", list(caminhos.ETAPAS), default=list(caminhos.ETAPAS),
        key="etapas_a_rodar",
        help="Uma etapa so, ou da etapa escolhida ate o fim. O orquestrador "
             "nao sabe dizer 'a primeira e a terceira'.")
    retomar = st.checkbox(
        "Retomar (pular etapa cuja saida ja confere com a entrada atual)",
        key="retomar",
        help="Sem isto, tudo roda. Pular por omissao entregaria resultado "
             "velho como novo.")
    if st.button("Rodar", key="rodar", disabled=emc is not None,
                 type="primary"):
        try:
            iniciada = execucao.iniciar(escolhidas, retomar=retomar)
            st.success(f"execucao iniciada (PID {iniciada['pid']})")
        except (ValueError, execucao.JaRodando) as erro:
            st.error(str(erro))
    if emc is not None:
        st.info(f"ja ha execucao de pe: PID {emc['pid']}, "
                f"{', '.join(emc['etapas']) or 'todas as etapas'}, "
                f"ha {_duracao(time.time() - emc['inicio'])}")
        if st.button("Encerrar", key="encerrar"):
            execucao.encerrar()
            st.rerun()


@st.fragment(run_every=2)
def bloco_acompanhamento():
    st.subheader("4. Em curso")
    # Colher o filho ANTES de perguntar se ha execucao: o `main.py` e filho do
    # processo do Streamlit, e um filho que saiu e ninguem colheu fica zumbi --
    # `os.kill(pid, 0)` responde sucesso, e o painel declararia de pe toda
    # execucao que terminou. E tambem o unico lugar de onde o codigo de saida
    # chega a tela.
    fim = execucao.colher()
    if fim is not None:
        _relatar_fim(fim)
    emc = execucao.em_curso()
    if emc is None:
        ultima = execucao.ultima_execucao()
        if ultima is None:
            st.caption("nenhuma execucao registrada por este painel")
        else:
            _relatar_fim(ultima)
    else:
        for etapa in emc["etapas"] or list(caminhos.ETAPAS):
            estado = leitura.estado_da_etapa(etapa)
            if estado["estado"] != "parcial":
                continue
            # A contagem do ritmo vem do TEMPORARIO inteiro, nao da copia: a
            # copia tem teto, e um contador que para no teto faz a taxa
            # observada nunca existir.
            tmp = leitura.temporario_da_etapa(etapa)
            escritas = leitura.contar_linhas(tmp) if tmp else None
            if escritas is not None:
                execucao.registrar_ritmo(escritas, etapa=etapa)
            taxa = execucao.ritmo_por_hora(etapa=etapa)
            st.write(f"**{etapa}**: {escritas:,} linhas escritas"
                     if escritas is not None else f"**{etapa}**: escrevendo")
            if taxa:
                st.caption(f"ritmo observado: {taxa:,.0f} linhas/h "
                           "(taxa medida, nao estimativa)")
            else:
                st.caption("ritmo: ainda sem duas amostras distantes o "
                           "bastante (leva ~10 min)")
            copia = caminhos.EXECUCAO / f"parcial_{etapa}.tsv"
            if leitura.copiar_parcial(etapa, copia, limite=200) is not None:
                try:
                    df = leitura.ler_tabela(copia, limite=200)
                except leitura.TabelaInvalida:
                    df = None  # o parcial ainda nao tem nem o cabecalho
                if df is not None and len(df):
                    st.dataframe(df, height=240, width="stretch")
    log = execucao.cauda_do_log()
    if log:
        with st.expander("log do orquestrador (fim)", expanded=emc is not None):
            st.code(log, language="text")


def _relatar_fim(ultima: dict):
    """O resultado da ultima execucao, pelo codigo de saida do `main.py`.

    `codigo` nulo significa que nao deu para saber -- painel reiniciado, e o
    processo ja colhido pelo sistema. A tela diz isso, em vez de mostrar 0 e
    chamar de completa uma execucao que pode ter falhado."""
    codigos = {0: ("success", "completa"),
               3: ("warning", "parcial (linhas de entrada descartadas)"),
               1: ("error", "nenhuma entrada utilizavel"),
               2: ("error", "erro (nada publicado)")}
    codigo = ultima.get("codigo")
    duracao = (f", {_duracao(ultima['fim'] - ultima['inicio'])}"
               if ultima.get("inicio") else "")
    etapas = ", ".join(ultima.get("etapas") or []) or "todas as etapas"
    if codigo is None:
        st.info(f"ultima execucao: terminou, e o codigo de saida nao foi "
                f"registrado (o painel foi reiniciado){duracao} — {etapas}. "
                f"O log e os manifestos de `tsv/` dizem o que foi publicado.")
        return
    tipo, texto = codigos.get(codigo, ("info", "codigo desconhecido"))
    getattr(st, tipo)(
        f"ultima execucao: {texto} (codigo {codigo}){duracao} — {etapas}")


def aba_operar():
    cfg = _cfg()
    bloco_entradas()
    bloco_alvos_enriquecidos(cfg)
    st.divider()
    bloco_parametros(cfg)
    st.divider()
    bloco_o_que_rodar(cfg)
    st.divider()
    bloco_acompanhamento()


# =========================================================================
# Analisar
# =========================================================================
TELAS = {
    "Ranking": ("ranking", "ranking.tsv"),
    "Nao ranqueados": ("ranking", "nao_ranqueados.tsv"),
    "Portao": ("portao", "portao.tsv"),
    "Propriedades": ("portao", "propriedades.tsv"),
    "Bula": ("portao", "alertas_bula.tsv"),
    "Alvos": ("portao", "alvos_anotados.tsv"),
    "Evidencias": ("etapa1", "evidencias.tsv"),
}


def _selo_do_manifesto(etapa: str):
    """De onde a tabela veio. Sem manifesto, aviso -- e nao a tabela
    escondida: tabela sem manifesto nao e citavel, e isso se diz, nao se
    esconde."""
    manifesto = leitura.estado_da_etapa(etapa)["manifesto"]
    if manifesto is None:
        st.warning(
            f"`tsv/{etapa}/manifesto.json` ausente ou ilegivel: esta tabela "
            "NAO e citavel. Rode a etapa de novo para publicar os dois juntos.")
        return
    versoes = manifesto.get("versoes", {})
    md5 = str(versoes.get("ligantes_md5", ""))[:8]
    with st.expander("procedencia desta tabela", expanded=False):
        st.caption(f"ligantes md5 `{md5}` · "
                   f"{manifesto.get('quando', 'data nao registrada')} · "
                   f"dispositivo {manifesto.get('dispositivo', 'nao registrado')}")
        st.json(manifesto, expanded=False)


def _tabela_ou_aviso(etapa: str, arquivo: str):
    """Devolve o DataFrame, ou `None` tendo dito na tela por que nao ha.

    Quatro estados, cada um com a sua mensagem: a etapa nao rodou; esta
    rodando AGORA (e a tabela, se existe, e a anterior); a tabela nao foi
    gerada nesta execucao; e a tabela existe mas nao e um TSV do projeto.
    Nenhum deles e um traceback."""
    estado = leitura.estado_da_etapa(etapa)
    rodando = estado["estado"] == "parcial"
    if estado["estado"] == "ausente":
        st.info(f"a etapa `{etapa}` ainda nao foi rodada: nao ha "
                f"`tsv/{etapa}/{arquivo}`. Rode pela aba Operar.")
        return None
    caminho = caminhos.TSV / etapa / arquivo
    if not caminho.is_file():
        if rodando:
            st.info(f"a etapa `{etapa}` esta **rodando agora**: "
                    f"`{arquivo}` so e publicada quando ela termina inteira. "
                    f"O andamento esta na aba Operar.")
        else:
            st.info(f"`tsv/{etapa}/{arquivo}` nao existe nesta execucao "
                    f"(a etapa pode ter rodado sem a opcao que o gera).")
        return None
    try:
        df = leitura.ler_tabela(caminho)
    except leitura.TabelaInvalida as erro:
        st.error(f"`tsv/{etapa}/{arquivo}` nao esta no formato do projeto "
                 f"(cabecalho `#` e uma linha por registro): {erro}")
        return None
    if rodando:
        st.warning(
            f"a etapa `{etapa}` esta rodando agora: a tabela abaixo e a "
            f"**anterior**, e vai ser substituida ao fim da execucao. O "
            f"manifesto ao lado e o dela.")
    _selo_do_manifesto(etapa)
    if df is None or not len(df):
        st.info(f"`{arquivo}` esta vazia: a etapa rodou e nenhuma linha "
                f"qualificou. E resultado, nao falha.")
        return None
    linhas = leitura.contar_linhas(caminho) or 0
    if linhas > len(df):
        st.caption(f"mostrando as primeiras {len(df):,} de {linhas:,} linhas; "
                   f"o TSV no disco esta inteiro")
    return df


def tela_ranking():
    df = _tabela_ou_aviso("ranking", "ranking.tsv")
    if df is None:
        return
    st.dataframe(vistas.resumo_das_faixas(df), hide_index=True,
                 width="stretch")
    colunas = st.columns(2)
    with colunas[0]:
        faixas = st.multiselect("Faixa", [f for f in FAIXAS if f != "3"],
                                key="f_faixa")
    with colunas[1]:
        alvos = st.multiselect(
            "Alvo", sorted(df["alvo"].dropna().unique())
            if "alvo" in df.columns else [], key="f_alvo")
    recorte = vistas.filtrar_alvo(vistas.filtrar_faixa(df, faixas), alvos)
    st.caption(f"{len(recorte):,} linha(s)")
    st.dataframe(recorte, hide_index=True, height=480, width="stretch")


def tela_nao_ranqueados():
    df = _tabela_ou_aviso("ranking", "nao_ranqueados.tsv")
    if df is None:
        return
    st.caption(f"`{MOTIVO_FAIXA_3}` e a faixa 3: "
               f"{vistas.legenda_da_faixa('3')}")
    st.dataframe(vistas.motivos_dos_nao_ranqueados(df), hide_index=True,
                 width="stretch")
    st.dataframe(df, hide_index=True, height=400, width="stretch")


def tela_portao():
    df = _tabela_ou_aviso("portao", "portao.tsv")
    if df is None:
        return
    st.caption("`indeterminado` nao e aprovacao: e criterio que nao foi "
               "pedido ou nao pode ser avaliado. Ausencia de documento nao e "
               "ausencia de risco.")
    st.dataframe(vistas.portao_para_tela(df), hide_index=True, height=480,
                 width="stretch")


def tela_bula():
    df = _tabela_ou_aviso("portao", "alertas_bula.tsv")
    if df is None:
        return
    st.caption("Um achado por linha, com o trecho literal. A negacao em texto "
               "livre ('no evidence of neuropathy') NAO esta resolvida: e por "
               "isso que o trecho fica a vista. `n_bulas_consultadas` e a "
               "janela; `n_bulas_com_o_termo` e o acervo inteiro.")
    st.dataframe(vistas.alertas_por_forca(df), hide_index=True, height=480,
                 width="stretch")


def tela_alvos():
    df = _tabela_ou_aviso("portao", "alvos_anotados.tsv")
    if df is None:
        return
    st.caption("Expressao e propriedade do ALVO: alvo de expressao baixa nao "
               "reprova ligante nenhum, porque o ligante pode ter outros alvos.")
    st.dataframe(vistas.alvos_com_expressao(df), hide_index=True, height=480,
                 width="stretch")


def tela_simples(etapa: str, arquivo: str):
    df = _tabela_ou_aviso(etapa, arquivo)
    if df is not None:
        st.dataframe(df, hide_index=True, height=480, width="stretch")


def aba_analisar():
    st.caption("Esta aba so le. O que ela mostra e o que esta publicado em "
               "`tsv/`, com o manifesto de cada etapa ao lado.")
    escolha = st.radio("Tabela", list(TELAS), horizontal=True, key="tela")
    if escolha == "Ranking":
        tela_ranking()
    elif escolha == "Nao ranqueados":
        tela_nao_ranqueados()
    elif escolha == "Portao":
        tela_portao()
    elif escolha == "Bula":
        tela_bula()
    elif escolha == "Alvos":
        tela_alvos()
    else:
        tela_simples(*TELAS[escolha])


def barra_lateral():
    """A estrutura do pipeline, **lida do proprio codigo**.

    As etapas vem da tabela `ETAPAS` do `main.py`, as faixas de
    `mods/ranque/faixa.py`, as pastas de `mods/caminhos.py`. Texto decorado
    aqui divergiria do pipeline na primeira mudanca, e documentacao que mente
    e pior que nenhuma."""
    with st.sidebar:
        st.markdown("## Estrutura do pipeline")
        st.caption(f"commit `{estrutura.versao_do_pipeline()}` · "
                   f"as informacoes abaixo sao lidas do codigo, nao escritas "
                   f"a mao")

        st.markdown("### As tres etapas")
        st.markdown(
            "Cada etapa roda como **processo separado** e publica a propria "
            "tabela com o proprio manifesto ao lado. Se a fase 3 falhar, a "
            "etapa 2 continua publicada e citavel — e por isso que o "
            "orquestrador chama por subprocesso.")
        for i, etapa in enumerate(estrutura.etapas(), start=1):
            rotulo = "fase 3" if etapa["id"] == "ranking" else f"etapa {i}"
            st.markdown(f"**{rotulo} · `{etapa['id']}`**")
            st.caption(etapa["descricao"])
            st.markdown(f"- roda: `{etapa['script']}`\n"
                        f"- escreve: `{etapa['saida']}`\n"
                        f"- com: `{etapa['manifesto']}`")

        st.markdown("### As faixas de evidencia")
        st.markdown(
            "A ordem e **lexicografica**: classe de evidencia, depois "
            "procedencia do alvo, depois valor. Medido vem antes de previsto, "
            "e alvo que voce declarou antes de alvo do enriquecimento — e por "
            "isso **coletar largo nao e ranquear largo**.")
        for faixa in estrutura.faixas():
            marca = "" if faixa["ranqueia"] else " · **nao ranqueia**"
            st.markdown(f"**`{faixa['faixa']}`**{marca} — {faixa['significado']}")

        st.markdown("### A arvore")
        st.markdown(
            "A regra que desempata `input/` de `tsv/`, porque ligantes e "
            "alvos tambem sao TSV: **procedencia, nao extensao.**")
        for pasta in estrutura.pastas():
            st.markdown(f"- **`{pasta['nome']}`** — {pasta['papel']}")

        st.markdown("### Codigos de saida")
        st.markdown(
            "O orquestrador devolve o **pior** entre as etapas que rodou. O 3 "
            "existe porque automacao nao le o `stderr`.")
        for codigo, texto in sorted(estrutura.codigos_de_saida().items()):
            st.markdown(f"- **{codigo}** — {texto}")

        st.markdown("### Duas regras que atravessam tudo")
        st.markdown(
            "1. **`indeterminado` nunca vira `passa`.** Criterio nao pedido ou "
            "que nao pode ser avaliado sai `indeterminado`, e ausencia de "
            "documento nao e ausencia de risco.\n"
            "2. **Tabela sem manifesto nao e citavel.** As saidas de uma etapa "
            "sao preparadas em temporarios e publicadas juntas, por troca de "
            "nome: uma falha no meio nao deixa nem o resultado velho nem o novo "
            "pela metade.")

        st.markdown("### Dependencias")
        # O painel NAO instala: ele e uma interface, e instalar por clique
        # mexe no ambiente de quem abriu a aba. Ele mostra o estado e o
        # comando -- quem instala e o orquestrador, na partida, ou a pessoa.
        estado = dependencias.conferir()
        if estado["completo"]:
            st.caption("exigidas: todas presentes — "
                       + ", ".join(f"`{d}`" for d in estrutura.dependencias()))
        else:
            st.error("faltam dependencias exigidas: "
                     + ", ".join(d["requisito"] for d in estado["faltando"])
                     + ". Rode `pip install -r requirements.txt`")
        for item in estado["opcionais_ausentes"]:
            st.warning(f"`{item['pacote']}` nao esta instalado (opcional, "
                       f"pesado). Sem ele, desliga {item['desliga']}")


st.title("Pipeline druglikeness")
barra_lateral()
operar, analisar = st.tabs(["Operar", "Analisar"])
with operar:
    aba_operar()
with analisar:
    aba_analisar()
