# -*- coding: utf-8 -*-
"""
painel.py -- painel Streamlit da etapa 1 (coletar evidencias).

    streamlit run painel.py

O painel e so uma casca: valida as entradas no upload e dispara a MESMA CLI que
o artigo cita (`python -m mods.coletar_evidencias`), em segundo plano. Nada de
coleta vive aqui; tudo o que o painel faz tambem se faz pela linha de comando.
"""

import hashlib
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
import streamlit as st

from mods.evidencias.entrada import ler_alvos, ler_ligantes, ler_sintomas
from mods.evidencias.esquema import escapar, ler_evidencias
from mods.evidencias.open_targets import COLUNAS_ENRIQUECIDOS, buscar_efo
from mods.evidencias.resumo import (
    assinatura_de_arquivo, chave_coleta, chave_do_editor, chave_enriquecimento,
    edicao_pendente, ler_tabela_de_alvos, montar_comando, pares_previstos,
    resumo_alvos, resumo_ligantes)

RAIZ = Path(__file__).resolve().parent
# Teto do que a tabela da tela le. O arquivo real pode ter ~23 milhoes de
# linhas: `ler_evidencias` le em fluxo e a leitura PARA aqui, em vez de
# materializar tudo e cortar depois. Quem quer o resto baixa o arquivo.
LINHAS_NA_TELA = 5_000

st.set_page_config(page_title="Evidências — etapa 1", layout="wide")


# ---------------------------------------------------------------- utilitarios
def salvar_envio(envio, pasta: Path) -> Path:
    pasta.mkdir(parents=True, exist_ok=True)
    destino = pasta / Path(envio.name).name
    dados = envio.getvalue()
    if not destino.exists() or hashlib.md5(destino.read_bytes()).digest() != \
            hashlib.md5(dados).digest():
        destino.write_bytes(dados)
    return destino


@st.cache_data(show_spinner="Lendo e validando...")
def ler_tudo(caminho_ligantes, caminho_alvos, caminho_sintomas, chave):
    # `chave` (MD5 dos tres arquivos) invalida o cache quando o conteudo muda.
    ligantes, p_lig = ler_ligantes(caminho_ligantes)
    alvos, p_alvo = ler_alvos(caminho_alvos)
    return ligantes, p_lig, alvos, p_alvo, ler_sintomas(caminho_sintomas)


@st.cache_data(show_spinner="Lendo os rótulos do modelo...")
def ler_rotulos(caminho, mtime):
    from mods.prever_alvos import carregar_rotulos
    return carregar_rotulos(caminho)


def tabela_de_alvos(caminho) -> pd.DataFrame:
    colunas, linhas = ler_tabela_de_alvos(caminho)
    return pd.DataFrame(linhas, columns=colunas)


def gravar_tabela_de_alvos(df: pd.DataFrame, caminho: Path) -> None:
    colunas = [c for c in COLUNAS_ENRIQUECIDOS if c in df.columns] + \
        [c for c in df.columns if c not in COLUNAS_ENRIQUECIDOS]
    with open(caminho, "w", encoding="utf-8", newline="\n") as f:
        f.write("#" + "\t".join(colunas) + "\n")
        for _, l in df.iterrows():
            f.write("\t".join(escapar("" if pd.isna(l[c]) else l[c])
                              for c in colunas) + "\n")


def iniciar(tipo: str, comando: list[str], saida: Path, pasta: Path) -> None:
    pasta.mkdir(parents=True, exist_ok=True)
    log = pasta / f"{tipo}.log"
    with open(log, "w", encoding="utf-8") as f:
        f.write("$ " + " ".join(comando) + "\n")
    # Processo proprio: sobrevive aos reruns do Streamlit e le/grava so em disco.
    proc = subprocess.Popen(
        comando, cwd=RAIZ, stdout=open(log, "a"), stderr=subprocess.STDOUT,
        start_new_session=True)
    st.session_state["trabalho"] = {
        "tipo": tipo, "proc": proc, "log": log, "inicio": time.time(),
        "saida": saida}


def edicao_em_curso() -> bool:
    """Ha edicao nao aplicada no editor QUE ESTA NA TELA. So o editor atual
    conta: estado de um editor de outra lista (ou de outra geracao, depois de
    um descarte) nao e trabalho que o pesquisador possa perder."""
    chave = st.session_state.get("editor_atual")
    return bool(chave) and edicao_pendente(st.session_state.get(chave))


@st.fragment(run_every=2)
def acompanhar():
    t = st.session_state.get("trabalho")
    if not t:
        return
    codigo = t["proc"].poll()
    decorrido = int(time.time() - t["inicio"])
    nome = "Enriquecimento" if t["tipo"] == "enriquecer" else "Coleta de evidências"
    texto = t["log"].read_text(encoding="utf-8", errors="replace")
    if codigo is None:
        st.info(f"{nome} em andamento há {decorrido // 60} min {decorrido % 60} s. "
                "A CLI não informa percentual; o relatório abaixo é o que ela "
                "escreveu até agora. Um enriquecimento sem cache leva 23 a 25 min.")
        st.progress(min(decorrido / 1500, 0.95) if t["tipo"] == "enriquecer" else 0.0,
                    text="tempo decorrido (estimativa, não progresso real)"
                    if t["tipo"] == "enriquecer" else "sem estimativa de duração")
        if st.button("Cancelar", key="cancelar"):
            try:  # o grupo todo: a previsao pode ter processos filhos
                os.killpg(t["proc"].pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
    else:
        if "concluido" not in t:
            t["concluido"] = codigo
            if edicao_em_curso():
                # Redesenho suprimido para não perder a edição em curso. Fica
                # registrado: o aviso e o botão de descarte só fazem sentido
                # quando um redesenho foi REALMENTE suprimido, e antes eles
                # apareciam a cada poll depois do fim do trabalho.
                t["rerun_suprimido"] = True
            else:
                st.rerun(scope="app")  # redesenha o app com o resultado
        if t.get("rerun_suprimido") and edicao_em_curso():
            st.warning("Execução terminada, mas há edições não aplicadas na lista "
                       "de alvos: a página não foi redesenhada para não as perder.")
            if st.button("Mostrar o resultado (descarta o que não foi aplicado)"):
                # Geração nova = widget novo, com estado vazio. Não depende de
                # o Streamlit esquecer o estado do widget anterior.
                st.session_state["geracao_editor"] = \
                    st.session_state.get("geracao_editor", 0) + 1
                st.rerun(scope="app")
        (st.success if codigo == 0 else st.error)(
            f"{nome} {'concluído' if codigo == 0 else f'falhou (código {codigo})'} "
            f"em {decorrido // 60} min {decorrido % 60} s.")
    st.code(texto[-6000:] or "(sem saída ainda)", language="text")


# ------------------------------------------------------------------ barra lateral
with st.sidebar:
    st.header("Modelo e pastas")
    rotulos_caminho = st.text_input(
        "target_labels.tsv (obrigatório)", os.environ.get("PP4CKD_ROTULOS", ""),
        help="Variável de ambiente PP4CKD_ROTULOS como padrão. Sem este arquivo o "
             "painel não sabe quantos alvos existem no modelo.")
    pesos_caminho = st.text_input(
        "Pesos do PP4CKD, .pt (opcional)", os.environ.get("PP4CKD_PESOS", ""),
        help="Variável PP4CKD_PESOS. O pacote publicado no Zenodo está truncado; "
             "use ecfp4_fold_4.pt (ver README).")
    medido_caminho = st.text_input(
        "evidencias_chembl.csv (opcional)", os.environ.get("EVIDENCIAS_CHEMBL", ""),
        help="Variável EVIDENCIAS_CHEMBL. Liga a evidência medida.")
    pasta_trabalho = Path(st.text_input("Pasta de trabalho", "painel_saida")).expanduser()
    cache_pasta = st.text_input("Cache de rede", str(pasta_trabalho / "cache_rede"))
    st.header("Parâmetros")
    corte = st.number_input("Corte de probabilidade", 0.0, 1.0, 0.1, 0.05)
    corte_assoc = st.number_input("Corte de associação (Open Targets)", 0.0, 1.0, 0.4, 0.05)
    teto_via = st.number_input("Teto de proteínas por via (Reactome)", 10, 1000, 150, 10)
    limite = st.number_input("Usar só os N primeiros ligantes (0 = todos)", 0, 10**7, 0)

rotulos = None
erro_rotulos = None
if rotulos_caminho:
    try:
        arq = Path(rotulos_caminho).expanduser()
        rotulos = ler_rotulos(str(arq), arq.stat().st_mtime if arq.exists() else 0)
    except (OSError, RuntimeError, ValueError) as e:
        erro_rotulos = str(e)
pesos_ok = bool(pesos_caminho) and Path(pesos_caminho).expanduser().is_file()

# ------------------------------------------------------------------------- corpo
st.title("Etapa 1 — coletar evidências")
st.caption("O painel dispara a mesma CLI que o artigo cita "
           "(`python -m mods.coletar_evidencias`); tudo o que ele faz também se faz "
           "pela linha de comando.")

if rotulos is None:
    st.error(erro_rotulos or "Informe o caminho do `target_labels.tsv` na barra "
             "lateral: o painel precisa dele para dizer quantos alvos existem no modelo.")
    st.stop()
st.sidebar.success(f"Modelo: {len(rotulos)} alvos")

if pesos_ok and medido_caminho:
    st.success("Evidência prevista (PP4CKD) e evidência medida (ChEMBL) ligadas.")
elif pesos_ok:
    st.info("Só evidência prevista: sem `evidencias_chembl.csv` não há evidência medida.")
elif medido_caminho:
    st.warning("**Modo degradado: só evidência medida.** Sem os pesos do PP4CKD "
               "(barra lateral) não há probabilidade prevista; as linhas virão só "
               "do ChEMBL.")
else:
    st.warning("**Modo degradado, sem nenhuma fonte de evidência.** Informe os pesos "
               "e/ou o `evidencias_chembl.csv` na barra lateral para coletar. "
               "O enriquecimento de alvos funciona assim mesmo.")

st.header("1. Entradas")
c1, c2, c3 = st.columns(3)
envio_lig = c1.file_uploader("Ligantes (TSV: nome, smiles)", type=["tsv", "csv", "txt"])
envio_alv = c2.file_uploader("Alvos (TSV: alvo_chembl_id, ou alvo + organismo)",
                             type=["tsv", "csv", "txt"])
envio_sin = c3.file_uploader("Sintomas (um por linha)", type=["txt", "tsv"])
if not (envio_lig and envio_alv and envio_sin):
    st.info("Envie os três arquivos; a validação roda em seguida, antes de qualquer execução.")
    st.stop()

entrada = pasta_trabalho / "entrada"
caminho_lig = salvar_envio(envio_lig, entrada / "ligantes")
caminho_alv = salvar_envio(envio_alv, entrada / "alvos")
caminho_sin = salvar_envio(envio_sin, entrada / "sintomas")
chave = hashlib.md5(b"".join(p.read_bytes() for p in (caminho_lig, caminho_alv, caminho_sin))
                    ).hexdigest()
ligantes, p_lig, alvos, p_alvo, sintomas = ler_tudo(
    str(caminho_lig), str(caminho_alv), str(caminho_sin), chave)
if limite:
    ligantes = ligantes[:int(limite)]
r_lig = resumo_ligantes(ligantes, p_lig)
r_alv = resumo_alvos(alvos, rotulos, problemas=p_alvo)

m = st.columns(5)
m[0].metric("Ligantes válidos", r_lig["validos"])
m[1].metric("SMILES ilegíveis", r_lig["invalidos"])
m[2].metric("Alvos no modelo", f"{r_alv['no_modelo']} de {r_alv['total']}")
m[3].metric("Alvos com pouco treino", r_alv["pouco_treino"])
m[4].metric("Sintomas", len(sintomas))

if r_alv["sem_organismo"]:
    st.error(f"{r_alv['sem_organismo']} alvo(s) dado(s) por nome **sem organismo** serão "
             "descartados: o mesmo nome existe em humano, camundongo e rato.")
if r_alv["fora_do_modelo"]:
    st.warning(f"{r_alv['fora_do_modelo']} alvo(s) **não existem no modelo**: não terão "
               "evidência prevista (zero linhas não significa 'sem atividade').")
    st.dataframe(pd.DataFrame(r_alv["lista_fora_do_modelo"]), hide_index=True)
if r_alv["ambiguo"]:
    st.error(f"{r_alv['ambiguo']} alvo(s) dado(s) por nome casam com **mais de um** "
             "alvo do modelo: nenhum é escolhido (escolher seria um palpite entre "
             "alvos diferentes). Informe o `alvo_chembl_id`.")
    st.dataframe(pd.DataFrame(r_alv["lista_ambiguo"]), hide_index=True)
if r_alv["pouco_treino"]:
    st.warning(f"{r_alv['pouco_treino']} alvo(s) treinados com menos de "
               f"{r_alv['minimo_treino']} moléculas ativas: a previsão neles é frágil.")
    st.dataframe(pd.DataFrame(r_alv["lista_pouco_treino"])[
        ["alvo_chembl_id", "alvo", "organismo", "n_ativos"]], hide_index=True)
if r_lig["invalidos"]:
    st.error(f"{r_lig['invalidos']} ligante(s) com SMILES ausente ou ilegível serão descartados.")
if not sintomas:
    st.warning("Nenhum sintoma: não há enriquecimento possível, só os alvos informados serão usados.")
with st.expander(f"Problemas de leitura ({len(p_lig) + len(p_alvo)})"):
    st.dataframe(pd.DataFrame([{"arquivo": a, **p.__dict__}
                               for a, ps in (("ligantes", p_lig), ("alvos", p_alvo))
                               for p in ps]), hide_index=True)
if not ligantes or not alvos:
    st.error("Sem ligantes ou sem alvos utilizáveis: nada a executar.")
    st.stop()

# --------------------------------------------------------------- 2. enriquecer
st.header("2. Enriquecer alvos (opcional, passo separado)")
st.markdown(
    "O enriquecimento **propõe** mais alvos (Open Targets e Reactome). **Os alvos "
    "enriquecidos NÃO são usados na coleta até que você os devolva** no passo 3: "
    "revise e edite a lista antes de gastar a coleta. Um caso real deu 235 alvos "
    "contra 108.935 ligantes, cerca de 23 milhões de pares previstos.")

efos: dict[str, str] = {}
pasta_enr = pasta_trabalho / "enriquecimento"
ponteiro = pasta_enr / f"ultimo_{chave}.json"  # EFO da ultima execucao destes arquivos
if sintomas:
    st.subheader("Sintoma → EFO")
    st.caption("O Open Targets sugere condições; a escolha errada gera uma lista de "
               "alvos errada sem nenhum sinal. Confira as alternativas e escolha.")
    if st.button("Buscar EFO dos sintomas"):
        try:
            st.session_state["efos_candidatos"] = {s: buscar_efo(s) for s in sintomas}
        except Exception as e:  # rede/servico: mostrar, nao esconder
            st.error(f"Busca de EFO falhou: {e}")
    candidatos = st.session_state.get("efos_candidatos", {})
    for s in sintomas:
        cands = candidatos.get(s)
        if cands is None:
            continue
        opcoes = [c["id"] for c in cands]
        nomes = {c["id"]: f"{c['id']}: {c['name']}" for c in cands}
        escolha = st.selectbox(f"EFO para “{s}”", opcoes, index=None,
                               format_func=nomes.get, key=f"efo_{s}",
                               placeholder="escolha (sem escolha, não enriquece)")
        if escolha:
            efos[s] = escolha
    ocupado = bool(st.session_state.get("trabalho")) and \
        st.session_state["trabalho"]["proc"].poll() is None

    def pasta_de(efos_):
        return pasta_enr / chave_enriquecimento(
            chave, efos_, corte_assoc, int(teto_via), int(limite) or None)

    if st.button("Enriquecer alvos", disabled=not efos or ocupado):
        saida_e = pasta_de(efos)
        # Nova execucao: nada da anterior (lista ou revisao) sobrevive a ela.
        for velho in (saida_e / "alvos_enriquecidos.tsv", saida_e / "alvos_revisados.tsv"):
            velho.unlink(missing_ok=True)
        saida_e.mkdir(parents=True, exist_ok=True)
        ponteiro.write_text(json.dumps(efos), encoding="utf-8")
        st.session_state.pop("alvos_da_coleta", None)
        iniciar("enriquecer", montar_comando(
            "enriquecer", python=sys.executable, ligantes=str(caminho_lig),
            alvos=str(caminho_alv), sintomas=str(caminho_sin),
            saida=str(saida_e / "evidencias.tsv"), rotulos=rotulos_caminho,
            alvos_enriquecidos=str(saida_e / "alvos_enriquecidos.tsv"), efos=efos,
            cache=cache_pasta, corte_associacao=corte_assoc, teto_via=int(teto_via),
            limite=int(limite) or None), saida_e, saida_e)
        st.rerun()

    # A lista exibida e a destes arquivos, destes EFO e destes parametros: a
    # pasta e nomeada por eles. Sem EFO escolhido na tela (pagina recarregada),
    # vale o da ultima execucao registrada para estes mesmos arquivos.
    efos_vigentes = efos
    if not candidatos and ponteiro.exists():
        try:
            efos_vigentes = json.loads(ponteiro.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            efos_vigentes = {}
    pasta_vigente = pasta_de(efos_vigentes) if efos_vigentes else None
else:
    pasta_vigente = None

st.session_state["editor_atual"] = None  # redefinido abaixo se houver editor
arq_enriq = pasta_vigente / "alvos_enriquecidos.tsv" if pasta_vigente else None
revisados = pasta_vigente / "alvos_revisados.tsv" if pasta_vigente else None
if arq_enriq and arq_enriq.exists():
    st.subheader("Lista proposta (edite à vontade)")
    st.caption("Calculada para estes arquivos, com os EFO: "
               + "; ".join(f"{k} = {v}" for k, v in efos_vigentes.items()))
    tabela = tabela_de_alvos(arq_enriq)
    if tabela.empty:
        st.info("O enriquecimento não acrescentou nenhum alvo.")
    chave_editor = chave_do_editor(pasta_vigente.name,
                                   st.session_state.get("geracao_editor", 0))
    st.session_state["editor_atual"] = chave_editor
    editado = st.data_editor(tabela, num_rows="dynamic", hide_index=True,
                             key=chave_editor)
    b1, b2 = st.columns(2)
    if b1.button("Usar esta lista na coleta"):
        gravar_tabela_de_alvos(editado, revisados)
        st.session_state["alvos_da_coleta"] = str(revisados)
        st.success(f"{len(editado)} alvos gravados em {revisados}")
    if b2.button("Voltar à lista enviada"):
        revisados.unlink(missing_ok=True)
        st.session_state.pop("alvos_da_coleta", None)

# ---------------------------------------------------------------- 3. coletar
st.header("3. Coletar evidências")
# A escolha sobrevive ao recarregar a pagina: ela mora em disco, ao lado da lista
# de que veio (mesma pasta por arquivos + EFO + parametros).
if revisados is not None and revisados.exists():
    st.session_state["alvos_da_coleta"] = str(revisados)
else:
    st.session_state.pop("alvos_da_coleta", None)
alvos_da_coleta = st.session_state.get("alvos_da_coleta") or str(caminho_alv)
st.write("Alvos da coleta: **" + ("lista enriquecida e revisada" if
         alvos_da_coleta != str(caminho_alv) else "lista enviada, sem enriquecimento")
         + f"** (`{alvos_da_coleta}`)")
n_alvos_modelo = r_alv["no_modelo"]
if alvos_da_coleta != str(caminho_alv):
    a2, _ = ler_alvos(alvos_da_coleta)
    n_alvos_modelo = resumo_alvos(a2, rotulos)["no_modelo"]
if pesos_ok:
    st.write(f"Evidência prevista: cerca de **{pares_previstos(r_lig['validos'], n_alvos_modelo):,}** "
             f"pares ({r_lig['validos']:,} ligantes × {n_alvos_modelo} alvos no modelo).")
pode = pesos_ok or bool(medido_caminho)
ocupado = bool(st.session_state.get("trabalho")) and \
    st.session_state["trabalho"]["proc"].poll() is None
# A pasta da coleta e nomeada pelo que a determina (arquivos enviados, lista de
# alvos usada, corte, limite e fontes de evidencia), como a do enriquecimento:
# um resultado so pode pertencer as entradas que estao na tela.
try:
    assinatura_alvos = hashlib.md5(Path(alvos_da_coleta).read_bytes()).hexdigest()
except OSError:
    assinatura_alvos = "ausente"
pasta_coleta = pasta_trabalho / "coleta" / chave_coleta(
    chave, assinatura_alvos, float(corte), int(limite) or None,
    assinatura_de_arquivo(pesos_caminho if pesos_ok else ""),
    assinatura_de_arquivo(medido_caminho))
if st.button("Coletar evidências", type="primary", disabled=not pode or ocupado):
    # Execucao nova: o resultado anterior DESTA mesma chave sai da frente, para
    # que uma coleta cancelada ou falha nao mostre a tabela da anterior como se
    # fosse a sua (a escrita atomica da CLI deixa o arquivo velho no lugar).
    for velho in (pasta_coleta / "evidencias.tsv", pasta_coleta / "manifesto.json"):
        velho.unlink(missing_ok=True)
    iniciar("coletar", montar_comando(
        "coletar", python=sys.executable, ligantes=str(caminho_lig),
        alvos=alvos_da_coleta, sintomas=str(caminho_sin),
        saida=str(pasta_coleta / "evidencias.tsv"), rotulos=rotulos_caminho,
        checkpoint=pesos_caminho if pesos_ok else None,
        evidencias_chembl=medido_caminho or None, cache=cache_pasta, corte=corte,
        limite=int(limite) or None), pasta_coleta / "evidencias.tsv", pasta_coleta)
    st.rerun()
if not pode:
    st.caption("Sem pesos e sem `evidencias_chembl.csv` não há o que coletar.")

acompanhar()

# ------------------------------------------------------------------ resultado
# So o resultado DESTAS entradas e DESTES parametros: a pasta e nomeada por
# eles. Resultado de outra escolha nunca aparece como se fosse o atual.
evidencias = pasta_coleta / "evidencias.tsv"
man = pasta_coleta / "manifesto.json"
if evidencias.exists():
    st.header("Resultado")
    manifesto_lido = None
    if man.exists():
        try:
            manifesto_lido = json.loads(man.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            manifesto_lido = None
    if manifesto_lido is None:
        st.warning("Tabela sem `manifesto.json` ao lado: o resultado não é "
                   "citável e pode estar incompleto.")
    else:
        contagens = manifesto_lido.get("contagens", {})
        st.caption(
            f"Execução de {manifesto_lido.get('versoes', {}).get('data', '?')} "
            f"sobre os arquivos e parâmetros desta tela ({pasta_coleta.name}); "
            f"{contagens.get('linhas_evidencia', '?')} linhas.")
        if contagens.get("execucao_parcial"):
            st.warning("**Execução PARCIAL:** o manifesto registra linhas de "
                       "entrada descartadas ou falhas de item. A tabela abaixo "
                       "está incompleta; veja `falhas_de_item` no manifesto.")
    linhas = []
    for i, l in enumerate(ler_evidencias(evidencias)):
        if i >= LINHAS_NA_TELA:
            st.warning(f"Mostrando só as primeiras {LINHAS_NA_TELA:,} linhas; "
                       "baixe o arquivo para ver tudo.")
            break
        linhas.append(l)
    df = pd.DataFrame(linhas)
    if df.empty:
        st.info("A tabela de evidências está vazia (só o cabeçalho).")
    else:
        f1, f2, f3 = st.columns(3)
        tipos = f1.multiselect("Tipo de evidência", sorted(df["tipo_evidencia"].unique()))
        lig = f2.text_input("Ligante contém")
        alv = f3.text_input("Alvo contém")
        if tipos:
            df = df[df["tipo_evidencia"].isin(tipos)]
        if lig:
            df = df[df["ligante"].str.contains(lig, case=False, regex=False)]
        if alv:
            df = df[df["alvo"].str.contains(alv, case=False, regex=False)]
        st.caption(f"{len(df):,} linhas. Medido e predito nunca se misturam: "
                   "`pchembl` é medida, probabilidade é do modelo.")
        st.dataframe(df, hide_index=True)
    d1, d2 = st.columns(2)
    # Download diferido: o arquivo e aberto quando o usuario clica, nao a cada
    # rerun. `read_bytes()` aqui punha a tabela inteira (ate ~23 milhoes de
    # linhas) na memoria do servidor em todo redesenho da pagina.
    d1.download_button("Baixar evidencias.tsv", lambda: open(evidencias, "rb"),
                       "evidencias.tsv")
    if manifesto_lido is not None:
        d2.download_button("Baixar manifesto.json", lambda: open(man, "rb"),
                           "manifesto.json")
        with st.expander("Manifesto"):
            st.json(manifesto_lido)
