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
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
import streamlit as st

from mods.evidencias.entrada import ler_alvos, ler_ligantes, ler_sintomas
from mods.evidencias.esquema import desescapar, escapar, ler_evidencias
from mods.evidencias.open_targets import COLUNAS_ENRIQUECIDOS, buscar_efo
from mods.evidencias.resumo import (
    montar_comando, pares_previstos, resumo_alvos, resumo_ligantes)

RAIZ = Path(__file__).resolve().parent
LINHAS_NA_TELA = 200_000

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
    linhas = Path(caminho).read_text(encoding="utf-8").splitlines()
    colunas = linhas[0].lstrip("#").split("\t")
    return pd.DataFrame([[desescapar(c) for c in l.split("\t")]
                         for l in linhas[1:] if l], columns=colunas)


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
            t["proc"].terminate()
    else:
        if "concluido" not in t:
            t["concluido"] = codigo
            st.rerun(scope="app")  # redesenha o app com o resultado
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
    efos: dict[str, str] = {}
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
    if st.button("Enriquecer alvos", disabled=not efos or ocupado):
        saida_e = pasta_trabalho / "enriquecimento"
        iniciar("enriquecer", montar_comando(
            "enriquecer", python=sys.executable, ligantes=str(caminho_lig),
            alvos=str(caminho_alv), sintomas=str(caminho_sin),
            saida=str(saida_e / "evidencias.tsv"), rotulos=rotulos_caminho,
            alvos_enriquecidos=str(saida_e / "alvos_enriquecidos.tsv"), efos=efos,
            cache=cache_pasta, corte_associacao=corte_assoc, teto_via=int(teto_via),
            limite=int(limite) or None), saida_e, saida_e)
        st.rerun()

arq_enriq = pasta_trabalho / "enriquecimento" / "alvos_enriquecidos.tsv"
revisados = pasta_trabalho / "alvos_revisados.tsv"
if arq_enriq.exists():
    st.subheader("Lista proposta (edite à vontade)")
    editado = st.data_editor(tabela_de_alvos(arq_enriq), num_rows="dynamic",
                             hide_index=True, key="editor_alvos")
    if st.button("Usar esta lista na coleta"):
        gravar_tabela_de_alvos(editado, revisados)
        st.session_state["alvos_da_coleta"] = str(revisados)
        st.success(f"{len(editado)} alvos gravados em {revisados}")

# ---------------------------------------------------------------- 3. coletar
st.header("3. Coletar evidências")
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
if st.button("Coletar evidências", type="primary", disabled=not pode or ocupado):
    saida_c = pasta_trabalho / "coleta"
    iniciar("coletar", montar_comando(
        "coletar", python=sys.executable, ligantes=str(caminho_lig),
        alvos=alvos_da_coleta, sintomas=str(caminho_sin),
        saida=str(saida_c / "evidencias.tsv"), rotulos=rotulos_caminho,
        checkpoint=pesos_caminho if pesos_ok else None,
        evidencias_chembl=medido_caminho or None, cache=cache_pasta, corte=corte,
        limite=int(limite) or None), saida_c / "evidencias.tsv", saida_c)
    st.rerun()
if not pode:
    st.caption("Sem pesos e sem `evidencias_chembl.csv` não há o que coletar.")

acompanhar()

# ------------------------------------------------------------------ resultado
evidencias = pasta_trabalho / "coleta" / "evidencias.tsv"
if evidencias.exists():
    st.header("Resultado")
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
    d1.download_button("Baixar evidencias.tsv", evidencias.read_bytes(), "evidencias.tsv")
    man = pasta_trabalho / "coleta" / "manifesto.json"
    if man.exists():
        d2.download_button("Baixar manifesto.json", man.read_bytes(), "manifesto.json")
        with st.expander("Manifesto"):
            st.json(json.loads(man.read_text(encoding="utf-8")))
