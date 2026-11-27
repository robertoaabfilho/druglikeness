# -*- coding: utf-8 -*-
"""
cache.py

Cache em disco das respostas de rede (ChEMBL, Open Targets, Reactome) e a
identificacao das versoes que elas representam.

Regras
------
* A chave e deterministica: sha256 de (url, corpo GraphQL) em JSON canonico.
  Cada entrada guarda a url e o corpo, e a leitura confere os dois: arquivo
  adulterado ou colisao vira falha, nunca resposta errada.
* Guarda a resposta CRUA da API, nunca algo calculado com parametros do
  usuario (corte, teto de via...). Mudar um parametro nao deixa nada velho
  para tras; so mudar o servico poderia, e para isso ver abaixo.
* NUNCA guarda: excecao (timeout, HTTP de erro, DNS: propagam sem gravar),
  resposta com sinal de erro (`errors`/`error_message`) ou pagina truncada
  (`page_meta.total_count` maior que a lista recebida). Elas voltam ao
  chamador como vieram, para que `ListaTruncada`/`RespostaInesperada` sigam
  sendo erros do ITEM, mas uma execucao futura tenta de novo.
* Invalidacao por versao: as entradas ficam em `<diretorio>/<versoes>/`, com a
  versao do Open Targets (`meta.dataVersion`), do ChEMBL e do Reactome no nome
  da pasta. Servico novo -> pasta nova -> nada velho e servido. As consultas
  de versao nunca passam pelo cache.
* Modo de repeticao (`somente_leitura=True`): so le; falta no cache e erro
  (`FaltaNoCache`), nunca cai para a rede.
* Gravacao atomica (temporario + os.replace): execucao interrompida nao deixa
  entrada pela metade.
"""

import hashlib
import json
import os
import re
import tempfile
from pathlib import Path
from typing import Callable

URL_OPEN_TARGETS = "https://api.platform.opentargets.org/api/v4/graphql"
URL_CHEMBL_STATUS = "https://www.ebi.ac.uk/chembl/api/data/status.json"
URL_REACTOME_VERSAO = "https://reactome.org/ContentService/data/database/version"
_BASE_REACTOME = "https://reactome.org/ContentService/"
# openFDA (etapa 2): nao tem numero de versao; `meta.last_updated` identifica
# a carga. Uma busca sem resultado responde 404, que aqui vira lista vazia --
# e uma resposta legitima e cacheavel ("este nome nao tem bula"), nao uma
# falha. Sem isso, a base inteira de produtos naturais, que quase nao tem
# bula, consultaria o openFDA de novo a cada execucao.
URL_OPENFDA_VERSAO = "https://api.fda.gov/drug/label.json?limit=1"
_BASE_OPENFDA = "https://api.fda.gov/"

_CONSULTA_VERSAO_OT = "query Versao { meta { dataVersion { year month iteration } } }"


class FaltaNoCache(Exception):
    """Repeticao pediu uma resposta que o cache nao tem; nao ha rede de reserva.
    Herda de Exception (nao de RuntimeError) de proposito: nenhum tratador de
    item, hoje ou depois de alargado, pode engoli-la; so a `main` a trata."""


class VersaoIndeterminada(Exception):
    """Nao foi possivel saber a versao de um servico (sem ela o cache nao e seguro)."""


def _buscar_padrao(url: str, corpo: dict | None = None):
    """Mesma escolha que cada modulo faz sozinho: POST no Open Targets, Reactome
    com 404 = lista vazia, o resto GET do ChEMBL."""
    from mods.evidencias import chembl, open_targets, reactome
    if corpo is not None:
        return open_targets._post_json_requests(url, corpo)
    if url.startswith(_BASE_REACTOME):
        return reactome._obter_json_requests(url)
    if url.startswith(_BASE_OPENFDA):
        return _obter_json_openfda(url)
    return chembl._obter_json_requests(url)


def _obter_json_openfda(url: str) -> dict:
    """GET no openFDA; 404 e "nenhuma bula com esse nome", nao erro."""
    import requests
    resposta = requests.get(url, timeout=60)
    if resposta.status_code == 404:
        return {"results": [], "openfda_sem_resultado": True}
    resposta.raise_for_status()
    return resposta.json()


def _chave(url: str, corpo: dict | None) -> str:
    canonico = json.dumps({"url": url, "corpo": corpo}, sort_keys=True,
                          ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(canonico.encode("utf-8")).hexdigest()


def e_resposta_confiavel(resposta) -> bool:
    """Falso para o que NAO pode ser tratado como resposta: vazio, com sinal de
    erro, ou pagina truncada."""
    if resposta is None:
        return False
    if not isinstance(resposta, dict):
        return True
    if resposta.get("errors") or resposta.get("error_message") or resposta.get("error"):
        return False
    total = (resposta.get("page_meta") or {}).get("total_count")
    if isinstance(total, int):
        listas = [v for k, v in resposta.items() if k != "page_meta" and isinstance(v, list)]
        if not listas or max(len(v) for v in listas) < total:
            return False
    return True


class Cache:
    def __init__(self, diretorio, buscar: Callable | None = None, *,
                 somente_leitura: bool = False):
        self.diretorio = Path(diretorio)
        self._buscar = buscar or _buscar_padrao
        self.somente_leitura = somente_leitura
        self.acertos = 0
        self.buscas = 0

    def _caminho(self, url: str, corpo: dict | None) -> Path:
        return self.diretorio / (_chave(url, corpo) + ".json")

    def obter_json(self, url: str, corpo: dict | None = None):
        """Serve do disco quando ha; senao busca e grava (se a resposta for
        confiavel). Mesma assinatura que `obter_json` dos modulos de fonte."""
        caminho = self._caminho(url, corpo)
        if caminho.is_file():
            try:
                entrada = json.loads(caminho.read_text(encoding="utf-8"))
                if entrada["url"] == url and entrada["corpo"] == corpo:
                    self.acertos += 1
                    return entrada["resposta"]
            except (OSError, ValueError, KeyError, TypeError):
                pass
            raise RuntimeError(f"entrada de cache corrompida: {caminho}; apague-a")
        if self.somente_leitura:
            raise FaltaNoCache(f"repeticao sem rede: resposta ausente do cache "
                               f"({url}{'' if corpo is None else ' + corpo'})")
        resposta = self._buscar(url) if corpo is None else self._buscar(url, corpo)
        self.buscas += 1
        if e_resposta_confiavel(resposta):
            self._gravar(caminho, {"url": url, "corpo": corpo, "resposta": resposta})
        return resposta

    def _gravar(self, caminho: Path, entrada: dict) -> None:
        caminho.parent.mkdir(parents=True, exist_ok=True)
        fd, nome = tempfile.mkstemp(dir=caminho.parent, suffix=".tmp")
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as f:
                json.dump(entrada, f, ensure_ascii=False)
            os.replace(nome, caminho)
        finally:
            Path(nome).unlink(missing_ok=True)


def _segmento(texto: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]", "_", texto)


_SERVICOS_BASE = ("open_targets", "chembl", "reactome")


def pasta_das_versoes(versoes: dict) -> str:
    """Nome da pasta do cache. Os tres servicos da etapa 1 vem primeiro e nessa
    ordem, para que o nome de um cache ja gravado nao mude; servico extra (o
    `openfda` da etapa 2) entra depois, em ordem alfabetica. Pasta diferente
    por conjunto de servicos e o certo: um cache sem a versao do openFDA nao
    pode servir respostas do openFDA."""
    extras = sorted(k for k in versoes if k not in _SERVICOS_BASE)
    chaves = [k for k in _SERVICOS_BASE if k in versoes] + extras
    return "__".join(f"{k}-{_segmento(str(versoes[k]))}" for k in chaves)


def _versao_valida(valor, servico: str) -> str:
    """Escalar com cara de versao (letras, digitos, '.', '_', '-', com ao menos
    um digito). Lista/dict/bool/vazio (ex.: o `[]` do 404 do Reactome) abortam."""
    if isinstance(valor, bool) or not isinstance(valor, (str, int)):
        raise VersaoIndeterminada(f"versao do {servico} nao e escalar: {valor!r}")
    texto = str(valor).strip()
    if not re.fullmatch(r"[A-Za-z0-9._-]+", texto) or not re.search(r"\d", texto):
        raise VersaoIndeterminada(f"versao do {servico} nao parece versao: {valor!r}")
    return texto


def consultar_versoes(buscar: Callable | None = None, *,
                      com_openfda: bool = False) -> dict:
    """Versoes ATUAIS dos servicos, direto da rede (nunca do cache). Falha com
    `VersaoIndeterminada`: sem versao nao ha como garantir que o cache vale.

    `com_openfda` acrescenta a data de carga do openFDA (etapa 2). Fica
    opcional para que a etapa 1, que nao consulta bula, nao pague a
    requisicao nem mude o nome das pastas que ela ja gravou."""
    buscar = buscar or _buscar_padrao
    try:
        ot = buscar(URL_OPEN_TARGETS, {"query": _CONSULTA_VERSAO_OT, "variables": {}})
        v = ot["data"]["meta"]["dataVersion"]
        if isinstance(v, dict):
            versao_ot = f"{int(v['year']) % 100:02d}.{int(v['month']):02d}"
            if v.get("iteration"):
                versao_ot += f".{v['iteration']}"
        else:
            versao_ot = v
        versao_ot = _versao_valida(versao_ot, "Open Targets")
        chembl = _versao_valida(buscar(URL_CHEMBL_STATUS)["chembl_db_version"], "ChEMBL")
        reactome = _versao_valida(buscar(URL_REACTOME_VERSAO), "Reactome")
        versoes = {"open_targets": versao_ot, "chembl": chembl,
                   "reactome": reactome}
        if com_openfda:
            carga = (buscar(URL_OPENFDA_VERSAO)["meta"]["last_updated"])
            versoes["openfda"] = _versao_valida(carga, "openFDA")
    except (KeyError, TypeError, ValueError) as e:
        raise VersaoIndeterminada(f"resposta de versao inesperada: {e!r}") from e
    except OSError as e:
        raise VersaoIndeterminada(f"consulta de versao falhou: {e}") from e
    return versoes


def abrir_ao_vivo(raiz, buscar: Callable | None = None, *,
                  com_openfda: bool = False) -> tuple[Cache, dict]:
    """Cache da execucao normal: descobre as versoes e usa a pasta delas."""
    versoes = consultar_versoes(buscar, com_openfda=com_openfda)
    pasta = Path(raiz) / pasta_das_versoes(versoes)
    pasta.mkdir(parents=True, exist_ok=True)
    (pasta / "versoes.json").write_text(
        json.dumps(versoes, indent=2) + "\n", encoding="utf-8")
    return Cache(pasta, buscar), versoes


def abrir_para_repeticao(raiz, pasta: str) -> tuple[Cache, dict]:
    """Cache so-leitura da pasta que o manifesto aponta; as versoes sao as que
    ela mesma declara."""
    caminho = Path(raiz) / pasta
    try:
        versoes = json.loads((caminho / "versoes.json").read_text(encoding="utf-8"))
    except (OSError, ValueError) as e:
        raise VersaoIndeterminada(f"cache {caminho} sem versoes.json legivel: {e}") from e
    return Cache(caminho, somente_leitura=True), versoes
