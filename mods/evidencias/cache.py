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

_CONSULTA_VERSAO_OT = "query Versao { meta { dataVersion { year month iteration } } }"


class FaltaNoCache(RuntimeError):
    """Repeticao pediu uma resposta que o cache nao tem; nao ha rede de reserva."""


class VersaoIndeterminada(RuntimeError):
    """Nao foi possivel saber a versao de um servico (sem ela o cache nao e seguro)."""


def _buscar_padrao(url: str, corpo: dict | None = None):
    """Mesma escolha que cada modulo faz sozinho: POST no Open Targets, Reactome
    com 404 = lista vazia, o resto GET do ChEMBL."""
    from mods.evidencias import chembl, open_targets, reactome
    if corpo is not None:
        return open_targets._post_json_requests(url, corpo)
    if url.startswith(_BASE_REACTOME):
        return reactome._obter_json_requests(url)
    return chembl._obter_json_requests(url)


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


def pasta_das_versoes(versoes: dict) -> str:
    return "__".join(f"{k}-{_segmento(str(versoes[k]))}"
                     for k in ("open_targets", "chembl", "reactome"))


def consultar_versoes(buscar: Callable | None = None) -> dict:
    """Versoes ATUAIS dos tres servicos, direto da rede (nunca do cache).
    Falha com `VersaoIndeterminada`: sem versao nao ha como garantir que o
    cache vale."""
    buscar = buscar or _buscar_padrao
    try:
        ot = buscar(URL_OPEN_TARGETS, {"query": _CONSULTA_VERSAO_OT, "variables": {}})
        v = ot["data"]["meta"]["dataVersion"]
        if isinstance(v, dict):
            versao_ot = f"{int(v['year']) % 100:02d}.{int(v['month']):02d}"
            if v.get("iteration"):
                versao_ot += f".{v['iteration']}"
        else:
            versao_ot = str(v)
        chembl = str(buscar(URL_CHEMBL_STATUS)["chembl_db_version"])
        reactome = buscar(URL_REACTOME_VERSAO)
        reactome = str(reactome).strip()
    except (KeyError, TypeError, ValueError) as e:
        raise VersaoIndeterminada(f"resposta de versao inesperada: {e!r}") from e
    except OSError as e:
        raise VersaoIndeterminada(f"consulta de versao falhou: {e}") from e
    if not (versao_ot and chembl and reactome):
        raise VersaoIndeterminada("versao vazia")
    return {"open_targets": versao_ot, "chembl": chembl, "reactome": reactome}


def abrir_ao_vivo(raiz, buscar: Callable | None = None) -> tuple[Cache, dict]:
    """Cache da execucao normal: descobre as versoes e usa a pasta delas."""
    versoes = consultar_versoes(buscar)
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
