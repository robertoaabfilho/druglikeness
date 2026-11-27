#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
openfda_massa.py - o acervo do openFDA baixado e respondido localmente.

Por que existe
-------------
`alerta_bula` e um dos DOIS unicos criterios que REPROVAM no portao. Com o teto
de taxa do servico publico (`--limite-bula`, padrao 500) sobre 108.935
ligantes, ele sai `indeterminado` para **99,5% da base**: o criterio existe e
quase nao age. O teto nunca foi uma decisao de desenho -- era a taxa da API.

Com o acervo local nao ha taxa: as 265 mil bulas ficam num indice SQLite e a
consulta e local. O teto pode ser removido, e o criterio passa a valer para a
base inteira.

A ideia que faz isto caber em pouco codigo
------------------------------------------
Este modulo **responde as mesmas URLs** que `mods/portao/bula.py` monta. O
`bula.py` nao muda uma linha: ele recebe outro `obter_json`, e toda a logica
medida continua valendo igual -- a janela de 50 bulas por data decrescente, a
consulta dirigida sobre o acervo inteiro, e o desconto da indicacao, que e o
que leva o controle de 12 a 0 de 6 falsos positivos no grupo que trata.

Tres URLs, e so elas
--------------------
1. `?limit=1`                          -> a versao do acervo
2. `?search=openfda.<campo>:"<nome>"&sort=...&limit=N`   -> a janela
3. `?search=openfda.<campo>:"<nome>" AND (secao:"termo" OR ...)&limit=1`
                                       -> a dirigida, com `meta.results.total`

Qualquer outra forma levanta erro que NOMEIA a URL. Devolver vazio em silencio
faria uma consulta nova sair "sem alerta" sem ninguem notar -- e ausencia de
documento nunca e ausencia de risco.

O que o indice guarda
---------------------
Somente os campos que o criterio le: os tres nomes do `openfda`, a data de
vigencia e as secoes de risco, de reacao adversa e de indicacao. A bula tem 34
campos e o acervo bruto 1,7 GB; guardar tudo multiplicaria o tamanho sem mudar
resposta nenhuma.

Uso
---
    python mods/openfda_massa.py --baixar --saida dados/openfda
    python mods/openfda_massa.py --indexar --acervo dados/openfda \\
        --indice dados/openfda/openfda.sqlite

    # e depois, na etapa 2, sem teto de taxa:
    python mods/aplicar_portao.py ... --acervo-bula dados/openfda/openfda.sqlite

O download sao ~1,7 GB em 14 particoes, e **retoma**: particao ja baixada e
pulada, particao truncada e refeita. Interromper e seguro.
"""
import argparse
import gzip
import json
import os
import re
import sqlite3
import sys
import time
import zipfile
from pathlib import Path
from urllib.parse import parse_qs, unquote, urlparse

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from mods.evidencias import atomico  # noqa: E402

CATALOGO = "https://api.fda.gov/download.json"

# As secoes que o criterio le, com o mesmo significado que `bula.py` lhes da.
# A lista vem de `bula.SECOES_*`, nao de literais repetidos aqui: se o criterio
# passar a ler outra secao, o indice tem de guarda-la.
def _secoes_do_criterio() -> tuple[str, ...]:
    from mods.portao import bula
    return tuple(dict.fromkeys(
        bula.SECOES_FORTES + bula.SECOES_FRACAS + bula.SECOES_INDICACAO))


# Campos de nome do `openfda`, na ordem em que `bula.buscar_bulas` os tenta.
def _campos_de_nome() -> tuple[str, ...]:
    from mods.portao import bula
    return tuple(bula.CAMPOS_NOME)


SAIDA_OK, SAIDA_ERRO = 0, 2


# =========================================================================
# o catalogo e o download
# =========================================================================
def _label(catalogo: dict) -> dict:
    alvo = ((catalogo or {}).get("results") or {}).get("drug", {}).get("label")
    if not alvo:
        raise ValueError(
            "download.json sem results.drug/label: o catalogo do openFDA mudou "
            "de forma, e baixar as particoes erradas daria um acervo que "
            "responde 'sem alerta' por omissao")
    return alvo


def data_do_export(catalogo: dict) -> str:
    """A data da carga -- a versao do acervo.

    O openFDA nao publica numero de versao; e a mesma regra de
    `bula.versao_openfda`, que usa `meta.last_updated`."""
    return str(_label(catalogo).get("export_date") or "")


def _nome_de_arquivo(url: str) -> str:
    """O nome do arquivo no disco, a partir da URL, higienizado.

    **O nome NUNCA vem do `display_name`.** O openFDA publica
    `drug/label (part 1 of 14)` ali -- com barra --, e usa-lo como nome de
    arquivo fazia o download escrever em `/drug/label (part 1 of 14).json.zip`,
    fora da pasta de destino. Um nome vindo de servico remoto nao decide
    caminho no disco.

    E nem a URL decide sozinha: fica so a ultima parte do caminho, e qualquer
    separador ou `..` que sobre e trocado por `_`."""
    bruto = Path(urlparse(url).path).name or "particao.json.zip"
    limpo = re.sub(r"[^A-Za-z0-9._-]", "_", bruto).lstrip(".")
    return limpo or "particao.json.zip"


def particoes(catalogo: dict) -> list[dict]:
    """As particoes do export.

    `arquivo` e o nome no disco (da URL, higienizado) e `nome` e o rotulo
    legivel que o servico publica -- este ultimo serve para a mensagem, nunca
    para o caminho."""
    saida = []
    for parte in _label(catalogo).get("partitions") or []:
        url = str(parte.get("file") or "")
        if not url:
            continue
        try:
            mb = float(parte.get("size_mb") or 0)
        except (TypeError, ValueError):
            mb = 0.0
        arquivo = _nome_de_arquivo(url)
        saida.append({
            "arquivo": arquivo,
            "nome": str(parte.get("display_name") or "").strip() or arquivo,
            "url": url, "mb": mb,
            "registros": int(parte.get("records") or 0)})
    if not saida:
        raise ValueError("download.json sem particoes em drug/label")
    return saida


def _zip_valido(caminho: Path) -> bool:
    """Um zip que nao abre e download interrompido, nao acervo.

    Aceita-lo na retomada daria acervo INCOMPLETO apresentado como completo --
    e `indeterminado` virando "sem alerta" para o ligante cuja bula estava na
    particao truncada."""
    try:
        with zipfile.ZipFile(caminho) as z:
            return z.testzip() is None and bool(z.namelist())
    except (zipfile.BadZipFile, OSError):
        return False


def baixar_url(url: str, destino: Path) -> None:
    """Baixa em fluxo, sem carregar o arquivo na memoria."""
    import urllib.request
    pedido = urllib.request.Request(
        url, headers={"User-Agent": "druglikeness-pipeline"})
    with urllib.request.urlopen(pedido, timeout=120) as resposta, \
            open(destino, "wb") as saida:
        while True:
            bloco = resposta.read(1 << 20)
            if not bloco:
                break
            saida.write(bloco)


def baixar(catalogo: dict, pasta, *, baixador=baixar_url,
           so_faltantes: bool = True) -> dict:
    """Baixa as particoes para `pasta` e devolve o relatorio.

    **Retoma.** Particao presente e valida e pulada; particao truncada e
    refeita. Interromper um download de 1,7 GB e seguro.

    **Escrita atomica.** Cada particao nasce num temporario ao lado do destino
    e so toma o lugar dele quando termina inteira -- uma falha no meio nao
    deixa um zip parcial que a retomada aceitaria depois."""
    pasta = Path(pasta)
    pasta.mkdir(parents=True, exist_ok=True)
    partes = particoes(catalogo)
    rel = {"baixadas": 0, "ja_tinha": 0, "mb": 0.0, "particoes": len(partes)}
    for parte in partes:
        destino = pasta / parte["arquivo"]
        if so_faltantes and destino.is_file() and _zip_valido(destino):
            rel["ja_tinha"] += 1
            continue
        tmp = atomico.temporario(destino)
        try:
            baixador(parte["url"], tmp)
            if not _zip_valido(tmp):
                raise OSError(
                    f"{parte['nome']}: o arquivo baixado nao e um zip valido")
            atomico.publicar(tmp, destino)
        finally:
            tmp.unlink(missing_ok=True)
        rel["baixadas"] += 1
        rel["mb"] += parte["mb"]
    atomico.escrever_texto(pasta / "procedencia.json", json.dumps({
        "fonte": CATALOGO,
        "export_date": data_do_export(catalogo),
        "total_records": int(_label(catalogo).get("total_records") or 0),
        "particoes": len(partes),
        "baixado_em": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }, ensure_ascii=False, indent=2) + "\n")
    return rel


# =========================================================================
# o indice
# =========================================================================
def _rotulos_do_arquivo(caminho: Path):
    """Gera as bulas de uma particao (.json.zip, .json.gz ou .json).

    Aceitar os tres evita exigir que quem baixou a mao recompacte 1,7 GB."""
    nome = caminho.name.lower()
    if nome.endswith(".zip"):
        with zipfile.ZipFile(caminho) as z:
            for interno in z.namelist():
                if not interno.lower().endswith(".json"):
                    continue
                with z.open(interno) as arq:
                    dados = json.load(arq)
                yield from dados.get("results") or []
        return
    abrir = gzip.open if nome.endswith(".gz") else open
    with abrir(caminho, "rt", encoding="utf-8") as arq:
        dados = json.load(arq)
    yield from dados.get("results") or []


def _texto(valor) -> str:
    """As secoes da bula vem como lista de strings; o criterio le texto."""
    if isinstance(valor, list):
        return "\n".join(str(v) for v in valor)
    return "" if valor is None else str(valor)


def _nomes(rotulo: dict, campo: str) -> list[str]:
    valor = (rotulo.get("openfda") or {}).get(campo) or []
    if not isinstance(valor, list):
        valor = [valor]
    return [str(v).strip() for v in valor if str(v).strip()]


def indexar(pasta, indice, *, export_date: str) -> dict:
    """Constroi o indice SQLite a partir das particoes em `pasta`.

    Uma linha por bula, com os nomes numa tabela propria (uma bula tem varios
    genericos e varias marcas), e so as secoes que o criterio le."""
    pasta, indice = Path(pasta), Path(indice)
    arquivos = sorted(
        p for p in pasta.iterdir()
        if p.is_file() and p.name.lower().endswith(
            (".json.zip", ".json.gz", ".json")) and p.name != "procedencia.json")
    if not arquivos:
        raise ValueError(
            f"nenhuma particao em {pasta}: baixe com "
            "`python mods/openfda_massa.py --baixar`")
    secoes = _secoes_do_criterio()
    indice.parent.mkdir(parents=True, exist_ok=True)
    tmp = atomico.temporario(indice)
    rel = {"bulas": 0, "nomes": 0, "arquivos": len(arquivos)}
    try:
        tmp.unlink(missing_ok=True)  # sqlite cria o proprio arquivo
        con = sqlite3.connect(tmp)
        try:
            con.execute("PRAGMA journal_mode=OFF")
            con.execute("PRAGMA synchronous=OFF")
            colunas = ", ".join(f'"{s}" TEXT' for s in secoes)
            con.execute(f"CREATE TABLE bulas (id INTEGER PRIMARY KEY, "
                        f"effective_time TEXT, {colunas})")
            con.execute("CREATE TABLE nomes (bula INTEGER, campo TEXT, "
                        "nome TEXT)")
            con.execute("CREATE TABLE meta (chave TEXT PRIMARY KEY, "
                        "valor TEXT)")
            for arquivo in arquivos:
                for rotulo in _rotulos_do_arquivo(arquivo):
                    rel["bulas"] += 1
                    ident = rel["bulas"]
                    con.execute(
                        f"INSERT INTO bulas (id, effective_time, "
                        f"{', '.join(chr(34) + s + chr(34) for s in secoes)}) "
                        f"VALUES ({', '.join('?' * (len(secoes) + 2))})",
                        [ident, str(rotulo.get("effective_time") or "")]
                        + [_texto(rotulo.get(s)) for s in secoes])
                    for campo in _campos_de_nome():
                        for nome in _nomes(rotulo, campo):
                            con.execute(
                                "INSERT INTO nomes VALUES (?, ?, ?)",
                                (ident, campo, nome.casefold()))
                            rel["nomes"] += 1
            con.execute("CREATE INDEX i_nomes ON nomes (nome, campo)")
            con.executemany("INSERT INTO meta VALUES (?, ?)", [
                ("export_date", str(export_date)),
                ("bulas", str(rel["bulas"])),
                ("construido_em", time.strftime("%Y-%m-%dT%H:%M:%S%z")),
                ("secoes", ",".join(secoes))])
            con.commit()
        finally:
            con.close()
        atomico.publicar(tmp, indice)
    finally:
        tmp.unlink(missing_ok=True)
    return rel


# =========================================================================
# responder as URLs do bula.py
# =========================================================================
_RE_NOME = re.compile(r'openfda\.(\w+):"([^"]*)"')
_RE_TERMO = re.compile(r'(\w+):"([^"]*)"')


class Acervo:
    """Responde as URLs do `bula.py` a partir do indice local.

    Entra no lugar de `cache.obter_json`, e por isso o criterio inteiro roda
    sem rede e sem teto de taxa."""

    def __init__(self, indice):
        self.indice = Path(indice)
        if not self.indice.is_file():
            raise FileNotFoundError(
                f"indice do openFDA nao encontrado: {self.indice}. "
                "Construa com `python mods/openfda_massa.py --indexar`")
        # `Path.as_uri()` e nao concatenacao: um caminho do Windows
        # (`C:\\x\\y.sqlite`) nao e URI valida -- a barra invertida e os
        # dois-pontos da unidade quebram. `as_uri()` da `file:///C:/x/y`.
        self._con = sqlite3.connect(
            f"{self.indice.resolve().as_uri()}?mode=ro", uri=True)
        self._con.row_factory = sqlite3.Row
        self.secoes = tuple(
            (self._meta("secoes") or "").split(",")) if self._meta("secoes") \
            else _secoes_do_criterio()

    def _meta(self, chave: str):
        linha = self._con.execute(
            "SELECT valor FROM meta WHERE chave=?", (chave,)).fetchone()
        return linha[0] if linha else None

    @property
    def versao(self) -> str:
        return f"openfda-massa {self._meta('export_date')}"

    # ---------------------------------------------------------------
    def obter_json(self, url: str) -> dict:
        """A resposta que o openFDA daria para aquela URL."""
        partes = urlparse(url)
        # O endpoint tem de ser o das BULAS. `?limit=1` de `drug/event` ou de
        # qualquer outro endpoint nao e a sonda de versao daqui, e responde-lo
        # como se fosse faria uma consulta nova sair "sem alerta" calada.
        if not partes.path.rstrip("/").endswith("/drug/label.json"):
            raise ValueError(
                f"openfda_massa: nao sei responder esta URL (o acervo e de "
                f"drug/label): {url}")
        consulta = parse_qs(partes.query)
        busca = unquote(consulta.get("search", [""])[0])
        limite = int(consulta.get("limit", ["1"])[0])
        if not busca:
            # `?limit=1`: a sonda de versao.
            return {"meta": {"last_updated": self.versao},
                    "results": []}
        achado_nome = _RE_NOME.search(busca)
        if not achado_nome:
            raise ValueError(
                f"openfda_massa: nao sei responder esta URL: {url}")
        campo, nome = achado_nome.group(1), achado_nome.group(2)
        resto = busca[achado_nome.end():]
        if " AND (" in resto:
            return self._dirigida(campo, nome, resto)
        if resto.strip() not in ("", ")"):
            raise ValueError(
                f"openfda_massa: nao sei responder esta URL: {url}")
        return self._janela(campo, nome, limite)

    # ---------------------------------------------------------------
    def _ids(self, campo: str, nome: str) -> list[int]:
        """As bulas cujo campo de nome contem `nome` como PALAVRA INTEIRA.

        MEDIDO contra a API: `openfda.generic_name:"ciprofloxacin"` responde
        **166** bulas. Casar por igualdade de string dava 76 -- a diferenca
        sao os nomes compostos (`ciprofloxacin hydrochloride`, `ciprofloxacin
        and dexamethasone`), que a API tokeniza e casa. Importa porque esse
        total e o **peso da evidencia** que vai para `n_bulas_com_o_termo`:
        63 de 76 e outra afirmacao que 118 de 166.

        Palavra inteira, nao trecho: `ciprofloxacin` nao casa `ciprofloxacina`
        nem `norciprofloxacin`. E a mesma regra da busca por termo nas secoes.

        O filtro grosso e feito no SQLite (`LIKE`, que usa o indice de nome) e
        a fronteira de palavra em Python: um `LIKE '%x%'` sozinho casaria
        dentro de outra palavra."""
        alvo = " ".join(nome.strip().casefold().split())
        if not alvo:
            return []
        padrao = re.compile(r"\b" + re.escape(alvo) + r"\b")
        vistos: list[int] = []
        visto = set()
        for ident, texto in self._con.execute(
                "SELECT bula, nome FROM nomes WHERE campo=? AND nome LIKE ?",
                (campo, f"%{alvo}%")):
            if ident not in visto and padrao.search(texto):
                visto.add(ident)
                vistos.append(ident)
        return vistos

    def _bulas(self, ids: list[int], limite: int | None) -> list[dict]:
        if not ids:
            return []
        marcas = ",".join("?" * len(ids))
        sql = (f"SELECT * FROM bulas WHERE id IN ({marcas}) "
               "ORDER BY effective_time DESC")
        if limite is not None:
            sql += f" LIMIT {int(limite)}"
        return [self._como_bula(r) for r in self._con.execute(sql, ids)]

    def _como_bula(self, linha: sqlite3.Row) -> dict:
        """A bula no formato do openFDA: secao como LISTA de strings, que e o
        que `bula.achados` espera."""
        bula = {"effective_time": linha["effective_time"]}
        for secao in self.secoes:
            texto = linha[secao]
            if texto:
                bula[secao] = [texto]
        nomes = {}
        for r in self._con.execute(
                "SELECT campo, nome FROM nomes WHERE bula=?", (linha["id"],)):
            nomes.setdefault(r[0], []).append(r[1])
        bula["openfda"] = nomes
        return bula

    def _janela(self, campo: str, nome: str, limite: int) -> dict:
        ids = self._ids(campo, nome)
        return {"meta": {"last_updated": self.versao,
                         "results": {"total": len(ids)}},
                "results": self._bulas(ids, limite)}

    def _dirigida(self, campo: str, nome: str, resto: str) -> dict:
        """"Existe alguma bula deste principio ativo com algum destes termos em
        alguma destas secoes?" -- sobre o acervo INTEIRO, nao sobre a janela.

        `meta.results.total` e quantas bulas casam: o peso da evidencia (118 de
        166 no ciprofloxacino contra 1 de 211 na indicacao da capsaicina), que
        uma contagem feita dentro da janela nao sabe dizer.

        O casamento e por PALAVRA INTEIRA, a mesma regra do criterio na rede:
        sem isso "neuropath" casaria dentro de "neuropathy" e o falso positivo
        nao seria nem conferivel."""
        pedidos: dict[str, list[str]] = {}
        for secao, termo in _RE_TERMO.findall(resto):
            if secao in self.secoes and termo.strip():
                pedidos.setdefault(secao, []).append(termo.strip())
        ids = self._ids(campo, nome)
        if not ids or not pedidos:
            return {"meta": {"last_updated": self.versao,
                             "results": {"total": 0}}, "results": []}
        padroes = {s: [re.compile(r"\b" + re.escape(t) + r"\b", re.IGNORECASE)
                       for t in ts] for s, ts in pedidos.items()}
        casaram = []
        for bula in self._bulas(ids, None):
            if any(p.search(_texto(bula.get(secao)))
                   for secao, ps in padroes.items() for p in ps):
                casaram.append(bula)
        return {"meta": {"last_updated": self.versao,
                         "results": {"total": len(casaram)}},
                "results": casaram[:1]}


# =========================================================================
# CLI
# =========================================================================
def _argumentos(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--baixar", action="store_true",
                    help="baixa as particoes do openFDA (~1,7 GB; retoma)")
    ap.add_argument("--indexar", action="store_true",
                    help="constroi o indice SQLite a partir das particoes")
    ap.add_argument("--saida", metavar="PASTA", default="dados/openfda",
                    help="pasta das particoes (padrao dados/openfda)")
    ap.add_argument("--acervo", metavar="PASTA",
                    help="pasta das particoes para indexar (padrao: --saida)")
    ap.add_argument("--indice", metavar="ARQUIVO",
                    help="indice a escrever (padrao <acervo>/openfda.sqlite)")
    ap.add_argument("--refazer", action="store_true",
                    help="baixa tudo de novo, mesmo o que ja esta no disco")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    args = _argumentos(argv)
    if not (args.baixar or args.indexar):
        print("[erro] escolha --baixar, --indexar, ou os dois", file=sys.stderr)
        return SAIDA_ERRO
    pasta = Path(args.acervo or args.saida)
    indice = Path(args.indice) if args.indice else pasta / "openfda.sqlite"
    catalogo = None
    try:
        if args.baixar:
            import urllib.request
            print(f"[info] catalogo: {CATALOGO}", file=sys.stderr)
            with urllib.request.urlopen(CATALOGO, timeout=60) as r:
                catalogo = json.loads(r.read().decode("utf-8"))
            print(f"[info] export do openFDA: {data_do_export(catalogo)}; "
                  f"{len(particoes(catalogo))} particoes", file=sys.stderr)
            rel = baixar(catalogo, pasta, so_faltantes=not args.refazer)
            print(f"[info] baixadas {rel['baixadas']}, ja tinha "
                  f"{rel['ja_tinha']} ({rel['mb']:.0f} MB novos)",
                  file=sys.stderr)
        if args.indexar:
            if catalogo is None:
                proc = pasta / "procedencia.json"
                export = (json.loads(proc.read_text(encoding="utf-8"))
                          .get("export_date", "desconhecido")
                          if proc.is_file() else "desconhecido")
            else:
                export = data_do_export(catalogo)
            rel = indexar(pasta, indice, export_date=export)
            tamanho = indice.stat().st_size / (1 << 20)
            print(f"[info] indice em {indice}: {rel['bulas']:,} bulas, "
                  f"{rel['nomes']:,} nomes, {tamanho:.0f} MB", file=sys.stderr)
            print(f"[info] use com: python mods/aplicar_portao.py ... "
                  f"--acervo-bula {indice}", file=sys.stderr)
    except (OSError, ValueError, json.JSONDecodeError) as erro:
        print(f"[erro] {erro}", file=sys.stderr)
        return SAIDA_ERRO
    return SAIDA_OK


if __name__ == "__main__":
    sys.exit(main())
