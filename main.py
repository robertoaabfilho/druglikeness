#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
main.py - orquestra as tres etapas do pipeline.

    etapa 1  COLETAR    o que se sabe           -> nunca decide
    etapa 2  O PORTAO   o que fica marcado      -> quase nunca reprova
    fase 3   RANQUEAR   a ordem usando os alvos -> faixas de evidencia

Cada etapa roda como processo separado e publica a propria tabela com o
proprio manifesto ao lado. Se a fase 3 falhar, a etapa 2 continua publicada e
citavel -- e por isso que o orquestrador chama por subprocesso em vez de rodar
tudo num processo so.

Uso:

    python main.py                       # as tres etapas, da primeira a ultima
    python main.py --retomar             # pula etapa cuja saida ja existe
    python main.py --a-partir-de portao  # comeca numa etapa
    python main.py --somente ranking     # roda so uma
    python main.py --config outro.json   # outra configuracao

DUAS COISAS QUE NAO SAO CONVENIENCIA
------------------------------------
**Nada e implicito.** Sem `--retomar`, tudo roda. Pular etapa cuja saida existe
e economia legitima, mas feita por omissao ela entrega resultado velho como
novo -- e num pipeline cuja execucao real sao ~23 milhoes de pares, ninguem
notaria.

**`--retomar` confere a ENTRADA, nao so a existencia da saida.** Se o
manifesto da saida aponta outro MD5 de ligantes, a entrada mudou desde que
aquela etapa rodou, e aproveita-la misturaria duas execucoes numa so. Nesse
caso o orquestrador REFAZ aquela etapa, em vez de aproveitar a saida -- e nao
aborta a execucao: mandar recomecar ~23 milhoes de pares por algo que ele pode
refazer sozinho seria pior. Manifesto que nao registra o MD5 sai
`indeterminado`, e a etapa roda de novo.

O codigo de saida e o PIOR entre as etapas que rodaram: 0 completa, 3 parcial
(linhas descartadas), 2 erro. Automacao nao le stderr.
"""
import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from mods import caminhos, config, dependencias  # noqa: E402

COMANDO = "python main.py"

# A tabela declarativa. `saida` e o arquivo que PROVA que a etapa rodou -- e o
# que `--retomar` consulta -- e `manifesto` e o que diz com que entrada ela
# rodou. Tabela e manifesto sao publicados juntos, em temporarios.
ETAPAS = [
    {"id": "etapa1", "script": "coletar_evidencias.py",
     "descricao": "Coletar: evidencia medida (ChEMBL) e prevista (PP4CKD) "
                  "sobre cada par ligante/alvo, mais o enriquecimento de "
                  "alvos (Open Targets e Reactome). Nunca decide.",
     "saida": "tsv/etapa1/evidencias.tsv",
     "manifesto": "tsv/etapa1/manifesto.json"},
    {"id": "portao", "script": "aplicar_portao.py",
     "descricao": "O portao: marca, anota e explica. Por padrao so reprovam "
                  "SMILES que o RDKit nao le e alerta forte de bula.",
     "saida": "tsv/portao/portao.tsv",
     "manifesto": "tsv/portao/manifesto.json"},
    {"id": "ranking", "script": "ranquear.py",
     "descricao": "Ranquear pelas faixas de evidencia (1a, 1b, 2a, 2b) e "
                  "mandar quem nao tem evidencia de alvo para nao_ranqueados.",
     "saida": "tsv/ranking/ranking.tsv",
     "manifesto": "tsv/ranking/manifesto.json"},
]

IDS = [e["id"] for e in ETAPAS]


class EtapaFalhou(Exception):
    def __init__(self, etapa_id: str, codigo: int, segundos: float):
        super().__init__(f"etapa {etapa_id} falhou com codigo {codigo}")
        self.etapa_id, self.codigo, self.segundos = etapa_id, codigo, segundos


def duracao(seg: float) -> str:
    m, s = divmod(int(seg), 60)
    h, m = divmod(m, 60)
    return f"{h}h{m:02d}m{s:02d}s" if h else f"{m}m{s:02d}s" if m else f"{s}s"


def _saida_existe(etapa: dict) -> bool:
    return (caminhos.RAIZ / etapa["saida"]).is_file()


def _md5_de(caminho: Path) -> str | None:
    import hashlib
    if not caminho.is_file():
        return None
    h = hashlib.md5()
    with open(caminho, "rb") as f:
        for bloco in iter(lambda: f.read(1 << 20), b""):
            h.update(bloco)
    return h.hexdigest()


def _entrada_confere(etapa: dict, caminho_ligantes: str | None) -> str:
    """A saida daquela etapa foi gerada com a lista de ligantes de agora?

    Devolve "confere", "difere" ou "indeterminado". As tres sao respostas
    diferentes e nao podem dar a mesma mensagem: "o manifesto registra outro
    MD5" e uma acusacao, e "o manifesto nao registra MD5" nao afirma nada --
    a mesma regra que faz `indeterminado` nunca penalizar um composto.

    As etapas gravam `versoes.ligantes_md5` = MD5 do arquivo que receberam em
    `--ligantes`. E com ESSE arquivo que se compara. (O defeito corrigido aqui
    comparava com o MD5 do CSV de ORIGEM, registrado em `procedencia.json`:
    dois arquivos diferentes, logo nunca iguais, logo `--retomar` recusava
    todas as saidas, sempre -- inclusive numa arvore perfeitamente coerente.)
    """
    manifesto = caminhos.RAIZ / etapa["manifesto"]
    if not manifesto.is_file() or not caminho_ligantes:
        return "indeterminado"
    try:
        dados = json.loads(manifesto.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return "indeterminado"
    gravado = (dados.get("versoes") or {}).get("ligantes_md5")
    if not gravado:
        return "indeterminado"
    atual = _md5_de(Path(caminho_ligantes))
    if atual is None:
        return "indeterminado"
    return "confere" if atual == gravado else "difere"


# Entradas que o pesquisador fornece e que o instantaneo NAO cria: ele converte
# a lista de ligantes, mas alvos e sintomas sao escolha de quem pesquisa.
ENTRADAS_OBRIGATORIAS = {
    "etapa1": ("ligantes", "alvos", "sintomas"),
    "portao": (),
    "ranking": (),
}

# Chaves que apontam para o modelo. Ausentes, a etapa 1 roda em MODO
# DEGRADADO: so com evidencia medida, que nao precisa de modelo nenhum.
CHAVES_DO_MODELO = ("checkpoint", "rotulos")


def _precheca(escolhidas: list[dict], cfg: dict) -> list[str]:
    """Confere os caminhos ANTES de rodar qualquer etapa, e devolve os erros.

    Sem isto a etapa 1 materializa os 107.123 ligantes antes de abrir o
    arquivo de alvos, e um "arquivo nao encontrado" chega depois de minutos de
    RDKit. A pre-checagem custa milissegundos.
    """
    erros = []
    for etapa in escolhidas:
        for chave in ENTRADAS_OBRIGATORIAS[etapa["id"]]:
            valor = (cfg.get(etapa["id"]) or {}).get(chave)
            if not valor:
                erros.append(
                    f"{etapa['id']}.{chave} nao esta em json/config.json")
                continue
            if not (caminhos.RAIZ / valor).exists():
                dica = ""
                if chave in ("alvos", "sintomas"):
                    dica = (f" -- este arquivo e SEU: escreva-o, ou copie um "
                            f"de exemplos/ para comecar")
                elif chave == "ligantes":
                    dica = (" -- traga a lista com "
                            "`python mods/instantaneo.py --ligantes-origem ...`")
                erros.append(f"{etapa['id']}.{chave}: {valor} nao existe{dica}")
    erros += _precheca_das_regras(escolhidas, cfg)
    return erros


# Regras ENTRE chaves que uma etapa recusa no proprio comeco. Conferi-las aqui
# custa milissegundos; deixa-las para a etapa significa descobrir o erro depois
# de a etapa 1 materializar 107.123 ligantes -- horas numa execucao real.
REGRAS_ENTRE_CHAVES = {
    "ranking": [(
        "checkpoint", "manifesto_etapa1",
        "sem ele nao ha como conferir que `p_afirmativo` (da etapa 1) e "
        "`p_negativo` (da fase 3) vem do MESMO modelo")],
}


# Chaves opcionais que, quando preenchidas, apontam um arquivo que TEM de
# existir -- e a conferencia vale antes de qualquer etapa, nao quando a etapa
# que o usa comeca. O acervo do openFDA e o caso: a etapa 2 so roda depois da
# etapa 1, que numa execucao real leva horas.
ARQUIVOS_SE_PEDIDOS = {
    "portao": {"acervo_bula": "baixe com `python mods/openfda_massa.py "
                              "--baixar --indexar`, ou rode `python main.py "
                              "--baixar-acervo`"},
}


def _precheca_das_regras(escolhidas: list[dict], cfg: dict) -> list[str]:
    """Chave que EXIGE outra chave, e chave opcional que aponta arquivo
    ausente -- as duas conferidas antes de rodar."""
    erros = []
    for etapa in escolhidas:
        secao = cfg.get(etapa["id"]) or {}
        for chave, remedio in ARQUIVOS_SE_PEDIDOS.get(etapa["id"], {}).items():
            valor = secao.get(chave)
            if valor and not (caminhos.RAIZ / valor).exists():
                erros.append(
                    f"{etapa['id']}.{chave}: {valor} nao existe -- {remedio}")
    for etapa in escolhidas:
        for exige, exigida, motivo in REGRAS_ENTRE_CHAVES.get(etapa["id"], ()):
            secao = cfg.get(etapa["id"]) or {}
            if secao.get(exige) and not secao.get(exigida):
                erros.append(
                    f"{etapa['id']}.{exige} esta na configuracao mas "
                    f"{etapa['id']}.{exigida} nao: {motivo}")
    return erros


def _baixar_acervo(cfg: dict) -> bool:
    """Baixa e indexa o acervo do openFDA, quando `--baixar-acervo` foi pedido.

    So por pedido explicito. Sao ~3,7 GB e ~20 min, e existe alternativa (a
    rede, com o teto de taxa) -- a mesma regra do `torch`: baixar GB sem pedir
    e decidir pelo pesquisador num ponto em que ele tem escolha. E, para uma
    base de produto natural, o acervo compra pouco: medido, 58 de 60 ligantes
    nao tem bula nenhuma no openFDA.

    O destino e o que a configuracao aponta em `portao.acervo_bula`; sem ele,
    o padrao do modulo."""
    from mods import openfda_massa
    indice = (cfg.get("portao") or {}).get("acervo_bula")
    indice = (caminhos.RAIZ / indice) if indice else None
    pasta = indice.parent if indice else caminhos.RAIZ / "dados/openfda"
    if indice and indice.is_file():
        print(f"[info] acervo do openFDA ja esta em {indice}; nada a baixar",
              file=sys.stderr)
        return True
    print(f"[info] baixando o acervo do openFDA em {pasta} "
          f"(~1,7 GB em 14 particoes, mais o indice; retoma se interromper)",
          file=sys.stderr)
    codigo = openfda_massa.main(
        ["--baixar", "--indexar", "--saida", str(pasta)]
        + (["--indice", str(indice)] if indice else []))
    if codigo != 0:
        print("[ERRO] o acervo nao foi construido; nada foi executado",
              file=sys.stderr)
        return False
    return True


def _degradar_sem_modelo(cfg: dict) -> list[str]:
    """Tira as chaves do modelo que apontam para arquivo ausente.

    `modelo/` esta no .gitignore, logo um clone nao tem pesos -- e a evidencia
    MEDIDA (ChEMBL) nao precisa deles. Passar um checkpoint inexistente fazia
    a etapa 1 sair com codigo 2 e derrubar o pipeline inteiro, em vez de
    coletar o que da e dizer o que faltou.
    """
    avisos = []
    valores = cfg.get("etapa1") or {}
    faltando = [c for c in CHAVES_DO_MODELO
                if valores.get(c) and not (caminhos.RAIZ / valores[c]).exists()]
    if not faltando:
        return avisos
    for chave in CHAVES_DO_MODELO:
        valores.pop(chave, None)
    avisos.append(
        "[AVISO] arquivo do modelo ausente (" + ", ".join(faltando) + "): a "
        "etapa 1 roda em MODO DEGRADADO, so com evidencia MEDIDA; nao havera "
        "evidencia prevista, e a fase 3 nao tera canal negativo (a coluna sai "
        "`indeterminado`, nunca 0). Para trazer os pesos: "
        "`python mods/instantaneo.py --modelo-origem PASTA_DOS_PESOS`.")
    return avisos


def _executar(etapa: dict, cfg: dict) -> tuple[int, float]:
    """Roda a etapa como subprocesso. Devolve (codigo, segundos)."""
    script = caminhos.MODS / etapa["script"]
    comando = [sys.executable, str(script)] + config.argumentos_de(etapa["id"], cfg)
    ambiente = os.environ.copy()
    ambiente.setdefault("PYTHONUNBUFFERED", "1")
    ambiente.setdefault("PYTHONUTF8", "1")
    inicio = time.time()
    r = subprocess.run(comando, cwd=str(caminhos.RAIZ), env=ambiente)
    return r.returncode, time.time() - inicio


def rodar(etapa: dict, cfg: dict) -> tuple[int, float]:
    """Anuncia, roda e devolve (codigo, segundos). Levanta EtapaFalhou em erro.

    A duracao volta porque numa execucao de horas o resumo e o unico registro
    de onde o tempo foi.
    """
    print(f"\n{'=' * 70}\nETAPA [{etapa['id']}]: {etapa['descricao']}\n"
          f"Script: mods/{etapa['script']}\n{'=' * 70}\n", flush=True)
    codigo, segundos = _executar(etapa, cfg)
    if codigo not in (0, 3):
        raise EtapaFalhou(etapa["id"], codigo, segundos)
    if codigo == 3:
        print(f"\nEtapa '{etapa['id']}' PARCIAL em {duracao(segundos)}: "
              f"linhas de entrada foram descartadas. O manifesto registra "
              f"`execucao_parcial`.", flush=True)
    else:
        print(f"\nEtapa '{etapa['id']}' concluida em {duracao(segundos)}.",
              flush=True)
    return codigo, segundos


def _resumo(tempos: dict, pior: int) -> None:
    if not tempos:
        return
    print(f"\n{'=' * 70}\nRESUMO")
    for etapa_id, seg in tempos.items():
        print(f"  {etapa_id:<10} {duracao(seg)}")
    print(f"  {'total':<10} {duracao(sum(tempos.values()))}")
    veredito = {0: "completa", 3: "PARCIAL (linhas descartadas)"}.get(
        pior, f"codigo {pior}")
    print(f"Execucao {veredito}.\n{'=' * 70}")


def _argumentos(argv):
    p = argparse.ArgumentParser(
        description="Orquestra as tres etapas do pipeline.")
    p.add_argument("--config", default=None,
                   help="Configuracao (padrao: json/config.json).")
    p.add_argument("--retomar", action="store_true",
                   help="Pula etapa cuja saida ja existe E foi gerada com a "
                        "entrada atual. Sem esta flag, tudo roda.")
    p.add_argument("--a-partir-de", dest="a_partir_de", choices=IDS,
                   help="Comeca nesta etapa.")
    p.add_argument("--somente", choices=IDS, help="Roda so esta etapa.")
    p.add_argument("--baixar-acervo", dest="baixar_acervo",
                   action="store_true",
                   help="Baixa e indexa o acervo do openFDA antes de rodar "
                        "(~3,7 GB com o indice; retoma). NAO e automatico: ha "
                        "alternativa (a rede, com o teto de taxa), e para uma "
                        "base de produto natural o acervo compra pouco -- a "
                        "maioria dos compostos nao tem bula. Pedir e seu.")
    p.add_argument("--sem-instalar", dest="sem_instalar", action="store_true",
                   help="Confere as dependencias e NAO instala nada. Para "
                        "ambiente compartilhado, automacao ou maquina sem "
                        "rede: o que falta e dito, e a execucao para.")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = _argumentos(argv)

    # ANTES de qualquer coisa: o ambiente. Uma execucao real sao ~23 milhoes
    # de pares, e um ModuleNotFoundError no meio dela custa horas por algo que
    # custa milissegundos conferir -- a mesma razao de `_precheca`.
    ambiente = dependencias.garantir(instalar=not args.sem_instalar)
    for linha in dependencias.relato(ambiente):
        print(linha, file=sys.stderr)
    if not ambiente["completo"]:
        print("[ERRO] o ambiente nao esta completo; nada foi executado. "
              "Nenhuma etapa roda com biblioteca faltando -- ela quebraria no "
              "meio, depois de horas.", file=sys.stderr)
        return 2

    try:
        cfg = config.ler(args.config)
    except config.ErroDeConfig as e:
        print(f"[ERRO] configuracao: {e}", file=sys.stderr)
        return 2

    if args.baixar_acervo and not _baixar_acervo(cfg):
        return 2

    if args.somente:
        escolhidas = [e for e in ETAPAS if e["id"] == args.somente]
    elif args.a_partir_de:
        escolhidas = ETAPAS[IDS.index(args.a_partir_de):]
    else:
        escolhidas = list(ETAPAS)

    for aviso in _degradar_sem_modelo(cfg):
        print(aviso, file=sys.stderr)
    erros = _precheca(escolhidas, cfg)
    if erros:
        print("[ERRO] a configuracao aponta para entradas que nao existem:",
              file=sys.stderr)
        for erro in erros:
            print(f"  - {erro}", file=sys.stderr)
        return 2

    # A conferencia de `--retomar` compara o arquivo que as etapas recebem em
    # `--ligantes`, e e a configuracao que diz qual e.
    ligantes_atual = (cfg.get("etapa1") or {}).get("ligantes")
    if ligantes_atual:
        ligantes_atual = str(caminhos.RAIZ / ligantes_atual)

    tempos, pior = {}, 0
    for etapa in escolhidas:
        if args.retomar and _saida_existe(etapa):
            veredito = _entrada_confere(etapa, ligantes_atual)
            if veredito == "confere":
                print(f"[Pulando '{etapa['id']}': {etapa['saida']} ja existe e "
                      f"confere com a lista de ligantes atual (--retomar)]")
                continue
            motivo = {
                "difere": "foi gerada com OUTRA lista de ligantes",
                "indeterminado": "nao da para saber com que entrada foi "
                                 "gerada (o manifesto nao registra o MD5 dos "
                                 "ligantes)",
            }[veredito]
            print(f"[--retomar NAO pula '{etapa['id']}': {etapa['saida']} "
                  f"existe, mas {motivo}. Refazendo esta etapa -- aproveitar "
                  f"aquela saida misturaria duas execucoes.]")
        try:
            codigo, segundos = rodar(etapa, cfg)
        except EtapaFalhou as e:
            tempos[e.etapa_id] = e.segundos
            print(f"\n[ERRO] a etapa '{e.etapa_id}' falhou (codigo {e.codigo}) "
                  f"apos {duracao(e.segundos)}.\nCorrija e retome com: "
                  f"{COMANDO} --a-partir-de {e.etapa_id}", file=sys.stderr)
            _resumo(tempos, e.codigo)
            return e.codigo
        tempos[etapa["id"]] = segundos
        pior = max(pior, codigo)
    _resumo(tempos, pior)
    return pior


if __name__ == "__main__":
    sys.exit(main())
