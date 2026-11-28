# -*- coding: utf-8 -*-
"""
atomico.py

Escrita atomica das saidas: o resultado anterior do usuario so e substituido
quando a escrita nova termina inteira.

Por que
-------
As tres saidas da etapa sao arquivos que o pesquisador guarda e cita:
`evidencias.tsv`, o `alvos_enriquecidos.tsv` que o README manda editar a mao, e
o `manifesto.json`, que e a chave da repeticao. Abrir o destino em modo "w"
destroi o anterior antes de a execucao nova saber se vai dar certo -- uma falha
no meio nao deixa o resultado velho nem o novo.

Como
----
O arquivo nasce num temporario NO MESMO diretorio do destino (mesmo sistema de
arquivos, para que publicar seja uma troca de nome, nao uma copia) e so toma o
lugar do destino ao fim. Preparar e publicar ficam separados porque a etapa
publica tres arquivos juntos: todos sao preparados antes de qualquer um ser
publicado, para que uma falha nao deixe tabela sem manifesto.
"""

import os
import tempfile
from pathlib import Path


def temporario(destino) -> Path:
    """Temporario vazio ao lado de `destino`, pronto para ser escrito."""
    destino = Path(destino)
    fd, nome = tempfile.mkstemp(dir=destino.parent, prefix=destino.name + ".",
                                suffix=".tmp")
    os.close(fd)
    return Path(nome)


def publicar(tmp, destino) -> None:
    """Troca o temporario pelo destino.

    `mkstemp` cria 0600; o arquivo publicado recebe o modo de uma escrita
    normal (0666 menos a umask), igual ao das outras saidas."""
    mascara = os.umask(0)
    os.umask(mascara)
    try:
        os.chmod(tmp, 0o666 & ~mascara)
    except OSError:
        pass  # ajuste de cortesia: nunca custa o resultado (vfat, FUSE, SMB)
    try:
        os.replace(tmp, destino)
    except PermissionError as erro:
        # No Windows, substituir um arquivo ABERTO por outro programa levanta
        # PermissionError -- e o caso provavel e o pesquisador com a tabela
        # aberta no Excel, rodando o pipeline de novo. No Linux a troca
        # funciona, e este caminho nunca aparece.
        raise OSError(
            f"nao foi possivel publicar {destino}: o arquivo parece estar "
            f"aberto em outro programa (no Windows, o Excel segura o arquivo "
            f"enquanto a planilha esta aberta). Feche-o e rode de novo; o "
            f"resultado anterior esta intacto. ({erro})") from erro


def conferir_publicavel(destinos) -> None:
    """Levanta se algum destino nao pode ser substituido AGORA.

    A etapa publica varias saidas de uma vez. Se o travamento aparecer no MEIO
    do laco, umas ja foram publicadas e outras nao -- e isso deixa tabela sem
    manifesto, que e exatamente o que a publicacao atomica existe para
    impedir. Entao se confere TODOS antes de trocar QUALQUER um.

    Abrir para escrita e o unico teste honesto: no Windows a permissao do
    arquivo nao diz nada sobre outro processo estar segurando-o."""
    for destino in destinos:
        destino = Path(destino)
        if not destino.exists():
            continue  # destino novo: nada para travar
        try:
            with open(destino, "r+b"):
                pass
        except OSError as erro:
            raise OSError(
                f"nao e possivel substituir {destino}: o arquivo parece estar "
                f"aberto em outro programa (no Windows, o Excel segura o "
                f"arquivo enquanto a planilha esta aberta). Feche-o e rode de "
                f"novo; NADA foi publicado. ({erro})") from erro


def escrever_texto(destino, texto: str) -> None:
    """Grava `texto` em `destino`; falha no meio deixa o anterior intacto."""
    tmp = temporario(destino)
    try:
        tmp.write_text(texto, encoding="utf-8")
        publicar(tmp, destino)
    finally:
        tmp.unlink(missing_ok=True)  # so o temporario, nunca o arquivo do usuario
