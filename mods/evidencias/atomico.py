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
    os.replace(tmp, destino)


def escrever_texto(destino, texto: str) -> None:
    """Grava `texto` em `destino`; falha no meio deixa o anterior intacto."""
    tmp = temporario(destino)
    try:
        tmp.write_text(texto, encoding="utf-8")
        publicar(tmp, destino)
    finally:
        tmp.unlink(missing_ok=True)  # so o temporario, nunca o arquivo do usuario
