# -*- coding: utf-8 -*-
"""
T8 -- o README diz o que o software faz, inclusive o painel.

Um README que descreve um painel que nao existe (ou que omite o que existe) e
a mesma falha: quem le acredita numa coisa e encontra outra. Estes testes
prendem as afirmacoes que, se ficarem erradas, mandam o pesquisador para o
lugar errado.
"""
from mods import caminhos

README = (caminhos.RAIZ / "README.md").read_text(encoding="utf-8")


def test_readme_documenta_como_abrir_o_painel():
    assert "streamlit run painel.py" in README


def test_readme_diz_que_operar_escreve_e_analisar_le():
    """A invariante que separa as abas. Sem ela documentada, alguem vai
    procurar o botao de rodar na aba de analise."""
    baixo = README.lower()
    assert "operar" in baixo and "analisar" in baixo
    assert "só lê" in baixo or "so le" in baixo


def test_readme_diz_que_a_execucao_sobrevive_ao_fechamento_da_aba():
    """E a propriedade que o PID compra, e a unica que o pesquisador precisa
    saber para nao perder uma execucao de horas fechando o navegador."""
    baixo = README.lower()
    assert "fechar" in baixo or "fechamento" in baixo


def test_readme_diz_onde_fica_o_estado_da_execucao():
    """Quem precisar destravar o painel tem de saber qual arquivo apagar."""
    assert "tsv/.execucao" in README


def test_readme_nao_chama_o_painel_de_plano_futuro():
    """O painel existe agora. `(plano D)` no README mandaria o leitor procurar
    uma versao futura do que ele tem na frente."""
    assert "(plano D)" not in README


def test_readme_descreve_o_retomar_como_ele_e_hoje():
    """`--retomar` REFAZ a etapa cuja saida nao confere; ele nao aborta a
    execucao. A revisao do ramo anterior mudou esse comportamento, e um README
    que ainda diga 'recusa' manda o pesquisador recomecar do zero."""
    # Ancora no lugar onde `--retomar` e EXPLICADO, nao na primeira vez que a
    # flag aparece: uma janela de N caracteres depois da primeira ocorrencia
    # quebra sempre que algo e inserido antes dela, e o que importa aqui e o
    # comportamento descrito, nao a posicao no arquivo.
    inicio = README.index("confere a entrada")
    trecho = README[inicio:inicio + 900].lower()
    assert "refaz" in trecho or "roda de novo" in trecho


def test_o_painel_e_citado_na_arvore():
    assert "painel.py" in README


def test_toda_tela_do_painel_esta_descrita_no_readme():
    """Cada tela de Analisar aparece no README, e o README nao promete tela que
    nao existe. Antes este teste extraia a tabela `TELAS` e nunca a comparava
    com nada -- o docstring prometia uma verificacao que o corpo nao fazia."""
    import re
    fonte = (caminhos.RAIZ / "painel.py").read_text(encoding="utf-8")
    telas = re.findall(r'^\s{4}"([^"]+)": \("(?:etapa1|portao|ranking)"',
                       fonte, re.M)
    assert telas, "nao achei a tabela TELAS no painel"
    secao = README[README.index("## O painel"):]
    secao = secao[:secao.index("\n## ", 3)] if "\n## " in secao[3:] else secao
    # O painel nomeia as telas sem acento (o resto do codigo tambem); o README
    # e escrito com acento. A comparacao e sobre o texto sem acento.
    import unicodedata

    def simples(texto: str) -> str:
        return "".join(c for c in unicodedata.normalize("NFD", texto.lower())
                       if unicodedata.category(c) != "Mn")
    plano = simples(secao)
    sem_mencao = [t for t in telas
                  if simples(t).split()[0].rstrip("s") not in plano]
    assert not sem_mencao, f"telas sem mencao no README: {sem_mencao}"


def test_readme_diz_como_encerrar_uma_execucao():
    """O README mandava destravar apagando `rodando.pid` e nunca mencionava o
    botao "Encerrar". Apagar a tranca com o `main.py` vivo e o caminho
    documentado para ter duas execucoes sobre as mesmas saidas -- o que a
    tranca existe para impedir."""
    assert "Encerrar" in README


# ===========================================================================
# Espectro de alvos: o README diz o que a decisao 2 mudou
# ===========================================================================
def test_readme_diz_que_alvo_fora_do_modelo_entra_marcado():
    assert "no_modelo" in README


def test_readme_diz_onde_esta_o_teto_do_modelo():
    """O numero que justifica retreinar o preditor tem endereco: a chave
    `teto_do_modelo` do manifesto da etapa 1. Sem o endereco, ele continua
    sendo citado de execucao antiga."""
    assert "teto_do_modelo" in README


def test_readme_nao_promete_evidencia_prevista_fora_do_modelo():
    """A honestidade que a mudanca exige: o alvo entra, recebe evidencia
    MEDIDA e NAO recebe prevista. Prometer previsao seria pior que o descarte
    antigo, porque zero previsoes pareceria zero atividade."""
    inicio = README.index("O espectro de alvos tem um teto")
    fim = README.index("teto_do_modelo", inicio) + 1200
    trecho = README[inicio:fim].lower()
    assert "medida" in trecho and "prevista" in trecho
    assert "nao foi pontuado" in trecho or "não foi pontuado" in trecho


def test_readme_nao_diz_mais_que_o_alvo_do_reactome_e_descartado():
    """A limitacao antiga afirmava um teto que a decisao 2 mudou de lugar: o
    alvo nao e mais descartado, ele entra sem previsao."""
    secao = README[README.index("## Limitações conhecidas"):]
    assert "só pode receber evidência **medida**" not in secao or \
        "entra na lista" in secao


# ===========================================================================
# O README como pagina de entrada do repositorio
# ===========================================================================
def test_o_readme_abre_dizendo_o_que_o_pipeline_responde():
    """Quem chega pelo GitHub nao sabe o que isto e. A primeira tela tem de
    dizer a pergunta que o pipeline responde."""
    inicio = README[:1200].lower()
    assert "lista de compostos" in inicio and "alvos" in inicio


def test_o_readme_tem_um_mapa_da_arquitetura():
    """Nao basta descrever as etapas em prosa: quem avalia o repositorio quer
    ver o fluxo de dados antes de ler."""
    assert "## Arquitetura" in README
    assert "```mermaid" in README


def test_o_mapa_nomeia_as_tres_etapas_e_os_arquivos_reais():
    mapa = README[README.index("## Arquitetura"):]
    mapa = mapa[:mapa.index("\n## ", 3)]
    for esperado in ("coletar_evidencias.py", "aplicar_portao.py",
                     "ranquear.py", "evidencias.tsv", "portao.tsv",
                     "ranking.tsv", "manifesto.json"):
        assert esperado in mapa, esperado


def test_o_readme_tem_um_inicio_rapido_de_clone_limpo():
    assert "## Início rápido" in README or "## Inicio rapido" in README
    bloco = README[README.lower().index("início rápido"):][:1500]
    assert "pip install -r requirements.txt" in bloco
    assert "streamlit run painel.py" in bloco or "python main.py" in bloco


def test_o_readme_documenta_o_acervo_do_openfda():
    """O plano C: sem isto o criterio que mais reprova sai indeterminado para
    99,5% da base, e quem clona nao sabe que ha remedio."""
    assert "openfda_massa.py" in README
    assert "--acervo-bula" in README


def test_o_readme_tem_secao_de_windows():
    assert "## Windows" in README
    bloco = README[README.index("## Windows"):][:1200]
    assert "pytest" in bloco


def test_o_readme_diz_o_que_nao_vai_no_repositorio():
    """Modelo e dados ficam fora do git. Quem clona precisa saber disso antes
    de rodar e achar que esta quebrado."""
    assert "fora do git" in README.lower()


def test_o_readme_mantem_as_medicoes():
    """O valor do documento sao os numeros medidos. Nenhuma reorganizacao pode
    perde-los."""
    for numero in ("0,484", "0,384", "107.123", "0 de 6", "99,5"):
        assert numero in README, numero


def test_o_readme_declara_as_licencas():
    assert "## Licenças" in README or "## Licencas" in README


def test_o_readme_distingue_o_teto_de_taxa_do_teto_de_COBERTURA():
    """MEDIDO com o acervo real, 60 ligantes de `exemplos/`: 58 saem
    `indeterminado` por **nao ter bula nenhuma** no openFDA, e 2 passam. O
    acervo remove o teto de TAXA; o teto de COBERTURA fica, e para uma base de
    produto natural ele e o que manda -- openFDA e rotulagem de farmaco dos
    EUA.

    Dizer so "o acervo remove o teto" prometeria ao pesquisador um criterio
    atuante que ele nao vai ter. `indeterminado` por falta de DOCUMENTO e a
    resposta certa, e tem de estar escrito."""
    assert "sem bula no openFDA" in README or "cobertura" in README.lower()
    i = README.lower().index("cobertura")
    trecho = README[max(0, i - 600):i + 900].lower()
    assert "indeterminado" in trecho
