# Etapa 1 — Coletar evidências · Plano de implementação

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produzir `evidencias.tsv` — uma linha por evidência de par ligante/alvo, medida ou prevista, rastreável — a partir de sintomas, alvos (do usuário, enriquecidos por Open Targets e expandidos por Reactome) e ligantes no padrão do prospector, com painel Streamlit e manifesto reproduzível.

**Architecture:** Um pacote `mods/evidencias/` com um módulo por responsabilidade, orquestrado por uma CLI autônoma `mods/coletar_evidencias.py`. A CLI fica executável já na tarefa 4 e cada tarefa seguinte acrescenta uma fonte de evidência sem quebrá-la — o usuário só consegue dar retorno rodando o pipeline completo, então ele precisa rodar cedo. Escrita em fluxo, porque a base real tem ~100 mil ligantes.

**Tech Stack:** Python 3, `requests`, `rdkit`, `torch`, `pandas`, `streamlit`, `pytest` (todos já instalados). Reaproveita `mods/prever_alvos.py` (previsão e rótulos) e `mods/importar_prospector.py` (química de identidade).

**Spec:** `docs/superpowers/specs/2026-11-25-etapa1-coletar-evidencias-design.md`

## Global Constraints

- **Formato:** TSV na entrada e na saída. `.csv` aceito na leitura, com aviso. 42,8% dos nomes de produto natural contêm vírgula.
- **Separador decimal:** ponto em toda coluna de máquina. A vírgula só aparece no texto de exibição de `atividade_relatada`.
- **Procedência nunca se mistura:** nenhuma coluna pode conter tanto `pchembl` quanto probabilidade de modelo. `pchembl` é a única coluna comparável entre linhas.
- **`torch.load` sempre com `weights_only=True`** — um `.pt` é um pickle.
- **Escala:** projetar para 108.935 ligantes. Escrita em fluxo; nada que materialize todas as linhas em memória.
- **Corte padrão de linhas previstas:** `--corte 0.1`, registrado no manifesto.
- **Teto de via do Reactome:** 150 proteínas.
- **Licença:** `LICENSE` MIT na raiz; `dados/LICENSE` CC BY-SA 3.0. Nunca usar os arquivos de *Functional Interaction* do Reactome.
- **Enriquecimento só acrescenta:** nenhum alvo do usuário é removido ou substituído, em nenhuma circunstância.
- **Nenhum import do `prospector`.** Só o formato dos dados dele é conhecido.

## Review Focus

Classes de entrada que o desenho implica mas que nenhuma tarefa exercitaria sozinha. Cada uma tem o seu teste atribuído à tarefa dona do código.

1. **Nome de ligante contendo tabulação** (não só vírgula) quebraria as colunas do TSV em silêncio. Esperado: escapar na escrita e devolver o mesmo texto na leitura. → Tarefa 2.
2. **Mesmo ligante repetido na entrada** com nomes diferentes e o mesmo `inchikey_parental` multiplicaria todas as linhas por alvo. Esperado: deduplicar pela parental, manter o primeiro nome, registrar os demais. → Tarefa 3.
3. **Alvo pedido por `alvo`+`organismo` que casa com mais de um `target_chembl_id`** (ortólogo, complexo, grupo de seletividade) seria escolhido em silêncio. Esperado: erro explícito listando os candidatos. → Tarefa 6.
4. **Registro com `unidade = '%'` e `pchembl` vazio** faria `valor` ser comparado entre unidades. Esperado: `pchembl` vazio, `valor`/`unidade` preenchidos, nenhuma ordenação cruzando unidades. → Tarefa 7.
5. **Open Targets ou Reactome fora do ar no meio da execução** deixaria saída parcial que parece completa. Esperado: falhar antes de escrever a primeira linha, ou marcar o enriquecimento como incompleto no manifesto. → Tarefa 9.

---

### Task 1: Licença, dependências e esqueleto de testes

Sem `LICENSE` o repositório público não autoriza reuso — é requisito direto de "disponível no GitHub e citável", e não depende de nenhum código.

**Files:**
- Create: `LICENSE`, `dados/LICENSE`, `requirements.txt`, `tests/test_licenca.py`

**Interfaces:**
- Consumes: nada.
- Produces: `pytest` roda no repositório; `requirements.txt` lista as dependências.

- [ ] **Step 1: Write the failing test** em `tests/test_licenca.py`

```python
from pathlib import Path
RAIZ = Path(__file__).resolve().parent.parent

def test_licenca_mit_na_raiz():
    texto = (RAIZ / "LICENSE").read_text(encoding="utf-8")
    assert "MIT License" in texto

def test_dados_derivados_do_chembl_sao_share_alike():
    texto = (RAIZ / "dados" / "LICENSE").read_text(encoding="utf-8")
    assert "CC BY-SA 3.0" in texto
    assert "ChEMBL" in texto

def test_requirements_lista_as_dependencias():
    texto = (RAIZ / "requirements.txt").read_text(encoding="utf-8").lower()
    for pacote in ("rdkit", "torch", "requests", "pandas", "streamlit", "pytest"):
        assert pacote in texto
```

- [ ] **Step 2: Run test to verify it fails** — `pytest tests/test_licenca.py -v`. Esperado: FAIL, `FileNotFoundError`.

- [ ] **Step 3: Criar os três arquivos.** `LICENSE` é o texto MIT padrão, detentor `Roberto Alexandre Filho`, ano 2026. `dados/LICENSE` declara que os arquivos daquela pasta derivam do ChEMBL e estão sob CC BY-SA 3.0, com a citação do ChEMBL. `requirements.txt` com os pacotes acima, sem fixar versão salvo `numpy>=2` (o `pp4ckd` pede `<2`, e este repositório deliberadamente não o importa).

- [ ] **Step 4: Run test to verify it passes** — `pytest tests/test_licenca.py -v`. Esperado: 3 passed.

- [ ] **Step 5: Commit**

```bash
git add LICENSE dados/LICENSE requirements.txt tests/test_licenca.py
git commit -m "Add MIT license, data license and dependency list"
```

---

### Task 2: Esquema e escrita em fluxo de `evidencias.tsv`

**Files:**
- Create: `mods/evidencias/__init__.py`, `mods/evidencias/esquema.py`, `tests/test_esquema.py`

**Interfaces:**
- Consumes: nada.
- Produces:
  - `COLUNAS_EVIDENCIA: list[str]` — ordem exata das colunas.
  - `class EscritorEvidencias` — gerenciador de contexto; `__init__(self, caminho: str | Path)`, `escrever(self, linha: dict) -> None`, escreve o cabeçalho (prefixado por `#`) na abertura e **não acumula linhas em memória**.
  - `ler_evidencias(caminho: str | Path) -> Iterator[dict]`
  - `escapar(valor) -> str` — normaliza um campo para TSV.
  - `formatar_atividade(medida: str, relacao: str, valor, unidade: str) -> str` — texto de exibição, ex. `IC50 = 120 nM`.

As colunas, nesta ordem: `ligante`, `ligante_inchikey_parental`, `ligante_chembl_id`, `alvo_chembl_id`, `alvo`, `organismo`, `tipo_alvo`, `origem_alvo`, `tipo_evidencia`, `fonte`, `atividade_relatada`, `medida`, `relacao`, `valor`, `unidade`, `pchembl`, `tipo_ensaio`, `ensaio_descricao`, `triagem_hts`, `confianca_ensaio`, `n_moleculas_treino`, `documento`, `pmid`, `ano`, `procedencia`.

- [ ] **Step 1: Write the failing test** em `tests/test_esquema.py`

```python
def test_cabecalho_tem_as_colunas_na_ordem(tmp_path):
    destino = tmp_path / "e.tsv"
    with EscritorEvidencias(destino) as esc:
        esc.escrever({"ligante": "quercetina", "alvo_chembl_id": "CHEMBL230"})
    primeira = destino.read_text(encoding="utf-8").splitlines()[0]
    assert primeira == "#" + "\t".join(COLUNAS_EVIDENCIA)

def test_nome_com_virgula_e_tabulacao_sobrevive_a_ida_e_volta(tmp_path):
    # Review Focus 1: 42,8% dos nomes têm vírgula, e tabulação quebraria a linha.
    nome = "1,8-cineol\tsin.: eucaliptol"
    destino = tmp_path / "e.tsv"
    with EscritorEvidencias(destino) as esc:
        esc.escrever({"ligante": nome, "alvo_chembl_id": "CHEMBL230"})
    linhas = list(ler_evidencias(destino))
    assert len(linhas) == 1
    assert linhas[0]["ligante"] == nome
    assert linhas[0]["alvo_chembl_id"] == "CHEMBL230"

def test_separador_decimal_e_ponto():
    assert formatar_atividade("IC50", "=", 120.5, "nM") == "IC50 = 120.5 nM"

def test_campo_ausente_vira_vazio_nao_a_palavra_none(tmp_path):
    destino = tmp_path / "e.tsv"
    with EscritorEvidencias(destino) as esc:
        esc.escrever({"ligante": "x"})
    assert "None" not in destino.read_text(encoding="utf-8")
```

- [ ] **Step 2: Run test to verify it fails** — `pytest tests/test_esquema.py -v`. Esperado: FAIL, `ModuleNotFoundError`.

- [ ] **Step 3: Implementar `esquema.py`.** `escapar` troca tabulação, retorno de carro e quebra de linha por espaço — nunca descarta o campo; se mudar o texto, isso é visível no teste de ida e volta, então a escolha aqui é **escapar preservando o conteúdo** (sequência de escape própria, revertida na leitura), não substituir.

- [ ] **Step 4: Run test to verify it passes** — `pytest tests/test_esquema.py -v`. Esperado: 4 passed.

- [ ] **Step 5: Commit**

```bash
git add mods/evidencias/ tests/test_esquema.py
git commit -m "Add evidence table schema with streaming TSV writer"
```

---

### Task 3: Leitura e validação das três entradas

**Files:**
- Create: `mods/evidencias/entrada.py`, `tests/test_entrada.py`, `tests/dados/ligantes_minimo.tsv`, `tests/dados/alvos_minimo.tsv`, `tests/dados/sintomas.txt`

**Interfaces:**
- Consumes: `mods.importar_prospector.forma_parental`, `inchikey_de`, `molecula`.
- Produces:
  - `@dataclass Problema: linha: int, campo: str, mensagem: str, gravidade: Literal["erro", "aviso"]`
  - `ler_ligantes(caminho) -> tuple[list[dict], list[Problema]]` — exige `nome` e `smiles`; recalcula `inchikey` e `inchikey_parental` quando ausentes; deduplica por `inchikey_parental`.
  - `ler_alvos(caminho) -> tuple[list[dict], list[Problema]]` — exige `alvo_chembl_id` **ou** `alvo`+`organismo`.
  - `ler_sintomas(caminho) -> list[str]` — um por linha, `#` comenta.

- [ ] **Step 1: Write the failing test** em `tests/test_entrada.py`

```python
def test_ligante_sem_smiles_e_erro():
    ligantes, problemas = ler_ligantes("tests/dados/ligantes_minimo.tsv")
    assert any(p.gravidade == "erro" and p.campo == "smiles" for p in problemas)

def test_inchikey_parental_e_recalculado_quando_ausente():
    ligantes, _ = ler_ligantes("tests/dados/ligantes_minimo.tsv")
    aspirina = next(l for l in ligantes if l["nome"] == "aspirina")
    assert aspirina["inchikey_parental"] == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"

def test_ligante_repetido_pela_parental_e_deduplicado():
    # Review Focus 2: sal e base livre, nomes diferentes, mesma molécula.
    ligantes, problemas = ler_ligantes("tests/dados/ligantes_minimo.tsv")
    chaves = [l["inchikey_parental"] for l in ligantes]
    assert len(chaves) == len(set(chaves))
    assert any("duplicad" in p.mensagem for p in problemas)

def test_alvo_por_nome_sem_organismo_e_erro():
    alvos, problemas = ler_alvos("tests/dados/alvos_minimo.tsv")
    assert any(p.gravidade == "erro" and p.campo == "organismo" for p in problemas)

def test_sintomas_ignora_comentario():
    assert ler_sintomas("tests/dados/sintomas.txt") == ["neuropathic pain", "peripheral neuropathy"]
```

- [ ] **Step 2: Criar as fixtures e rodar o teste.** `ligantes_minimo.tsv` tem: aspirina só com `nome`+`smiles`; uma linha sem `smiles`; cloridrato de duloxetina e duloxetina base com nomes diferentes. `alvos_minimo.tsv` tem uma linha com `alvo_chembl_id`, uma com `alvo`+`organismo`, e uma com `alvo` e organismo vazio. Run: `pytest tests/test_entrada.py -v`. Esperado: FAIL.

- [ ] **Step 3: Implementar `entrada.py`.** A deduplicação mantém o primeiro nome visto e emite um `Problema` de gravidade `aviso` nomeando os descartados. Leitor aceita `.tsv` e `.csv`; se o arquivo for `.csv` e alguma célula contiver vírgula não citada, emite aviso.

- [ ] **Step 4: Run test to verify it passes** — `pytest tests/test_entrada.py -v`. Esperado: 5 passed.

- [ ] **Step 5: Commit**

```bash
git add mods/evidencias/entrada.py tests/test_entrada.py tests/dados/
git commit -m "Read and validate the three stage-1 inputs"
```

---

### Task 4: CLI e manifesto — primeiro ponto executável

Depois desta tarefa o pipeline roda ponta a ponta e produz arquivos válidos, ainda sem nenhuma linha de evidência. É o ponto em que o usuário consegue executar e reagir.

**Files:**
- Create: `mods/coletar_evidencias.py`, `mods/evidencias/manifesto.py`, `tests/test_cli.py`

**Interfaces:**
- Consumes: Tarefas 2 e 3.
- Produces:
  - `escrever_manifesto(caminho, *, parametros: dict, versoes: dict, contagens: dict) -> None`
  - `ler_manifesto(caminho) -> dict`
  - `main(argv: list[str] | None = None) -> int` em `coletar_evidencias.py`, com `--ligantes`, `--alvos`, `--sintomas`, `--saida`, `--manifesto`, `--corte`, `--limite`.

- [ ] **Step 1: Write the failing test** em `tests/test_cli.py`

```python
def test_execucao_minima_produz_tabela_valida_e_manifesto(tmp_path):
    saida = tmp_path / "evidencias.tsv"
    codigo = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
                   "--alvos", "tests/dados/alvos_minimo.tsv",
                   "--sintomas", "tests/dados/sintomas.txt",
                   "--saida", str(saida), "--sem-rede"])
    assert codigo == 0
    assert saida.read_text(encoding="utf-8").splitlines()[0].startswith("#ligante\t")
    m = ler_manifesto(tmp_path / "manifesto.json")
    assert m["parametros"]["corte"] == 0.1
    assert "pipeline_commit" in m["versoes"]

def test_erro_de_entrada_interrompe_antes_de_escrever(tmp_path):
    saida = tmp_path / "evidencias.tsv"
    codigo = main(["--ligantes", "tests/dados/ligantes_sem_coluna.tsv",
                   "--alvos", "tests/dados/alvos_minimo.tsv",
                   "--sintomas", "tests/dados/sintomas.txt",
                   "--saida", str(saida), "--sem-rede"])
    assert codigo != 0
    assert not saida.exists()
```

- [ ] **Step 2: Run test to verify it fails** — `pytest tests/test_cli.py -v`. Esperado: FAIL.

- [ ] **Step 3: Implementar a CLI e o manifesto.** `--sem-rede` pula enriquecimento e evidência medida por API — é o que torna a tarefa testável sozinha e permanece útil depois. O manifesto grava `pipeline_commit` via `git rev-parse HEAD`, a data, os parâmetros e as contagens. Validação de entrada roda inteira **antes** de abrir o arquivo de saída.

- [ ] **Step 4: Run test to verify it passes** — `pytest tests/test_cli.py -v`. Esperado: 2 passed.

- [ ] **Step 5: Commit**

```bash
git add mods/coletar_evidencias.py mods/evidencias/manifesto.py tests/
git commit -m "Add stage-1 CLI and run manifest"
```

---

### Task 5: Evidência prevista (PP4CKD)

Primeira fonte real de linhas, e funciona sem rede — só precisa dos pesos.

**Files:**
- Create: `mods/evidencias/predito.py`, `tests/test_predito.py`
- Modify: `mods/coletar_evidencias.py` (ligar a fonte)

**Interfaces:**
- Consumes: `mods.prever_alvos.carregar_modelo`, `carregar_rotulos`, `fingerprint_de`; Tarefa 2.
- Produces: `linhas_preditas(ligantes: list[dict], alvos: list[dict], *, checkpoint: str, rotulos: list[dict], corte: float = 0.1) -> Iterator[dict]`

- [ ] **Step 1: Write the failing test** em `tests/test_predito.py` (usa um modelo falso de pesos aleatórios gravado com `torch.save`, não os pesos reais, para o teste rodar sem download)

```python
def test_respeita_o_corte(modelo_falso, rotulos_falsos):
    linhas = list(linhas_preditas(LIGANTES, ALVOS, checkpoint=modelo_falso,
                                  rotulos=rotulos_falsos, corte=0.99))
    assert all(float(l["pchembl"] or 0) == 0 for l in linhas)
    assert all(l["tipo_evidencia"] == "predito" for l in linhas)

def test_alvo_fora_do_modelo_nao_gera_linha(modelo_falso, rotulos_falsos):
    alvos = ALVOS + [{"alvo_chembl_id": "CHEMBL_INEXISTENTE", "no_modelo": False}]
    linhas = list(linhas_preditas(LIGANTES, alvos, checkpoint=modelo_falso,
                                  rotulos=rotulos_falsos, corte=0.0))
    assert not any(l["alvo_chembl_id"] == "CHEMBL_INEXISTENTE" for l in linhas)

def test_nunca_preenche_pchembl(modelo_falso, rotulos_falsos):
    linhas = list(linhas_preditas(LIGANTES, ALVOS, checkpoint=modelo_falso,
                                  rotulos=rotulos_falsos, corte=0.0))
    assert all(l["pchembl"] == "" for l in linhas)
    assert all(l["n_moleculas_treino"] != "" for l in linhas)
```

- [ ] **Step 2: Run test to verify it fails** — `pytest tests/test_predito.py -v`. Esperado: FAIL.

- [ ] **Step 3: Implementar `predito.py`** e ligá-la na CLI atrás de `--checkpoint`. `atividade_relatada` recebe `p = 0.93 (ECFP4)`; `procedencia` recebe `predito (PP4CKD ECFP4 <nome do checkpoint>)`. Gera em lotes e entrega por `yield`, sem materializar.

- [ ] **Step 4: Run test to verify it passes** — `pytest tests/test_predito.py -v`. Esperado: 3 passed.

- [ ] **Step 5: Commit**

```bash
git add mods/evidencias/predito.py tests/test_predito.py mods/coletar_evidencias.py
git commit -m "Add predicted evidence rows from PP4CKD"
```

---

### Task 6: Resolução no ChEMBL — alvo por UniProt, ligante por InChIKey

**Files:**
- Create: `mods/evidencias/chembl.py`, `tests/test_chembl.py`, `tests/dados/respostas/` (respostas de rede gravadas)

**Interfaces:**
- Consumes: Tarefa 9 é posterior, então esta tarefa recebe um `obter_json: Callable[[str], dict]` injetável, cujo padrão usa `requests`.
- Produces:
  - `alvos_por_uniprot(acc: str, organismo: str, *, obter_json=...) -> list[dict]` — ordenados preferindo `SINGLE PROTEIN`.
  - `resolver_alvo_por_nome(nome: str, organismo: str, *, obter_json=...) -> dict` — levanta `AlvoAmbiguo` quando há mais de um candidato.
  - `resolver_ligante(inchikey: str, inchikey_parental: str, *, obter_json=...) -> str | None` — tenta InChIKey exato e depois o bloco de conectividade (14 primeiros caracteres).
  - `uniprot_do_alvo(chembl_id: str, *, obter_json=...) -> str | None` — o caminho inverso, necessário porque o alvo do usuário chega como ChEMBL ID e a expansão do Reactome parte de UniProt.
  - `classe_de_proteina(chembl_id: str, *, obter_json=...) -> str | None` — a classificação que o colapso de família usa como definição de "família".
  - `class AlvoAmbiguo(Exception)` com atributo `candidatos: list[dict]`.

- [ ] **Step 1: Write the failing test** em `tests/test_chembl.py`

```python
def test_uniprot_mapeia_para_o_single_protein(respostas):
    achados = alvos_por_uniprot("P23975", "Homo sapiens", obter_json=respostas)
    assert achados[0]["target_chembl_id"] == "CHEMBL222"
    assert achados[0]["target_type"] == "SINGLE PROTEIN"

def test_alvo_ambiguo_levanta_com_os_candidatos(respostas):
    # Review Focus 3: nunca escolher sozinho entre ortólogo, complexo e grupo.
    with pytest.raises(AlvoAmbiguo) as erro:
        resolver_alvo_por_nome("noradrenaline transporter", "", obter_json=respostas)
    assert len(erro.value.candidatos) > 1

def test_ligante_sem_chembl_id_resolve_pelo_bloco_de_conectividade(respostas):
    # O ChEMBL indexa por forma salina; a parental só encontra pelo bloco.
    assert resolver_ligante("", "ZEUITGRIYCTCEM-KRWDZBQOSA-N", obter_json=respostas) is not None

def test_caminho_inverso_alvo_para_uniprot(respostas):
    assert uniprot_do_alvo("CHEMBL222", obter_json=respostas) == "P23975"

def test_classe_de_proteina_para_o_colapso(respostas):
    assert classe_de_proteina("CHEMBL222", obter_json=respostas)
```

- [ ] **Step 2: Gravar as respostas de rede** em `tests/dados/respostas/` (um JSON por URL consultada, gerado uma vez com `curl` contra a API real) e rodar o teste. Run: `pytest tests/test_chembl.py -v`. Esperado: FAIL.

- [ ] **Step 3: Implementar `chembl.py`.** Endpoint de alvo: `/chembl/api/data/target.json?target_components__accession=<acc>&organism=<org>`. Endpoint de molécula: `/chembl/api/data/molecule.json?molecule_structures__standard_inchi_key=<ik>` e, para o bloco, `__istartswith`.

- [ ] **Step 4: Run test to verify it passes** — `pytest tests/test_chembl.py -v`. Esperado: 5 passed.

- [ ] **Step 5: Commit**

```bash
git add mods/evidencias/chembl.py tests/test_chembl.py tests/dados/respostas/
git commit -m "Resolve ChEMBL targets by UniProt and ligands by InChIKey"
```

---

### Task 7: Evidência medida a partir de `evidencias_chembl.csv`

**Files:**
- Create: `mods/evidencias/medido.py`, `tests/test_medido.py`, `dados/evidencias_controle.tsv`
- Modify: `mods/coletar_evidencias.py`

**Interfaces:**
- Consumes: Tarefas 2, 3, 6.
- Produces: `linhas_medidas(caminho_csv: str, ligantes: list[dict], alvos: list[dict]) -> Iterator[dict]`

- [ ] **Step 1: Write the failing test** em `tests/test_medido.py`

```python
def test_registro_em_porcentagem_nao_preenche_pchembl():
    # Review Focus 4: 120 nM e 78% não se ordenam; pchembl é a única coluna comparável.
    linhas = list(linhas_medidas("tests/dados/evidencias_amostra.csv", LIGANTES, ALVOS))
    pct = [l for l in linhas if l["unidade"] == "%"]
    assert pct and all(l["pchembl"] == "" for l in pct)
    assert all(l["valor"] != "" for l in pct)

def test_triagem_hts_e_marcada():
    linhas = list(linhas_medidas("tests/dados/evidencias_amostra.csv", LIGANTES, ALVOS))
    assert any(l["triagem_hts"] == "sim" for l in linhas)
    assert any(l["triagem_hts"] == "nao" for l in linhas)

def test_relacao_de_limite_inferior_entra_marcada():
    linhas = list(linhas_medidas("tests/dados/evidencias_amostra.csv", LIGANTES, ALVOS))
    assert any(l["relacao"] == ">" for l in linhas)

def test_uma_linha_por_ensaio_nao_por_par():
    linhas = list(linhas_medidas("tests/dados/evidencias_amostra.csv", LIGANTES, ALVOS))
    pares = [(l["ligante"], l["alvo_chembl_id"]) for l in linhas]
    assert len(pares) > len(set(pares))
```

- [ ] **Step 2: Criar `tests/dados/evidencias_amostra.csv`** (extrato de ~30 linhas reais do `evidencias_chembl.csv`, cobrindo nM, `%`, relação `>`, uma linha `PUBCHEM_BIOASSAY` e um par com três ensaios) e rodar. Run: `pytest tests/test_medido.py -v`. Esperado: FAIL.

- [ ] **Step 3: Implementar `medido.py`.** Casa o alvo por `(alvo, organismo_alvo)` e o ligante pela `InChIKey`. `triagem_hts` reusa a heurística de `PUBCHEM_BIOASSAY`/`qHTS` que já existe em `mods/validar_alvos_medidos.py`. Gravar também `dados/evidencias_controle.tsv`: o resultado para o conjunto-controle, versionado sob CC BY-SA 3.0, para os testes de aceite rodarem sem rede.

- [ ] **Step 4: Run test to verify it passes** — `pytest tests/test_medido.py -v`. Esperado: 4 passed.

- [ ] **Step 5: Commit**

```bash
git add mods/evidencias/medido.py tests/ dados/evidencias_controle.tsv mods/coletar_evidencias.py
git commit -m "Add measured evidence rows from ChEMBL bioactivity"
```

---

### Task 8: Enriquecimento — Open Targets e expansão por via do Reactome

**Files:**
- Create: `mods/evidencias/open_targets.py`, `mods/evidencias/reactome.py`, `tests/test_enriquecimento.py`
- Modify: `mods/coletar_evidencias.py`

**Interfaces:**
- Consumes: Tarefa 6 (`alvos_por_uniprot`); `obter_json` injetável.
- Produces:
  - `buscar_efo(sintoma: str, *, obter_json=...) -> list[dict]` — candidatos com `id` e `name`, **nunca escolhe sozinho**.
  - `alvos_associados(efo_id: str, *, corte: float, obter_json=...) -> list[dict]` — com `score`, `datatypes`, `uniprot`.
  - `colapsar_familia(alvos: list[dict]) -> list[dict]` — um representante por classe de proteína do ChEMBL; os demais em `familia_colapsados`.
  - `expandir_vias(sementes: list[dict], *, teto: int = 150, obter_json=...) -> list[dict]` — com `via_origem` e `alvo_semente`. Recebe as sementes já com `uniprot` preenchido por `uniprot_do_alvo` (Tarefa 6).
  - `escrever_alvos_enriquecidos(caminho, alvos: list[dict]) -> None` — grava a segunda saída da etapa, com as colunas do desenho: `alvo_chembl_id`, `alvo`, `organismo`, `tipo_alvo`, `gene`, `origem`, `efo_id`, `sintoma`, `escore_associacao`, `familia_colapsados`, `via_origem`, `alvo_semente`, `n_moleculas_treino`, `no_modelo`, `observacao`.

`colapsar_familia` exige que cada alvo já traga `classe`, preenchida por `classe_de_proteina` (Tarefa 6).

- [ ] **Step 1: Write the failing test** em `tests/test_enriquecimento.py`

```python
def test_efo_devolve_candidatos_sem_escolher(respostas):
    cands = buscar_efo("neuropathic pain", obter_json=respostas)
    assert cands[0]["id"] == "EFO_0005762"
    assert len(cands) > 1

def test_colapso_de_familia_guarda_os_descartados():
    alvos = [{"alvo_chembl_id": "A", "classe": "VGCC", "score": 0.605},
             {"alvo_chembl_id": "B", "classe": "VGCC", "score": 0.603},
             {"alvo_chembl_id": "C", "classe": "SLC6", "score": 0.583}]
    saida = colapsar_familia(alvos)
    assert [a["alvo_chembl_id"] for a in saida] == ["A", "C"]
    assert saida[0]["familia_colapsados"] == "B"

def test_via_grande_demais_e_descartada(respostas):
    # 'G alpha (i) signalling events' tem 317 proteínas; 'Opioid Signalling', 90.
    expandidos = expandir_vias([{"uniprot": "P35372", "alvo_chembl_id": "CHEMBL233"}],
                               teto=150, obter_json=respostas)
    vias = {a["via_origem"] for a in expandidos}
    assert "R-HSA-111885" in vias
    assert "R-HSA-418594" not in vias

def test_expansao_registra_a_semente(respostas):
    expandidos = expandir_vias([{"uniprot": "P35372", "alvo_chembl_id": "CHEMBL233"}],
                               teto=150, obter_json=respostas)
    assert all(a["alvo_semente"] == "CHEMBL233" for a in expandidos)
    assert all(a["origem"] == "reactome" for a in expandidos)

def test_alvo_do_usuario_nunca_desaparece_do_enriquecido(tmp_path):
    # Restrição global: o enriquecimento só acrescenta.
    destino = tmp_path / "alvos_enriquecidos.tsv"
    escrever_alvos_enriquecidos(destino, [
        {"alvo_chembl_id": "CHEMBL222", "origem": "usuario"},
        {"alvo_chembl_id": "CHEMBL233", "origem": "open_targets"},
    ])
    texto = destino.read_text(encoding="utf-8")
    assert "CHEMBL222" in texto and "usuario" in texto
```

- [ ] **Step 2: Gravar as respostas** do Open Targets e do Reactome em `tests/dados/respostas/` e rodar. Run: `pytest tests/test_enriquecimento.py -v`. Esperado: FAIL.

- [ ] **Step 3: Implementar os dois módulos e ligá-los na CLI.** Open Targets: GraphQL em `https://api.platform.opentargets.org/api/v4/graphql`, campos `search` (para EFO), `disease.associatedTargets` (com `score` e `datatypeScores`) e `target.proteinIds` (UniProt). Reactome: `GET /ContentService/data/mapping/UniProt/{acc}/pathways?species=9606` e `GET /ContentService/data/participants/{stId}/referenceEntities` — **nunca** os arquivos de Functional Interaction. A classe de proteína para o colapso vem do ChEMBL (`target_components__protein_classifications`). O alvo expandido só entra se existir no espaço do PP4CKD. **Medir e registrar no README quanto esse segundo filtro corta** — o desenho deixa isso em aberto de propósito.

- [ ] **Step 4: Run test to verify it passes** — `pytest tests/test_enriquecimento.py -v`. Esperado: 5 passed.

- [ ] **Step 5: Commit**

```bash
git add mods/evidencias/open_targets.py mods/evidencias/reactome.py tests/ mods/coletar_evidencias.py
git commit -m "Enrich targets from Open Targets and expand pathways from Reactome"
```

---

### Task 9: Cache de rede e modo de repetição do manifesto

**Files:**
- Create: `mods/evidencias/cache.py`, `tests/test_repeticao.py`
- Modify: `mods/evidencias/manifesto.py`, `mods/coletar_evidencias.py`

**Interfaces:**
- Consumes: Tarefas 4, 6, 8.
- Produces:
  - `class Cache` — `__init__(self, diretorio)`, `obter_json(self, url) -> dict` (serve do disco quando houver, senão busca e grava).
  - `verificar_manifesto(anterior: dict, atual: dict) -> list[str]` — divergências encontradas; lista vazia significa que pode repetir.

- [ ] **Step 1: Write the failing test** em `tests/test_repeticao.py`

```python
def test_repeticao_recusa_versao_divergente():
    anterior = {"versoes": {"open_targets": "24.09", "checkpoint_md5": "abc"}}
    atual = {"versoes": {"open_targets": "25.03", "checkpoint_md5": "abc"}}
    assert "open_targets" in " ".join(verificar_manifesto(anterior, atual))

def test_queda_de_rede_nao_deixa_saida_parcial(tmp_path, rede_que_falha):
    # Review Focus 5: falhar antes de escrever, nunca a meio arquivo.
    saida = tmp_path / "evidencias.tsv"
    codigo = main(["--ligantes", "tests/dados/ligantes_minimo.tsv",
                   "--alvos", "tests/dados/alvos_minimo.tsv",
                   "--sintomas", "tests/dados/sintomas.txt",
                   "--saida", str(saida)])
    assert codigo != 0
    assert not saida.exists()

def test_cache_serve_a_segunda_chamada_sem_rede(tmp_path):
    cache = Cache(tmp_path)
    primeira = cache.obter_json("https://exemplo/teste")   # grava
    segunda = cache.obter_json("https://exemplo/teste")    # sem rede
    assert primeira == segunda
```

- [ ] **Step 2: Run test to verify it fails** — `pytest tests/test_repeticao.py -v`. Esperado: FAIL.

- [ ] **Step 3: Implementar o cache e o modo de repetição.** A versão do Open Targets vem do campo `meta.dataVersion` da própria API. Escrita atômica: gravar num arquivo temporário e renomear só no fim, para que uma falha nunca deixe saída parcial. `--manifesto <arquivo>` entra em modo de repetição e aborta se `verificar_manifesto` devolver divergências.

- [ ] **Step 4: Run test to verify it passes** — `pytest tests/test_repeticao.py -v`. Esperado: 3 passed.

- [ ] **Step 5: Commit**

```bash
git add mods/evidencias/cache.py mods/evidencias/manifesto.py tests/test_repeticao.py mods/coletar_evidencias.py
git commit -m "Add network cache and manifest replay mode"
```

---

### Task 10: Painel Streamlit

**Files:**
- Create: `painel.py`, `mods/evidencias/resumo.py`, `tests/test_resumo.py`
- Modify: `README.md`

**Interfaces:**
- Consumes: todas as anteriores.
- Produces:
  - `resumo_ligantes(ligantes, problemas) -> dict` — contagens que o painel mostra no upload.
  - `resumo_alvos(alvos, rotulos) -> dict` — inclui `sem_organismo`, `fora_do_modelo`, `pouco_treino`.

As funções de resumo ficam fora do `painel.py` justamente para serem testáveis sem subir o Streamlit.

- [ ] **Step 1: Write the failing test** em `tests/test_resumo.py`

```python
def test_resumo_de_alvos_conta_os_problemas(rotulos_falsos):
    alvos = [{"alvo_chembl_id": "CHEMBL222", "organismo": "Homo sapiens"},
             {"alvo": "algum alvo", "organismo": ""},
             {"alvo_chembl_id": "CHEMBL_X", "organismo": "Homo sapiens"}]
    r = resumo_alvos(alvos, rotulos_falsos)
    assert r["sem_organismo"] == 1
    assert r["fora_do_modelo"] == 1

def test_resumo_de_ligantes_separa_validos_de_invalidos():
    r = resumo_ligantes([{"smiles": "CCO"}, {"smiles": "xyz!"}], [])
    assert r["validos"] == 1 and r["invalidos"] == 1
```

- [ ] **Step 2: Run test to verify it fails** — `pytest tests/test_resumo.py -v`. Esperado: FAIL.

- [ ] **Step 3: Implementar `resumo.py` e `painel.py`.** O painel: três uploads; resumo e avisos na hora do upload; botão "Enriquecer alvos" separado de "Coletar evidências", para o usuário editar a lista antes de gastar a coleta; EFO mostrado com as alternativas e trocável; execução em segundo plano com progresso; tabela filtrável e botões de descarregar. Caminho dos pesos e do `target_labels.tsv` na barra lateral, com variável de ambiente como padrão; sem os pesos, roda em **modo degradado** (só evidência medida) e diz isso na tela. Documentar no README como rodar (`streamlit run painel.py`), a nota sobre o tarball truncado do Zenodo e o contorno com `ecfp4_fold_4.pt`.

- [ ] **Step 4: Run test to verify it passes** — `pytest tests/ -v` (a suíte inteira). Esperado: tudo passa.

- [ ] **Step 5: Commit**

```bash
git add painel.py mods/evidencias/resumo.py tests/test_resumo.py README.md
git commit -m "Add Streamlit panel for stage-1 evidence collection"
```

---

## Critérios de aceite da etapa (do desenho)

1. `sintomas.txt` com "neuropathic pain" devolve alvos que incluam os canais de cálcio α2δ e os transportadores de monoamina.
2. O conjunto-controle produz evidências em que duloxetina e amitriptilina têm linhas `medido` para o transportador de noradrenalina.
3. Nenhuma linha mistura `pchembl` com probabilidade de modelo.
4. Rodar duas vezes com o mesmo manifesto produz o mesmo arquivo.
