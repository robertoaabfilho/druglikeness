# Fase 3 — Ranqueamento com alvos: plano de implementação

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `mods/ranquear.py` ordena compostos por evidência de alvo, em faixas lexicográficas, com o canal negativo deduzido da varredura dos 7.676 alvos do PP4CKD.

**Architecture:** Três módulos puros em `mods/ranque/` (classificação em faixa, canal negativo, esquema das tabelas) e um CLI que costura, no mesmo molde de `mods/aplicar_portao.py`: toda validação antes de abrir saída, saídas em temporários publicadas juntas, manifesto com MD5/versões/parâmetros. Nenhuma constante nova: `MINIMO_TREINO = 50` vem de `mods/evidencias/resumo.py`, as origens de enriquecimento de `mods/evidencias/entrada.py:204`.

**Tech Stack:** Python 3, RDKit, torch (só para o canal negativo), pytest. Sem rede em nenhum teste.

**Spec:** `docs/superpowers/specs/2026-11-26-fase3-ranqueamento-design.md` (requisitos em `2026-11-25-fase3-alvos-no-ranqueamento-requisitos.md`, incluindo o adendo medido de 2026-11-26)

## Global Constraints

- `pchembl` é a **única** coluna comparável entre linhas de `evidencias.tsv`. `valor` + `unidade` nunca se comparam entre unidades.
- Medido e previsto **nunca** no mesmo número. O canal negativo não multiplica a faixa 1.
- Zero constante aprendida ou arbitrada. `MINIMO_TREINO = 50` de `mods/evidencias/resumo.py`; `ORIGENS_DE_ENRIQUECIMENTO = ("open_targets", "reactome")` de `mods/evidencias/entrada.py:204`.
- Ausência de dado nunca vira penalidade: canal negativo sem modelo é `indeterminado`, **nunca** 0.
- Nada desaparece em silêncio: quem não ranqueia sai em `nao_ranqueados.tsv` com o motivo.
- Códigos de saída: `0` completa, `3` parcial, `1` nenhum composto utilizável, `2` erro (nada publicado).
- As três saídas são publicadas juntas, via `mods/evidencias/atomico.temporario`.
- Nenhum teste usa rede nem exige arquivo fora do repositório (nem checkpoint, nem `target_labels.tsv` reais).

## Review Focus

1. **Composto cujo SMILES o RDKit não lê** — sem fingerprint não há forward: o canal negativo sai `indeterminado` e o composto ranqueia pela evidência que tem. Nunca `p_negativo = 0`. (Tarefa 4)
2. **Alvo afirmativo fora do modelo PP4CKD** — tem evidência medida mas nunca previsão. Cai na faixa 1 normalmente; não pode virar faixa 3. (Tarefa 1)
3. **Empate exato de valor dentro da faixa** — ordem estável por `ligante_inchikey_parental`, para que duas execuções deem o mesmo arquivo. (Tarefa 3)
4. **Composto reprovado no portão que também é faixa 3** — sai **uma vez só** em `nao_ranqueados.tsv`, com os dois motivos na mesma linha. (Tarefa 3)
5. **`--portao` ausente** — ranquear sem a etapa 2 é permitido, e a coluna diz que o portão não foi consultado; nunca `aprovado = sim` por omissão. (Tarefa 3)

---

### Task 1: `faixa.py` — classificar o composto e dar a chave de ordem

**Files:**
- Create: `mods/ranque/__init__.py`, `mods/ranque/faixa.py`
- Test: `tests/test_ranque_faixa.py`

**Interfaces:**
- Consumes: `mods.evidencias.entrada.ORIGENS_DE_ENRIQUECIMENTO`
- Produces:
  - `FAIXAS = ("1a", "1b", "2a", "2b", "3")`
  - `classificar(linhas: list[dict], *, p_negativo: float | None) -> dict` → `{"faixa", "valor", "alvo", "p_afirmativo", "pchembl", "medida_sem_pchembl", "observacao"}`. `linhas` são as linhas de `evidencias.tsv` de **um** composto.
  - `chave(resultado: dict, inchikey: str) -> tuple` — chave lexicográfica decrescente usada no `sort`.

- [ ] **Step 1: Write the failing tests**

```python
def test_medido_em_alvo_do_usuario_e_faixa_1a():
    r = classificar([_med(origem="usuario", pchembl="8.2")], p_negativo=0.9)
    assert r["faixa"] == "1a" and r["valor"] == 8.2

def test_medido_em_alvo_do_enriquecimento_e_faixa_1b():
    for origem in ("open_targets", "reactome"):
        assert classificar([_med(origem=origem, pchembl="7.0")],
                           p_negativo=None)["faixa"] == "1b"

def test_canal_negativo_nao_entra_no_valor_da_faixa_1():
    # pchembl e medido, p_negativo e previsto: multiplicar misturaria os dois
    a = classificar([_med(origem="usuario", pchembl="8.2")], p_negativo=0.0)
    b = classificar([_med(origem="usuario", pchembl="8.2")], p_negativo=0.99)
    assert a["valor"] == b["valor"] == 8.2

def test_previsto_multiplica_pelo_canal_negativo():
    r = classificar([_prev(origem="usuario", valor="0.90")], p_negativo=0.5)
    assert r["faixa"] == "2a" and r["valor"] == pytest.approx(0.45)

def test_previsto_sem_canal_negativo_ordena_pelo_p_sozinho():
    r = classificar([_prev(origem="usuario", valor="0.90")], p_negativo=None)
    assert r["valor"] == pytest.approx(0.90)

def test_sem_evidencia_em_alvo_afirmativo_e_faixa_3_sem_valor():
    r = classificar([], p_negativo=0.1)
    assert r["faixa"] == "3" and r["valor"] is None

def test_agrega_varios_alvos_por_maximo_nao_por_soma():
    r = classificar([_med(origem="usuario", pchembl="6.0"),
                     _med(origem="usuario", pchembl="7.5")], p_negativo=None)
    assert r["valor"] == 7.5

def test_medida_sem_pchembl_nao_ordena_e_nao_rebaixa_para_faixa_3():
    # 1 das 49 linhas medidas do evidencias_controle.tsv e assim
    r = classificar([_med(origem="usuario", pchembl=""),
                     _prev(origem="usuario", valor="0.80")], p_negativo=None)
    assert r["faixa"] == "2a"
    assert r["medida_sem_pchembl"] == "sim"
    assert r["valor"] == pytest.approx(0.80)

def test_medida_sem_pchembl_sozinha_nao_e_faixa_1():
    r = classificar([_med(origem="usuario", pchembl="")], p_negativo=None)
    assert r["faixa"] == "3" and r["medida_sem_pchembl"] == "sim"

def test_outra_linha_medida_com_pchembl_mantem_na_faixa_1():
    r = classificar([_med(origem="usuario", pchembl=""),
                     _med(origem="usuario", pchembl="9.1")], p_negativo=None)
    assert r["faixa"] == "1a" and r["valor"] == 9.1

def test_alvo_afirmativo_fora_do_modelo_continua_faixa_1():
    # sem previsao nenhuma: medido basta (Review Focus 2)
    r = classificar([_med(origem="usuario", pchembl="8.0", no_modelo="nao")],
                    p_negativo=None)
    assert r["faixa"] == "1a"

def test_chave_ordena_as_faixas_na_ordem_do_desenho():
    ordem = [chave(classificar(l, p_negativo=None), "K") for l in _uma_de_cada()]
    assert ordem == sorted(ordem, reverse=True)
```

`_med`/`_prev` montam uma linha de `evidencias.tsv` com `tipo_evidencia` `medido`/`previsto`; `_uma_de_cada` devolve as listas de evidência de 1a, 1b, 2a, 2b, 3 nessa ordem.

- [ ] **Step 2: Run to verify they fail**

Run: `python -m pytest tests/test_ranque_faixa.py -v`
Expected: FAIL — `ModuleNotFoundError: mods.ranque`

- [ ] **Step 3: Implement `classificar` e `chave` em `mods/ranque/faixa.py`**

`classificar` separa as linhas por `tipo_evidencia` e por `origem_alvo` ∈ `("usuario",)` vs `ORIGENS_DE_ENRIQUECIMENTO`, e escolhe a primeira faixa não vazia na ordem `1a, 1b, 2a, 2b`. Uma linha medida só conta para a faixa 1 se `pchembl` converter para float. Faixa 3 sai com `valor = None`.

`chave` devolve `(-FAIXAS.index(faixa), valor ou -inf, inchikey)` com o `inchikey` invertido para empate estável decrescente — a docstring diz por que o `inchikey` entra: duas execuções têm de dar o mesmo arquivo.

- [ ] **Step 4: Run to verify they pass**

Run: `python -m pytest tests/test_ranque_faixa.py -v`
Expected: PASS (12 testes)

- [ ] **Step 5: Commit**

```bash
git add mods/ranque/ tests/test_ranque_faixa.py
git commit -m "Fase 3: faixa de evidencia e chave de ordem"
```

---

### Task 2: `negativo.py` — o canal negativo deduzido

**Files:**
- Create: `mods/ranque/negativo.py`
- Create: `tests/dados/rotulos_minimo.tsv` (rótulos sintéticos: humano, não humano com 1.000 ativos, não humano com 12 ativos, ortólogo de nome igual a um alvo pedido)
- Test: `tests/test_ranque_negativo.py`

**Interfaces:**
- Consumes: `mods.evidencias.resumo.MINIMO_TREINO`, `mods.prever_alvos.carregar_rotulos`
- Produces:
  - `mascara(rotulos: list[dict], idx_afirmativos: list[int], *, minimo_treino: int = MINIMO_TREINO) -> numpy.ndarray` (bool)
  - `agregar(probabilidades, mascara) -> float | None` — `max` sobre a máscara; `None` quando a máscara é vazia

- [ ] **Step 1: Write the failing tests**

```python
def test_alvo_humano_nunca_entra_no_canal():
    assert not mascara(ROTULOS, [])[IDX_HUMANO]

def test_alvo_pedido_nao_entra_no_canal():
    assert not mascara(ROTULOS, [IDX_NAO_HUMANO_FORTE])[IDX_NAO_HUMANO_FORTE]

def test_alvo_pouco_treinado_nao_acusa():
    # requisito 1: esparso e ausencia de evidencia, nao evidencia de atividade
    assert not mascara(ROTULOS, [])[IDX_NAO_HUMANO_FRACO]

def test_ortologo_de_alvo_pedido_fica_fora():
    # CHEMBL304 (NET de rato) e ortologo do CHEMBL222 pedido
    assert not mascara(ROTULOS, [IDX_NET_HUMANO])[IDX_NET_RATO]

def test_nao_humano_treinado_e_nao_pedido_entra():
    assert mascara(ROTULOS, [])[IDX_NAO_HUMANO_FORTE]

def test_agregar_e_maximo():
    assert agregar(np.array([0.1, 0.9, 0.4]), np.array([True, True, False])) == 0.9

def test_mascara_vazia_devolve_none_nao_zero():
    # zero e "nenhuma atividade indesejada", que e afirmacao, nao ausencia
    assert agregar(np.array([0.9]), np.array([False])) is None
```

- [ ] **Step 2: Run to verify they fail**

Run: `python -m pytest tests/test_ranque_negativo.py -v`
Expected: FAIL — `cannot import name 'mascara'`

- [ ] **Step 3: Implement `mascara` e `agregar` em `mods/ranque/negativo.py`**

`mascara` = não afirmativo & `organismo != "Homo sapiens"` & `n_ativos >= minimo_treino` & `pref_name` diferente (casefold) de todo `pref_name` dos afirmativos. A docstring registra o número medido: 1.169 alvos no canal para os 9 alvos do controle, e que a exclusão de ortólogo custa 9 alvos sem mudar a AUC.

- [ ] **Step 4: Run to verify they pass**

Run: `python -m pytest tests/test_ranque_negativo.py -v`
Expected: PASS (7 testes)

- [ ] **Step 5: Commit**

```bash
git add mods/ranque/negativo.py tests/test_ranque_negativo.py tests/dados/rotulos_minimo.tsv
git commit -m "Fase 3: canal negativo -- nao pedido, nao humano, treinado, nao ortologo"
```

---

### Task 3: o CLI, as faixas e as duas tabelas (sem canal negativo)

**Files:**
- Create: `mods/ranquear.py`, `mods/ranque/esquema.py`
- Test: `tests/test_ranquear_cli.py`

**Interfaces:**
- Consumes: `mods.ranque.faixa`, `mods.portao.esquema.EscritorTSV`/`ler_tsv`, `mods.evidencias.atomico.temporario`, `mods.evidencias.procedencia`
- Produces:
  - `mods/ranque/esquema.py`: `COLUNAS_RANKING`, `COLUNAS_NAO_RANQUEADOS`
  - `mods/ranquear.py`: `main(argv=None) -> int`; nomes `NOME_RANKING = "ranking.tsv"`, `NOME_NAO_RANQUEADOS = "nao_ranqueados.tsv"`, `NOME_MANIFESTO = "manifesto.json"`

Flags desta tarefa: `--evidencias` (obrigatória), `--saida` (diretório, obrigatória), `--portao`, `--alvos`, `--limite`, `--sem-canal-negativo`.

- [ ] **Step 1: Write the failing tests**

```python
def test_ordena_por_faixa_e_publica_as_tres_saidas(tmp_path):
    código = main(["--evidencias", EVID, "--saida", str(tmp_path),
                   "--sem-canal-negativo"])
    assert código == 0
    assert {p.name for p in tmp_path.iterdir()} >= {
        "ranking.tsv", "nao_ranqueados.tsv", "manifesto.json"}
    assert [l["faixa"] for l in ler(tmp_path / "ranking.tsv")] == \
        sorted([l["faixa"] for l in ler(tmp_path / "ranking.tsv")])

def test_faixa_3_sai_em_nao_ranqueados_sem_posicao(tmp_path):
    ...
    fora = ler(tmp_path / "nao_ranqueados.tsv")
    assert any(l["motivo"] == "sem evidencia em alvo afirmativo" for l in fora)
    assert "posicao" not in fora[0]

def test_reprovado_no_portao_nao_ranqueia_mas_sai_com_o_criterio(tmp_path):
    ...
    assert l["motivo"].startswith("reprovado no portao") and "alerta_bula" in l["motivo"]

def test_reprovado_e_faixa_3_aparece_uma_vez_so_com_os_dois_motivos(tmp_path):
    # Review Focus 4
    nomes = [l["#ligante"] for l in ler(tmp_path / "nao_ranqueados.tsv")]
    assert len(nomes) == len(set(nomes))

def test_marcado_no_portao_ranqueia_carregando_as_marcas(tmp_path):
    # a etapa 2 marca porque reprovar apaga os positivos; o ranking nao desfaz
    assert l["portao_aprovado"] == "sim" and l["portao_motivos"] == ""

def test_sem_portao_a_coluna_diz_que_nao_foi_consultado(tmp_path):
    # Review Focus 5: nunca `aprovado = sim` por omissao
    assert ler(tmp_path / "ranking.tsv")[0]["portao_aprovado"] == "nao consultado"

def test_empate_de_valor_ordena_estavel_pelo_inchikey(tmp_path):
    # Review Focus 3: duas execucoes dao o mesmo arquivo
    assert primeira_execucao == segunda_execucao

def test_saida_que_aponta_para_arquivo_existente_e_erro_antes_de_escrever(tmp_path):
    assert main(["--evidencias", EVID, "--saida", str(arquivo)]) == 2

def test_nenhum_composto_utilizavel_sai_1_e_nao_escreve_nada(tmp_path):
    assert main(["--evidencias", VAZIO, "--saida", str(tmp_path)]) == 1
    assert not list(tmp_path.iterdir())
```

- [ ] **Step 2: Run to verify they fail**

Run: `python -m pytest tests/test_ranquear_cli.py -v`
Expected: FAIL — `No module named 'mods.ranquear'`

- [ ] **Step 3: Implement `mods/ranque/esquema.py` e `mods/ranquear.py`**

`COLUNAS_RANKING` = `posicao`, `#ligante`, `ligante_inchikey_parental`, `faixa`, `valor`, `alvo`, `p_afirmativo`, `pchembl`, `p_negativo`, `alvo_negativo`, `medida_sem_pchembl`, `portao_aprovado`, `portao_motivos`, `procedencia`. `COLUNAS_NAO_RANQUEADOS` = `#ligante`, `ligante_inchikey_parental`, `motivo`, `faixa`, `portao_aprovado`, `portao_motivos`, `procedencia`.

O CLI agrupa `evidencias.tsv` por `ligante_inchikey_parental`, chama `faixa.classificar` com `p_negativo=None`, ordena por `faixa.chave`, e escreve. Um composto reprovado no portão **e** de faixa 3 produz uma linha só, com os dois motivos unidos por `;`. Toda validação de flag acontece antes de abrir qualquer saída (o padrão de `aplicar_portao.py`, que existe porque a validação rodava depois de carregar dez modelos).

- [ ] **Step 4: Run to verify they pass**

Run: `python -m pytest tests/test_ranquear_cli.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add mods/ranquear.py mods/ranque/esquema.py tests/test_ranquear_cli.py
git commit -m "Fase 3: CLI, faixas e as duas tabelas"
```

---

### Task 4: o canal negativo ligado, com o guarda de checkpoint

**Files:**
- Modify: `mods/ranquear.py`
- Modify: `tests/test_ranquear_cli.py`

**Interfaces:**
- Consumes: `mods.ranque.negativo.mascara`/`agregar`, `mods.prever_alvos.carregar_modelo`/`carregar_rotulos`/`fingerprint_de`, `mods.evidencias.manifesto.ler_manifesto`
- Produces: flags `--checkpoint`, `--rotulos`, `--manifesto-etapa1`, `--lote`, `--gpu`

- [ ] **Step 1: Write the failing tests**

```python
def test_recusa_checkpoint_diferente_do_que_a_etapa_1_registrou(tmp_path):
    # os dois lados de `p x (1 - p_negativo)` tem de vir do MESMO modelo
    código = main([..., "--checkpoint", OUTRO, "--manifesto-etapa1", MAN])
    assert código == 2
    assert "md5" in capsys.readouterr().err.lower()

def test_sem_pesos_registra_canal_negativo_falso_e_nao_usa_zero(tmp_path):
    main(["--evidencias", EVID, "--saida", str(tmp_path)])   # sem --checkpoint
    assert json.load(...)["parametros"]["canal_negativo"] is False
    assert ler(tmp_path / "ranking.tsv")[0]["p_negativo"] == "indeterminado"

def test_smiles_ilegivel_da_canal_indeterminado_e_ranqueia_pelo_resto(tmp_path):
    # Review Focus 1
    l = linha_de(tmp_path, "ilegivel")
    assert l["p_negativo"] == "indeterminado" and l["posicao"]

def test_canal_negativo_multiplica_so_a_faixa_prevista(tmp_path):
    assert faixa_1["valor"] == pytest.approx(PCHEMBL_ORIGINAL)

def test_dispositivo_registrado_e_o_que_rodou_de_fato(tmp_path):
    assert json.load(...)["parametros"]["dispositivo"] in ("cpu", "cuda")
```

- [ ] **Step 2: Run to verify they fail**

Run: `python -m pytest tests/test_ranquear_cli.py -k canal -v`
Expected: FAIL — flag `--checkpoint` inexistente

- [ ] **Step 3: Ligar o canal em `mods/ranquear.py`**

Antes do forward, comparar o MD5 do `--checkpoint` com o que `--manifesto-etapa1` registra e **recusar com código 2** nomeando os dois, quando diferirem. O forward roda em lotes de `--lote` compostos, guardando por composto apenas `agregar(probabilidades, mascara)` e o nome do alvo do argmax — nunca a matriz. Composto sem fingerprint e execução sem modelo dão `p_negativo = "indeterminado"`; o padrão é CPU, como na etapa 2, porque a GPU não reproduz.

- [ ] **Step 4: Run to verify they pass**

Run: `python -m pytest tests/test_ranquear_cli.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add mods/ranquear.py tests/test_ranquear_cli.py
git commit -m "Fase 3: canal negativo ligado, com guarda de checkpoint"
```

---

### Task 5: as duas AUC do adendo, como teste

**Files:**
- Create: `tests/dados/controle_previsoes.tsv` (por fármaco do controle: `p` dos 9 alvos pedidos e `p_negativo`, gravados do `ecfp4_fold_4.pt` real)
- Modify: `tests/dados/LICENSE` (procedência: saída do PP4CKD sobre o controle; o que a fixture fixa é a **aritmética**, não o modelo)
- Test: `tests/test_ranque_auc.py`

**Interfaces:**
- Consumes: `mods.ranque.faixa.classificar`/`chave`

- [ ] **Step 1: Write the failing test**

```python
def test_reproduz_as_duas_auc_do_adendo():
    assert auc(ordenar(previsoes, canal=False)) == pytest.approx(0.750, abs=0.001)
    assert auc(ordenar(previsoes, canal=True)) == pytest.approx(0.806, abs=0.001)

def test_contar_promiscuidade_seria_pior_e_a_amitriptilina_zeraria():
    # Desenho DESCARTADO, fixado aqui para nao voltar. A contagem nao existe em
    # `mods/`: o teste a calcula da coluna `n_promiscuo` da fixture, para que
    # nenhum codigo de producao precise suportar um desenho rejeitado.
    escores = {p["nome"]: p["p_afirmativo"] * (1 - p["n_promiscuo"] / PIOR)
               for p in previsoes}
    assert auc(escores) == pytest.approx(0.583, abs=0.001)
    assert escores["amitriptyline"] == 0.0
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/test_ranque_auc.py -v`
Expected: FAIL — fixture inexistente

- [ ] **Step 3: Gravar a fixture e implementar o teste**

O gravador vive no scratchpad, não no repositório, e o comentário do arquivo diz qual checkpoint e quais 9 alvos. A função `auc` é a contagem de pares concordantes com empate valendo meio — não entra em `mods/`, é do teste.

- [ ] **Step 4: Run to verify it passes**

Run: `python -m pytest tests/test_ranque_auc.py -v`
Expected: PASS, com 0,750 e 0,806

- [ ] **Step 5: Commit**

```bash
git add tests/test_ranque_auc.py tests/dados/controle_previsoes.tsv tests/dados/LICENSE
git commit -m "Fase 3: as duas AUC do adendo viram teste"
```

---

### Task 6: manifesto e repetição

**Files:**
- Modify: `mods/ranquear.py`
- Modify: `tests/test_ranquear_cli.py`

**Interfaces:**
- Consumes: `mods.evidencias.manifesto.texto_do_manifesto`/`verificar_manifesto`/`hashes_de_entrada_ausentes`
- Produces: flag `--manifesto`; `VERSOES_DA_REPETICAO`, `PARAMETROS_DA_REPETICAO`

- [ ] **Step 1: Write the failing tests**

```python
def test_repeticao_da_o_mesmo_arquivo_byte_a_byte(tmp_path):
    main([..., "--saida", str(a)])
    main([..., "--saida", str(b), "--manifesto", str(a / "manifesto.json")])
    assert (a / "ranking.tsv").read_bytes() == (b / "ranking.tsv").read_bytes()

def test_repeticao_recusa_se_o_md5_das_evidencias_mudou(tmp_path):
    (tmp_path / "evidencias.tsv").write_text(OUTRO_CONTEUDO)
    assert main([..., "--manifesto", MAN]) == 2
    assert "evidencias_md5" in capsys.readouterr().err

def test_repeticao_recusa_se_o_minimo_treino_mudou(tmp_path):
    assert main([..., "--manifesto", MAN, "--minimo-treino", "10"]) == 2

def test_repeticao_recusa_se_o_canal_negativo_foi_desligado(tmp_path):
    assert main([..., "--manifesto", MAN, "--sem-canal-negativo"]) == 2

def test_manifesto_registra_canal_negativo_minimo_treino_e_dispositivo(tmp_path):
    p = json.loads((tmp_path / "manifesto.json").read_text())["parametros"]
    assert set(p) >= {"canal_negativo", "minimo_treino", "dispositivo", "lote"}

def test_manifesto_ausente_no_destino_nunca_fica_sem_tabela(tmp_path):
    # tabela sem manifesto nao e citavel: as tres saem juntas ou nenhuma sai
    assert not list(tmp_path.iterdir()) or (tmp_path / "manifesto.json").exists()
```

- [ ] **Step 2: Run to verify they fail**

Run: `python -m pytest tests/test_ranquear_cli.py -k repeticao -v`
Expected: FAIL — flag `--manifesto` inexistente

- [ ] **Step 3: Implement o modo de repetição**

`VERSOES_DA_REPETICAO` = `("rdkit", "pp4ckd_md5", "evidencias_md5", "portao_md5", "alvos_md5")`; `PARAMETROS_DA_REPETICAO` = `("canal_negativo", "minimo_treino", "limite", "dispositivo", "lote")`. O manifesto novo sai como `manifesto_repeticao.json` e nunca sobrescreve o original, como na etapa 2.

- [ ] **Step 4: Run to verify they pass**

Run: `python -m pytest tests/ -q`
Expected: toda a suíte verde

- [ ] **Step 5: Commit**

```bash
git add mods/ranquear.py tests/test_ranquear_cli.py
git commit -m "Fase 3: manifesto e repeticao"
```

---

### Task 7: README

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Escrever a seção 13 — Fase 3**

O que ela ordena e o que **não** ordena (faixa 3), a tabela das cinco faixas com o motivo medido de cada decisão (por que o canal negativo não multiplica a faixa 1; por que a procedência é sub-faixa e não peso; por que `max` e não soma), as duas AUC (0,750 e 0,806) com a **ressalva de n = 12 e dos seis antibióticos**, as três saídas e suas colunas, e o comando completo.

- [ ] **Step 2: Marcar `ranking_multibase.py` como substituído**

Na seção 2 e na tabela das 6 etapas, com o número que decide: os dois QED mais altos do controle são ciprofloxacino 0,893 e linezolida 0,892, ambos do grupo que **causa**.

- [ ] **Step 3: Atualizar "Limitações conhecidas"**

Trocar a linha "a fase 3 não existe ainda" pelo que fica em aberto: faixa permeável não medida, canal negativo validado só em n = 12, e o lado humano do canal medido como prejudicial (AUC 0,389 sozinho e invertido).

- [ ] **Step 4: Commit**

```bash
git add README.md
git commit -m "README: documenta a fase 3 e marca o ranqueamento antigo como substituido"
```

---

## Ordem e verificação

Tarefas 1→6 em sequência (cada uma deixa a suíte verde), 7 no fim. Verificação final: suíte completa verde, as duas AUC reproduzidas, e uma execução real sobre o controle de 12 com o `evidencias_controle.tsv` e o `portao.tsv` da etapa 2.

## Fora de escopo

Remover `ranking_multibase.py`, `farmaco_completo.py` e `checar_alertas_adversos.py` — saem quando esta fase estiver validada contra conjunto independente. Faixa permeável. Aprender pesos.
