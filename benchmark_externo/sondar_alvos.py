"""Fase 1: quantos ativos pos-2022 cada alvo do painel ChEMBL 28 tem.

Uma requisicao por alvo (limit=1, so para ler page_meta.total_count), entao
custa ~616 requisicoes pequenas. O resultado decide quais alvos valem a coleta
completa.
"""
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

AQUI = Path(__file__).parent
MT28 = AQUI / "mt/trained_models/chembl_28_model"
SAIDA = AQUI / "sondagem_alvos.json"
AGENTE = "druglikeness-benchmark/1.0 (research script)"
ANO_TESTE_MIN = 2022
PCHEMBL_MIN = 6.0


def contar(alvo, ano_min=None):
    p = {"target_chembl_id": alvo, "pchembl_value__gte": PCHEMBL_MIN,
         "only": "activity_id", "limit": 1}
    if ano_min:
        p["document_year__gte"] = ano_min
    url = ("https://www.ebi.ac.uk/chembl/api/data/activity.json?"
           + urllib.parse.urlencode(p))
    req = urllib.request.Request(url, headers={"User-Agent": AGENTE})
    for t in range(5):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read())["page_meta"]["total_count"]
        except Exception:
            if t == 4:
                return None
            time.sleep(2 * (t + 1))


def main():
    alvos = json.load(open(MT28 / "targets_28.json"))
    print(f"{len(alvos)} alvos no painel ChEMBL 28", flush=True)
    feito = {}
    if SAIDA.is_file():
        feito = json.loads(SAIDA.read_text())
        print(f"  {len(feito)} ja sondados (cache)", flush=True)
    for i, a in enumerate(alvos):
        if a in feito:
            continue
        # So a contagem de recentes: o filtro por total era redundante
        # (alvo com >= 30 ativos pos-2022 tem folgadamente >= 100 no total)
        # e dobrava o numero de requisicoes.
        feito[a] = {"recentes": contar(a, ANO_TESTE_MIN)}
        if (i + 1) % 25 == 0:
            SAIDA.write_text(json.dumps(feito))
            print(f"  {i+1}/{len(alvos)}", flush=True)
        time.sleep(0.15)
    SAIDA.write_text(json.dumps(feito))

    bons = {a: v for a, v in feito.items()
            if v.get("recentes") and v["recentes"] >= 30}
    print(f"\nalvos com >= 30 ativos pos-{ANO_TESTE_MIN}: {len(bons)}")
    soma = sum(v["recentes"] for v in bons.values())
    print(f"soma de ativos recentes nesses alvos: {soma:,}")
    for lim in (30, 50, 100, 200, 500):
        q = [a for a, v in feito.items()
             if v.get("recentes") and v["recentes"] >= lim]
        print(f"  >= {lim:>3} recentes: {len(q):>3} alvos")


if __name__ == "__main__":
    main()
