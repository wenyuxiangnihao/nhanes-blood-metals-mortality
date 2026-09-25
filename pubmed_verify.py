#!/usr/bin/env python3
"""PubMed 文献核验：esearch -> esummary，输出真实 PMID/DOI"""
import urllib.request, urllib.parse, json, time, csv, os

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT = "/sandbox/workspace/heavymetal/results/references.csv"

QUERIES = [
 ("cadmium",      "blood cadmium mortality NHANES cohort"),
 ("lead",         "blood lead level mortality cardiovascular adults"),
 ("mercury",      "blood mercury fish cardiovascular mortality"),
 ("selenium",     "serum selenium mortality U-shaped nonlinear"),
 ("manganese",    "blood manganese mortality adults"),
 ("qgcomp",       "quantile g-computation mixture exposures Keil"),
 ("bkmr",         "Bayesian kernel machine regression metal mixtures"),
 ("nhanes_design","NHANES complex survey design weighting analysis"),
 ("ckdepi",       "CKD-EPI 2021 creatinine equation race-free"),
 ("mixture_cvd",  "metal mixtures cardiovascular mortality NHANES"),
 ("cd_mech",      "cadmium cardiovascular disease mechanism endothelial"),
 ("hg_review",    "mercury exposure cardiovascular risk review toxicology"),
]

def eget(url):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode())

rows = []
for tag, q in QUERIES:
    try:
        u = f"{BASE}/esearch.fcgi?db=pubmed&retmax=3&retmode=json&term=" + urllib.parse.quote(q)
        ids = eget(u)["esearchresult"]["idlist"]
        if not ids: continue
        s = eget(f"{BASE}/esummary.fcgi?db=pubmed&retmode=json&id=" + ",".join(ids))["result"]
        for pid in ids:
            it = s.get(pid, {})
            doi = ""
            for aid in it.get("articleids", []):
                if aid.get("idtype") == "doi": doi = aid.get("value", "")
            rows.append(dict(tag=tag, pmid=pid, year=it.get("pubdate","")[:4],
                             journal=it.get("source",""), title=it.get("title","")[:150], doi=doi))
        time.sleep(0.5)
    except Exception as e:
        print(f"  ERR {tag}: {e}")

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["tag","pmid","year","journal","title","doi"])
    w.writeheader(); w.writerows(rows)
for r in rows:
    print(f"{r['tag']:<13} {r['pmid']:<9} {r['year']}  {r['journal'][:22]:<22} {r['title'][:80]}")
print(f"\ntotal {len(rows)} refs -> {OUT}")
