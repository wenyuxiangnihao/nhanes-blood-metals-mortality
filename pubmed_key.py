#!/usr/bin/env python3
"""定向核验关键文献（作者/标题/年份），输出真实 PMID+DOI"""
import urllib.request, urllib.parse, json, time, csv, os
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT = "/sandbox/workspace/heavymetal/results/references_key.csv"

TARGETS = {
 "lead_mortality_LancetPH":   "Lanphear[Author] AND lead[Title] AND mortality[Title]",
 "cadmium_mortality_NHANES":  "cadmium[Title] AND mortality[Title] AND (NHANES[Title/Abstract] OR national health and nutrition[Title/Abstract])",
 "cadmium_CVD_review":        "cadmium[Title] AND cardiovascular[Title] AND (review[pt] OR mechanism*)",
 "mercury_fish_JAMA":         "fish[Title] AND (mercury[Title] OR contaminant*[Title]) AND (risk[Title] OR health[Title])",
 "selenium_mortality":        "selenium[Title] AND mortality[Title] AND serum[Title/Abstract]",
 "manganese_mortality":       "manganese[Title] AND mortality[Title]",
 "qgcomp_keil_2020":          "quantile[Title] AND g-computation[Title]",
 "bkmr_bobb_2015":            "Bayesian kernel machine regression[Title]",
 "ckdepi2021_inker":          "creatinine[Title] AND cystatin[Title] AND equation[Title] AND 2021[dp]",
 "nhanes_analytic_guidelines":"National Health and Nutrition Examination Survey[Title] AND (analytic guidelines[Title] OR sample design[Title])",
 "metal_mixture_mortality":   "metal[Title] AND mixture[Title] AND mortality[Title]",
 "metals_mortality_recent":   "heavy metal*[Title] AND (all-cause[Title/Abstract] OR cardiovascular[Title/Abstract]) AND mortality[Title] AND 2023:2026[dp]",
 "selenium_nhanes":           "selenium[Title] AND (NHANES[Title/Abstract] OR national health and nutrition[Title/Abstract]) AND mortality[Title/Abstract]",
 "mercury_mortality_meta":    "mercury[Title] AND mortality[Title] AND (meta-analysis[Title/Abstract] OR cohort[Title/Abstract])",
}

def eget(u):
    req = urllib.request.Request(u, headers={"User-Agent":"Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as r: return json.loads(r.read().decode())

rows=[]
for tag,q in TARGETS.items():
    try:
        ids = eget(f"{BASE}/esearch.fcgi?db=pubmed&retmax=3&retmode=json&sort=relevance&term="+urllib.parse.quote(q))["esearchresult"]["idlist"]
        if not ids: print(f"  (none) {tag}"); continue
        s = eget(f"{BASE}/esummary.fcgi?db=pubmed&retmode=json&id="+",".join(ids))["result"]
        for pid in ids:
            it=s.get(pid,{}); doi=""
            for a in it.get("articleids",[]):
                if a.get("idtype")=="doi": doi=a.get("value","")
            rows.append(dict(tag=tag,pmid=pid,year=it.get("pubdate","")[:4],journal=it.get("source",""),title=it.get("title","")[:160],doi=doi))
        time.sleep(0.5)
    except Exception as e: print(f"  ERR {tag}: {e}")

with open(OUT,"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["tag","pmid","year","journal","title","doi"]); w.writeheader(); w.writerows(rows)
for r in rows: print(f"{r['tag']:<26}{r['pmid']:<10}{r['year']} {r['journal'][:20]:<20}{r['title'][:78]}")
print(f"\n{len(rows)} refs -> {OUT}")
