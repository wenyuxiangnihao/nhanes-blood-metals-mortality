#!/usr/bin/env python3
"""从 results/*.csv 生成 markdown 表格 (Table 1-4)."""
import csv, os
R = "/sandbox/workspace/heavymetal/results"
OUT = "/sandbox/workspace/heavymetal/tables.md"

def rd(fn):
    with open(os.path.join(R, fn)) as f: return list(csv.DictReader(f))

def f2(x): return f"{float(x):.2f}"
def ci(hr, lo, hi): return f"{float(hr):.2f} ({float(lo):.2f}–{float(hi):.2f})"
def pv(p):
    p = float(p)
    return "<0.001" if p < 0.001 else f"{p:.3f}"

META = {"LBXBPB": "Lead", "LBXBCD": "Cadmium", "LBXTHG": "Mercury", "LBXBSE": "Selenium", "LBXBMN": "Manganese"}
L = []

# ---------------- Table 1 ----------------
t1 = rd("table1_by_leadQ.csv")
label = {"RIDAGEYR":"Age, years","BMXBMI":"BMI, kg/m2","INDFMPIR":"Income-to-poverty ratio",
 "sbp":"Systolic BP, mmHg","dbp":"Diastolic BP, mmHg","egfr":"eGFR, mL/min/1.73m2",
 "tc":"Total cholesterol, mg/dL","hdl":"HDL cholesterol, mg/dL","LBXBPB":"Blood lead, ug/dL",
 "LBXBCD":"Blood cadmium, ug/L","LBXTHG":"Blood mercury, ug/L",
 "RIAGENDR=2(%)":"Female, %","smoke=2(%)":"Current smoker, %","dm=1(%)":"Diabetes, %",
 "htn=1(%)":"Hypertension, %","cvd=1(%)":"Baseline CVD, %",
 "event_all(%)":"All-cause death, %","event_cvd(%)":"CVD death, %"}
L.append("### Table 1. Baseline characteristics by quartile of blood lead (weighted)\n")
L.append("| Characteristic | Q1 | Q2 | Q3 | Q4 |")
L.append("|---|---|---|---|---|")
for r in t1:
    v = r["var"]
    if v in ("event_all(%)","event_cvd(%)"):
        row = [f"{float(r[q]):.1f}" for q in ("Q1","Q2","Q3","Q4")]
    elif v in ("RIAGENDR=2(%)","smoke=2(%)","dm=1(%)","htn=1(%)","cvd=1(%)"):
        row = [f"{100*float(r[q]):.1f}" for q in ("Q1","Q2","Q3","Q4")]
    else:
        row = [f2(r[q]) for q in ("Q1","Q2","Q3","Q4")]
    L.append(f"| {label.get(v,v)} | " + " | ".join(row) + " |")
L.append("")

# ---------------- Table 2 ----------------
sm = rd("single_metal_cox.csv"); qt = rd("single_metal_quartiles.csv")
L.append("### Table 2. Blood metals and mortality (weighted Cox, fully adjusted)\n")
L.append("| Metal | Outcome | Per SD (ln) HR (95% CI) | P | Q4 vs Q1 HR (95% CI) | P-trend |")
L.append("|---|---|---|---|---|---|")
for m in ("Lead","Cadmium","Mercury"):
    key = [k for k,v in META.items() if v==m][0]
    for ev,evl in (("event_all","All-cause"),("event_cvd","CVD")):
        a = [x for x in sm if x["metal"]==m and x["outcome"]==ev][0]
        q4 = [x for x in qt if x["metal"]==key and x["outcome"]==ev and x["level"]=="Q4"][0]
        pt = [x for x in qt if x["metal"]==key and x["outcome"]==ev and x["level"]=="P-trend"][0]
        L.append(f"| {m} | {evl} | {ci(a['HR'],a['lo'],a['hi'])} | {pv(a['p'])} | {ci(q4['HR'],q4['lo'],q4['hi'])} | {pv(pt['p'])} |")
L.append("")

# ---------------- Table 3 ----------------
qz = rd("qgcomp_mixture.csv"); qw = rd("qgcomp_weights.csv")
L.append("### Table 3. Mixture effects (quantile g-computation, weighted Cox)\n")
L.append("| Mixture | Outcome | N | HR per one-quartile increase in all metals (95% CI) | P |")
L.append("|---|---|---|---|---|")
for r in qz:
    evl = "All-cause" if r["outcome"]=="event_all" else "CVD"
    L.append(f"| {r['set']} | {evl} | {int(r['n']):,} | {ci(r['HR'],r['lo'],r['hi'])} | {pv(r['p'])} |")
L.append("")
L.append("Component weights (3-metal mixture):")
for ev, evl in (("event_all", "all-cause"), ("event_cvd", "CVD")):
    pos = [x for x in qw if x["set"]=="3-metal" and x["outcome"]==ev and x["weight_pos"]!="NA"]
    neg = [x for x in qw if x["set"]=="3-metal" and x["outcome"]==ev and x["weight_neg"]!="NA"]
    pstr = ", ".join(f"{META[x['exposure'].replace('qc_','')]} {float(x['weight_pos']):.2f}" for x in sorted(pos, key=lambda y:-float(y["weight_pos"])))
    nstr = ", ".join(f"{META[x['exposure'].replace('qc_','')]} {float(x['weight_neg']):.2f}" for x in neg)
    L.append(f"- {evl}: positive direction — {pstr}; negative direction — {nstr}.")
L.append("")

# ---------------- Table 4 ----------------
se = rd("selenium_manganese.csv"); ti = rd("subgroup_interaction_fdr.csv")
n3 = rd("mercury_n3_sensitivity.csv"); sens = rd("sensitivity.csv"); ph = rd("ph_test.csv")
L.append("### Table 4. Essential elements, sensitivity analyses and interactions\n")
L.append("**4a. Selenium and manganese (2011–2018, n=13,460; 805 all-cause and 240 CVD deaths)**\n")
L.append("| Element | Outcome | Per SD (ln) HR (95% CI) | P for non-linearity |")
L.append("|---|---|---|---|")
for m in ("Selenium","Manganese"):
    key = [k for k,v in META.items() if v==m][0]
    for ev,evl in (("event_all","All-cause"),("event_cvd","CVD")):
        a = [x for x in se if x["metal"]==key and x["outcome"]==ev and x["term"]=="per-SD(ln)"][0]
        b = [x for x in se if x["metal"]==key and x["outcome"]==ev and x["term"]=="p_nonlinear"][0]
        L.append(f"| {m} | {evl} | {ci(a['HR'],a['lo'],a['hi'])} | {pv(b['hi'])} |")
L.append("")
L.append("**4b. Mercury adjusted for marine n-3 (subsample n=3,258) and for selenium**\n")
L.append("| Model | All-cause HR (95% CI) |")
L.append("|---|---|")
hg = [x for x in n3 if x["metal"]=="LBXTHG" and x["outcome"]=="event_all"]
for x in hg:
    L.append(f"| Mercury {x['model']} | {ci(x['HR'],x['lo'],x['hi'])} |")
for x in rd("mercury_se_sensitivity.csv"):
    if x["outcome"]=="event_all":
        L.append(f"| {x['model']} | {ci(x['HR'],x['lo'],x['hi'])} |")
L.append("")
L.append("**4c. Sensitivity analyses (per SD, ln)**\n")
L.append("| Model | Lead | Cadmium | Mercury |")
L.append("|---|---|---|---|")
for model,lab in (("main","Main"),("excl_24m","Excluding deaths <24 months"),("excl_baselineCVD","Excluding baseline CVD")):
    row=[]
    for m in ("LBXBPB","LBXBCD","LBXTHG"):
        x = [y for y in sens if y["outcome"]=="event_all" and y["model"]==model and y["metal"]==m][0]
        row.append(ci(x["HR"],x["lo"],x["hi"]))
    L.append(f"| {lab} (n={int(x['n']):,}) | " + " | ".join(row) + " |")
L.append("")
L.append("**4d. Interaction tests (all-cause mortality; BH-FDR across 9 tests)**\n")
L.append("| Metal | Stratifier | P interaction | P (FDR) |")
L.append("|---|---|---|---|")
for x in sorted(ti, key=lambda y: float(y["p_interaction"])):
    sl = {"RIAGENDR":"Sex","smoke":"Smoking","agegrp":"Age (<60/60+)"}[x["strat"]]
    L.append(f"| {META[x['metal']]} | {sl} | {float(x['p_interaction']):.3f} | {float(x['p_fdr']):.3f} |")
L.append("")
L.append("**4e. Proportional-hazards assumption (P for exposure term)**\n")
for x in ph:
    L.append(f"- {META[x['metal']]}: P = {float(x['p_PH_exposure']):.3f}")
L.append("")

open(OUT, "w").write("\n".join(L))
print("\n".join(L))
