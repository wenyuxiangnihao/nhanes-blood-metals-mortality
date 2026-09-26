#!/usr/bin/env python3
import os
ROOT = os.environ.get("NHANES_ROOT") or os.path.dirname(os.path.abspath(__file__))
import matplotlib; matplotlib.use("Agg")
matplotlib.rcParams["ps.fonttype"] = 42   # embed editable text in EPS
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np, csv, os
R=os.path.join(ROOT, "results"); F=os.path.join(ROOT, "figures")
SUBMIT=os.path.join(ROOT, "figures_submission")
os.makedirs(F,exist_ok=True); os.makedirs(SUBMIT,exist_ok=True)
def savefig3(name):
    """Submission-quality export: PNG(600dpi) + TIFF(600dpi) + vector EPS.
    Legacy 300-dpi PNG preview in figures/ is kept unchanged."""
    plt.savefig(f"{F}/{name}.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{SUBMIT}/{name}.png", dpi=600, bbox_inches="tight")
    plt.savefig(f"{SUBMIT}/{name}.tiff", dpi=600, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.savefig(f"{SUBMIT}/{name}.eps", bbox_inches="tight")
def rc(fn):
    with open(os.path.join(R,fn)) as f: return list(csv.DictReader(f))
def fnum(x):
    try: return float(x)
    except: return np.nan

# ---- Fig1: forest per-SD M3 ----
d=rc("single_metal_cox_M1M2M3_10cyc.csv")
mets=[("LBXBPB","Lead"),("LBXBCD","Cadmium"),("LBXTHG","Mercury")]
fig,axes=plt.subplots(1,2,figsize=(11,3.6),sharey=True)
for ax,ev,ttl in zip(axes,["event_all","event_cvd"],["All-cause mortality","CVD mortality"]):
    ys=np.arange(len(mets))[::-1]
    for y,(mk,nm) in zip(ys,mets):
        rr=[x for x in d if x["metal"]==mk and x["outcome"]==ev and x["model"]=="M3"][0]
        hr,lo,hi=fnum(rr["HR"]),fnum(rr["lo"]),fnum(rr["hi"])
        col="#c0392b" if lo>1 else ("#2471a3" if hi<1 else "#7f8c8d")
        ax.errorbar(hr,y,xerr=[[hr-lo],[hi-hr]],fmt="o",color=col,capsize=4,ms=7)
        ax.text(1.55,y,f"{hr:.2f} ({lo:.2f}\u2013{hi:.2f})",va="center",fontsize=9)
    ax.axvline(1,color="k",ls="--",lw=0.8); ax.set_yticks(ys); ax.set_yticklabels([m[1] for m in mets])
    ax.set_xlim(0.75,2.05); ax.set_xlabel("Hazard ratio per SD (log-concentration)"); ax.set_title(ttl,fontsize=11)
fig.suptitle("Blood metals and mortality (weighted Cox, fully adjusted, 10 cycles incl. 2003\u20132004)",y=1.02,fontsize=12)
plt.tight_layout(); savefig3("fig1_forest"); plt.close()

# ---- Fig2: qgcomp weights 3-metal ----
w=rc("qgcomp_weights_10cyc.csv"); nmap={"LBXBPB":"Lead","LBXBCD":"Cadmium","LBXTHG":"Mercury","LBXBSE":"Selenium","LBXBMN":"Manganese"}
fig,axes=plt.subplots(1,2,figsize=(10,3.4))
for ax,ev,ttl in zip(axes,["event_all","event_cvd"],["All-cause","CVD"]):
    rows=[x for x in w if x["set"]=="3-metal" and x["outcome"]==ev]
    names=[nmap.get(x["exposure"].replace("qc_",""),x["exposure"]) for x in rows]
    pos=[fnum(x["weight_pos"]) if x["weight_pos"]!="NA" else 0 for x in rows]
    neg=[-(fnum(x["weight_neg"])) if x["weight_neg"]!="NA" else 0 for x in rows]
    y=np.arange(len(names))
    ax.barh(y,pos,color="#c0392b",label="Positive weight"); ax.barh(y,neg,color="#2471a3",label="Negative weight")
    ax.set_yticks(y); ax.set_yticklabels(names); ax.axvline(0,color="k",lw=0.8)
    ax.set_xlim(-1.1,1.1); ax.set_xlabel("Weight"); ax.set_title(f"3-metal mixture \u2013 {ttl}",fontsize=11); ax.legend(fontsize=8,loc="lower right")
plt.tight_layout(); savefig3("fig2a_qgcomp_weights_3m"); plt.close()

# ---- Fig3: quartiles ----
qt=rc("single_metal_quartiles_10cyc.csv")
fig,axes=plt.subplots(1,3,figsize=(13,3.6),sharey=True)
for ax,m in zip(axes,["LBXBPB","LBXBCD","LBXTHG"]):
    for ev,col,mk in [("event_all","#c0392b","o"),("event_cvd","#2471a3","s")]:
        rows=[x for x in qt if x["metal"]==m and x["outcome"]==ev and x["level"] in ("Q2","Q3","Q4")]
        rows.sort(key=lambda x:x["level"])
        hr=[1.0]+[fnum(x["HR"]) for x in rows]; lo=[1.0]+[fnum(x["lo"]) for x in rows]; hi=[1.0]+[fnum(x["hi"]) for x in rows]
        xs=np.arange(4)
        ax.errorbar(xs,hr,yerr=[np.array(hr)-np.array(lo),np.array(hi)-np.array(hr)],fmt=mk,color=col,capsize=3,label={"event_all":"All-cause mortality","event_cvd":"CVD mortality"}[ev])
    ax.axhline(1,color="k",ls="--",lw=0.8); ax.set_xticks(range(4)); ax.set_xticklabels(["Q1","Q2","Q3","Q4"])
    ax.set_title({"LBXBPB":"Lead","LBXBCD":"Cadmium","LBXTHG":"Mercury"}[m],fontsize=11)
axes[0].set_ylabel("Hazard ratio"); axes[0].legend(fontsize=8)
fig.suptitle("Quartiles of blood metals and mortality (10 cycles)",y=1.02,fontsize=12)
plt.tight_layout(); savefig3("fig3_quartiles"); plt.close()

# ---- Fig4: RCS Pb/Cd/Hg ----
rcv=rc("rcs_curves_10cyc.csv")
fig,axes=plt.subplots(1,3,figsize=(13,3.6))
for ax,m,nm in zip(axes,["LBXBPB","LBXBCD","LBXTHG"],["Lead","Cadmium","Mercury"]):
    rows=sorted([x for x in rcv if x["metal"]==m],key=lambda x:fnum(x["z"]))
    z=np.array([fnum(x["z"]) for x in rows]);hr=np.array([fnum(x["HR"]) for x in rows]);lo=np.array([fnum(x["lo"]) for x in rows]);hi=np.array([fnum(x["hi"]) for x in rows])
    ax.fill_between(z,lo,hi,color="#c0392b",alpha=0.18); ax.plot(z,hr,color="#c0392b",lw=1.8)
    ax.axhline(1,color="k",ls="--",lw=0.8); ax.axvline(0,color="grey",ls=":",lw=0.8)
    ax.set_title(nm,fontsize=11); ax.set_xlabel("SD of log-concentration")
axes[0].set_ylabel("Hazard ratio (all-cause, ref = median)")
fig.suptitle("Restricted cubic spline dose\u2013response (fully adjusted, 10 cycles)",y=1.02,fontsize=12)
plt.tight_layout(); savefig3("fig4_rcs"); plt.close()

# ---- Fig5 (NEW): RCS Se/Mn ----
rcv=rc("rcs_se_mn_10cyc.csv")
fig,axes=plt.subplots(1,2,figsize=(10,3.8))
for ax,m,nm in zip(axes,["LBXBSE","LBXBMN"],["Selenium","Manganese"]):
    rows=sorted([x for x in rcv if x["metal"]==m],key=lambda x:fnum(x["z"]))
    z=np.array([fnum(x["z"]) for x in rows]);hr=np.array([fnum(x["HR"]) for x in rows]);lo=np.array([fnum(x["lo"]) for x in rows]);hi=np.array([fnum(x["hi"]) for x in rows])
    ax.fill_between(z,lo,hi,color="#1e8449",alpha=0.18); ax.plot(z,hr,color="#1e8449",lw=1.8)
    ax.axhline(1,color="k",ls="--",lw=0.8); ax.axvline(0,color="grey",ls=":",lw=0.8)
    ax.set_title(nm,fontsize=11); ax.set_xlabel("SD of log-concentration")
axes[0].set_ylabel("Hazard ratio (all-cause, ref = median)")
fig.suptitle("Restricted cubic spline: selenium & manganese, all-cause mortality (2011\u20132018)",y=1.02,fontsize=12)
plt.tight_layout(); savefig3("fig5_rcs_se_mn"); plt.close()

# ---- FigS1: flow (updated counts) ----
fig,ax=plt.subplots(figsize=(8.8,8.4)); ax.set_xlim(0,10); ax.set_ylim(0,11); ax.axis("off")
def box(y,txt,h=0.72,fc="#eef3fb",ec="#2c5f9e"):
    ax.add_patch(FancyBboxPatch((1.0,y),8.0,h,boxstyle="round,pad=0.06",fc=fc,ec=ec,lw=1.2)); ax.text(5.0,y+h/2,txt,ha="center",va="center",fontsize=9.2,linespacing=1.35)
def arrow(y1,y2): ax.add_patch(FancyArrowPatch((5.0,y1),(5.0,y2),arrowstyle="-|>",mutation_scale=12,lw=1.1,color="#333"))
box(10.0,"NHANES 1999\u20132018 (10 cycles incl. 2003\u20132004)\nblood metals measured + mortality linkage\nn = 85,591")
arrow(10.0,9.6)
box(8.88,"Excluded: age <20 years\nn = 38,707",fc="#fdeaea",ec="#b03a2e")
arrow(8.88,8.48)
box(7.76,"Excluded: pregnant women n = 1,394\nnot eligible for linkage n = 99",fc="#fdeaea",ec="#b03a2e")
arrow(7.76,7.36)
box(6.64,"Excluded: missing one or more of\nblood lead / cadmium / mercury n = 8,640",fc="#fdeaea",ec="#b03a2e")
arrow(6.64,6.24)
box(5.52,"Excluded: no follow-up time n = 8\nmissing covariates n = 3,639",fc="#fdeaea",ec="#b03a2e")
arrow(5.52,5.12)
box(4.10,"Analytic sample n = 33,104\nall-cause deaths 4,260 | cardiovascular deaths 1,332",h=0.95,fc="#e8f6ee",ec="#1e8449")
arrow(4.10,3.70)
box(2.20,"Sub-analyses\n\u2022 selenium / manganese (2011\u20132018): n = 13,460 (805 / 240 deaths)\n\u2022 BKMR (exploratory, supplementary): random subsample n = 800 (106 deaths)",h=1.30,fc="#fef9e7",ec="#b9770e")
plt.tight_layout(); savefig3("figS1_flow"); plt.close()
print("figures:",sorted(os.listdir(F)))
