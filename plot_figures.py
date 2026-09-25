#!/usr/bin/env python3
"""生成图表：森林图 / qgcomp权重 / 四分位剂量反应"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, csv, os

R = "/sandbox/workspace/heavymetal/results"
F = "/sandbox/workspace/heavymetal/figures"
os.makedirs(F, exist_ok=True)

def readcsv(fn):
    with open(os.path.join(R, fn)) as f:
        return list(csv.DictReader(f))

# ---------- Fig1: forest plot (single metal per-SD) ----------
sm = readcsv("single_metal_cox.csv")
mets = ["Lead", "Cadmium", "Mercury"]
fig, axes = plt.subplots(1, 2, figsize=(11, 3.6), sharey=True)
for ax, ev, ttl in zip(axes, ["event_all", "event_cvd"], ["All-cause mortality", "CVD mortality"]):
    ys = np.arange(len(mets))[::-1]
    for y, m in zip(ys, mets):
        r = [x for x in sm if x["metal"] == m and x["outcome"] == ev][0]
        hr, lo, hi = float(r["HR"]), float(r["lo"]), float(r["hi"])
        col = "#c0392b" if lo > 1 else ("#2471a3" if hi < 1 else "#7f8c8d")
        ax.errorbar(hr, y, xerr=[[hr-lo], [hi-hr]], fmt="o", color=col, capsize=4, ms=7)
        ax.text(1.62, y, f"{hr:.2f} ({lo:.2f}–{hi:.2f})", va="center", fontsize=9)
    ax.axvline(1, color="k", ls="--", lw=0.8)
    ax.set_yticks(ys); ax.set_yticklabels(mets)
    ax.set_xlim(0.75, 2.1); ax.set_xlabel("Hazard ratio per SD (log-concentration)")
    ax.set_title(ttl, fontsize=11)
    ax.axvspan(1.5, 2.1, color="white")
fig.suptitle("Blood metals and mortality (weighted Cox, fully adjusted)", y=1.02, fontsize=12)
plt.tight_layout(); plt.savefig(f"{F}/fig1_forest.png", dpi=300, bbox_inches="tight"); plt.close()

# ---------- Fig2: qgcomp weights ----------
w = readcsv("qgcomp_weights.csv")
for st, fname in [("3-metal", "fig2a_qgcomp_weights_3m.png")]:
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.4))
    for ax, ev, ttl in zip(axes, ["event_all", "event_cvd"], ["All-cause", "CVD"]):
        rows = [x for x in w if x["set"] == st and x["outcome"] == ev]
        names = [{"LBXBPB":"Lead","LBXBCD":"Cadmium","LBXTHG":"Mercury",
                  "LBXBSE":"Selenium","LBXBMN":"Manganese"}.get(x["exposure"].replace("qc_",""), x["exposure"]) for x in rows]
        pos = [float(x["weight_pos"]) if x["weight_pos"] != "NA" else 0 for x in rows]
        neg = [-(float(x["weight_neg"])) if x["weight_neg"] != "NA" else 0 for x in rows]
        y = np.arange(len(names))
        ax.barh(y, pos, color="#c0392b", label="Positive weight")
        ax.barh(y, neg, color="#2471a3", label="Negative weight")
        ax.set_yticks(y); ax.set_yticklabels(names)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_xlim(-1.1, 1.1); ax.set_xlabel("Weight"); ax.set_title(f"{st} mixture – {ttl}", fontsize=11)
        ax.legend(fontsize=8, loc="lower right")
    plt.tight_layout(); plt.savefig(f"{F}/{fname}", dpi=300, bbox_inches="tight"); plt.close()

# ---------- Fig3: quartile dose-response ----------
qt = readcsv("single_metal_quartiles.csv")
fig, axes = plt.subplots(1, 3, figsize=(13, 3.6), sharey=True)
for ax, m in zip(axes, ["LBXBPB", "LBXBCD", "LBXTHG"]):
    for ev, col, mk in [("event_all", "#c0392b", "o"), ("event_cvd", "#2471a3", "s")]:
        rows = [x for x in qt if x["metal"] == m and x["outcome"] == ev and x["level"].startswith("Q")]
        rows.sort(key=lambda x: x["level"])
        hr = [1.0]+[float(x["HR"]) for x in rows]; lo=[1.0]+[float(x["lo"]) for x in rows]; hi=[1.0]+[float(x["hi"]) for x in rows]
        xs = np.arange(4)
        ax.errorbar(xs, hr, yerr=[np.array(hr)-np.array(lo), np.array(hi)-np.array(hr)], fmt=mk, color=col, capsize=3, label=ev)
    ax.axhline(1, color="k", ls="--", lw=0.8); ax.set_xticks(range(4)); ax.set_xticklabels(["Q1","Q2","Q3","Q4"])
    ax.set_title({"LBXBPB":"Lead","LBXBCD":"Cadmium","LBXTHG":"Mercury"}[m], fontsize=11)
axes[0].set_ylabel("Hazard ratio"); axes[0].legend(fontsize=8)
fig.suptitle("Quartiles of blood metals and mortality", y=1.02, fontsize=12)
plt.tight_layout(); plt.savefig(f"{F}/fig3_quartiles.png", dpi=300, bbox_inches="tight"); plt.close()

print("figures written:"); [print(" ", x) for x in sorted(os.listdir(F))]

# ---------- Fig4: RCS dose-response curves ----------
try:
    rc = readcsv("rcs_curves.csv")
    if rc:
        fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
        for ax, m, nm in zip(axes, ["LBXBPB","LBXBCD","LBXTHG"], ["Lead","Cadmium","Mercury"]):
            rows = sorted([x for x in rc if x["metal"]==m], key=lambda x: float(x["z"]))
            z = np.array([float(x["z"]) for x in rows]); hr=np.array([float(x["HR"]) for x in rows])
            lo=np.array([float(x["lo"]) for x in rows]); hi=np.array([float(x["hi"]) for x in rows])
            ax.fill_between(z, lo, hi, color="#c0392b", alpha=0.18)
            ax.plot(z, hr, color="#c0392b", lw=1.8)
            ax.axhline(1, color="k", ls="--", lw=0.8); ax.axvline(0, color="grey", ls=":", lw=0.8)
            ax.set_title(nm, fontsize=11); ax.set_xlabel("SD of log-concentration")
        axes[0].set_ylabel("Hazard ratio (all-cause, ref = median)")
        fig.suptitle("Restricted cubic spline dose–response (fully adjusted)", y=1.02, fontsize=12)
        plt.tight_layout(); plt.savefig(f"{F}/fig4_rcs.png", dpi=300, bbox_inches="tight"); plt.close()
        print("  fig4_rcs.png")
except Exception as e:
    print("fig4 error:", e)

