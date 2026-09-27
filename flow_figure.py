import matplotlib; matplotlib.use("Agg")
matplotlib.rcParams["ps.fonttype"] = 42   # embed editable text in EPS
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import os
ROOT = os.environ.get("NHANES_ROOT") or os.path.dirname(os.path.abspath(__file__))
F = os.path.join(ROOT, "figures"); SUBMIT = os.path.join(ROOT, "figures_submission")
os.makedirs(F,exist_ok=True); os.makedirs(SUBMIT,exist_ok=True)
def savefig3(name):
    """Submission-quality export: PNG(600dpi) + TIFF(600dpi) + vector EPS."""
    plt.savefig(f"{F}/{name}.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{SUBMIT}/{name}.png", dpi=600, bbox_inches="tight")
    plt.savefig(f"{SUBMIT}/{name}.tiff", dpi=600, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.savefig(f"{SUBMIT}/{name}.eps", bbox_inches="tight")
fig, ax = plt.subplots(figsize=(8.6, 8.2)); ax.set_xlim(0,10); ax.set_ylim(0,11); ax.axis("off")
def box(y, txt, h=0.72, fc="#eef3fb", ec="#2c5f9e"):
    ax.add_patch(FancyBboxPatch((1.0,y),8.0,h, boxstyle="round,pad=0.06", fc=fc, ec=ec, lw=1.2))
    ax.text(5.0, y+h/2, txt, ha="center", va="center", fontsize=9.2, linespacing=1.35)
def arrow(y1,y2):
    ax.add_patch(FancyArrowPatch((5.0,y1),(5.0,y2), arrowstyle="-|>", mutation_scale=12, lw=1.1, color="#333"))
box(10.0, "NHANES 1999\u20132018 (10 cycles incl. 2003\u20132004)\nblood metals measured + mortality linkage\nn = 85,591")
arrow(10.0, 9.6)
box(8.88, "Excluded: age <20 years\nn = 38,707", h=0.72, fc="#fdeaea", ec="#b03a2e")
arrow(8.88, 8.48)
box(7.76, "Excluded: pregnant women  n = 1,394\nnot eligible for linkage  n = 99", h=0.72, fc="#fdeaea", ec="#b03a2e")
arrow(7.76, 7.36)
box(6.64, "Excluded: missing one or more of\nblood lead / cadmium / mercury  n = 8,640", h=0.72, fc="#fdeaea", ec="#b03a2e")
arrow(6.64, 6.24)
box(5.52, "Excluded: no follow-up time  n = 8\nmissing covariates  n = 3,639", h=0.72, fc="#fdeaea", ec="#b03a2e")
arrow(5.52, 5.12)
box(4.10, "Analytic sample  n = 33,104\nall-cause deaths 4,260  |  cardiovascular deaths 1,332\nmedian follow-up 9.3 years", h=0.95, fc="#e8f6ee", ec="#1e8449")
arrow(4.10, 3.70)
box(2.35, "Sub-analyses\n\u2022 selenium / manganese (2011\u20132018): n = 13,460 (805 / 240 deaths)\n\u2022 marine n-3 sensitivity (2011\u20132014): n = 3,258\n\u2022 further sensitivity and subgroup analyses: Tables S1\u2013S24", h=1.30, fc="#fef9e7", ec="#b9770e")
plt.tight_layout(); savefig3("figS1_flow")
print("ok")
