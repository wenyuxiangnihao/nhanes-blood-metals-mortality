# -*- coding: utf-8 -*-
"""Assemble a single review package from the submission materials."""
import io, re, os

BASE = "/sandbox/workspace/heavymetal"
OUT = os.path.join(BASE, "REVIEW_PACKAGE.md")

def read(p):
    return io.open(os.path.join(BASE, p), encoding="utf-8").read().rstrip("\n")

def demote(text):
    """Push every markdown heading down one level so the package keeps a clean hierarchy."""
    out = []
    for line in text.split("\n"):
        m = re.match(r"^(#{1,5})(\s)", line)
        if m:
            line = "#" + line
        out.append(line)
    return "\n".join(out)

ms   = demote(read("manuscript.md"))
cl   = demote(read("cover_letter_EnvironmentalResearch.md"))
rev  = demote(read("suggested_reviewers.md"))
strb = demote(read("STROBE_checklist.md"))

figs = [
 ("fig1_forest", "Figure 1", "Hazard ratios (95% CI) for all-cause and cardiovascular mortality per SD of log-transformed blood lead, cadmium and mercury (survey-weighted Cox, fully adjusted)."),
 ("fig2a_qgcomp_weights_3m", "Figure 2", "Component weights of the three-metal mixture from quantile g-computation, shown separately for the positive and negative directions."),
 ("fig3_quartiles", "Figure 3", "Hazard ratios (95% CI) across quartiles of blood lead, cadmium and mercury (reference Q1)."),
 ("fig4_rcs", "Figure 4", "Restricted cubic spline curves for lead, cadmium and mercury (SD of log-concentration) and all-cause mortality."),
 ("fig5_rcs_se_mn", "Figure 5", "Restricted cubic spline curves for blood selenium and manganese and all-cause mortality."),
 ("figS1_flow", "Figure S1", "Flow diagram of participant selection."),
 ("bkmr_overall", "Figure S2", "Exploratory BKMR: overall joint effect of the three-metal mixture (probit scale) with 95% credible bands."),
 ("bkmr_single", "Figure S3", "Exploratory BKMR: univariate exposure-response functions for lead, cadmium and mercury."),
]

ms_txt = read("manuscript.md"); tb_txt = read("tables.md")
NREF = max(int(x) for x in re.findall(r"^(\d+)\. ", ms_txt, re.M))
NSA  = len(re.findall(r"^\*\*Table S\d+", tb_txt, re.M))
_body = "\n".join(l for l in ms_txt[ms_txt.index("## Introduction"):ms_txt.index("## References")].split("\n") if not l.strip().startswith("|"))
NWORDS = len([w for w in re.split(r"\s+", re.sub(r"[#*_`]", " ", _body)) if w.strip()])

PREAMBLE = f"""# REVIEW PACKAGE
## Opposing directions in blood metal mixtures and mortality in US adults: a nationally representative sample with mortality linkage (NHANES 1999-2018)

**Revision 7** - prepared 25 September 2026 for external pre-submission review.
Author: Yuxiang Wen (sole author), Department of Cardiology, The First Affiliated Hospital of Yangtze University, Jingzhou, Hubei, China.
Intended journal: *Environmental Research* (Elsevier).

**Contents of this package**
1. Manuscript (full text: abstract, methods, results, discussion, conclusion, Tables 1-5, figure legends, {NREF} references, declarations; ~{NWORDS} words of body text)
2. Cover letter
3. Suggested reviewers
4. STROBE checklist
5. Figures 1-5 and S1-S3 (separate 600 dpi TIFF and vector EPS files supplied on submission; PNG previews inside this package) with legends
6. Supplementary tables S1-S{NSA} (separate file `tables.md` / `submission/supplementary_tables.docx`)
7. Revision log (changes made since the previous review round)

**Reviewer request.** Please assess scientific soundness, clarity and fit for *Environmental Research*, and flag any errors, over-statements or missing analyses. All numerical results are traceable to the project's analysis outputs (shown in the tables); please check internal consistency between the text and the tables.
"""
parts = [PREAMBLE]
parts.append("\n---\n\n# Part 1. Manuscript\n\n" + ms)
parts.append("\n---\n\n# Part 2. Cover letter\n\n" + cl)
parts.append("\n---\n\n# Part 3. Suggested reviewers\n\n" + rev)
parts.append("\n---\n\n# Part 4. STROBE checklist\n\n" + strb)

fig_md = "\n---\n\n# Part 5. Figures (separate image files) with legends\n\n"
for fn, label, leg in figs:
    fig_md += "**%s** (%s.tif / %s.eps) - %s\n\n" % (label, fn, fn, leg)
parts.append(fig_md)

parts.append("""
---

# Part 6. Revision log (changes since the previous review round)

**Statistical and graphical**

1. **Restricted cubic spline confidence bands.** The variance of the spline contrast was computed incorrectly (the reference-point variance was added instead of the variance of the contrast). At the reference point the hazard ratio is now exactly 1.000 with a zero-width confidence interval, and the upper confidence limits of the curves narrowed substantially (lead 24.7 -> 3.17). Figures 4 and 5 were re-drawn and re-exported at 600 dpi.
2. **Assignment of metals to directions.** A sample-splitting analysis was added (quantile g-computation with sample splitting; directions determined in one random half of the data, partial effects estimated in the other half, 200 splits; new Table S19). Mercury was assigned to the negative direction in 200 of 200 splits and the partial effects matched the full-sample estimates (all-cause 1.21 and 0.87).
3. **Missing covariate data.** Multiple imputation in the full eligible sample (n=36,743) and inverse probability of inclusion weighting in the complete-case sample were added (new Table S20). Neither changed any estimate by more than 4%.
4. **Absolute risks.** Table S6 now reports person-years and absolute rates per 1,000 person-years alongside the hazard ratios.
5. **Flow diagram.** Figure S1 previously stated that BKMR had been skipped; it now describes the analysis that was performed (random subsample of 800 participants).

**Corrections of internal inconsistencies**

6. The period-specific summary of the mercury association now matches Table S12 (in 2015-2018 the mercury estimate was at the null).
7. The note to Table S8 no longer claims that excluding observations at the assay floor left the results unchanged; it now refers to the mercury cardiovascular estimate moving to the null (Table S10).
8. The note to Table 1 no longer states that weighted percentages are lower than the crude percentages, which was incorrect for the fourth quartile.
9. The positive- and negative-direction formulae in the Methods now carry the correct subscripts.
10. Figure 3 legends now use the exposure names rather than internal variable names; the Figure 2 legend states that component weights are normalised within each direction.
11. Duplicate artificial-intelligence declarations were merged; the cover letter uses the first person singular throughout.
12. The delayed-entry description now states precisely what happens (participants are retained but their earlier person-time does not enter the risk sets).

**References**

13. Five methodological references were added (Kamenetsky 2025; Carrico 2015; Gennings 2021 comment; Renzetti 2023; Yoshizawa 2002), bringing the list to 40.
""" + "\n".join([]))
io.open(OUT, "w", encoding="utf-8").write("\n".join(parts) + "\n")
print("written:", OUT)
print("words:", len(re.findall(r"[A-Za-z]+", io.open(OUT, encoding="utf-8").read())))
