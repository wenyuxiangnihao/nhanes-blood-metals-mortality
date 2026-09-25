# -*- coding: utf-8 -*-
"""Apply the merged (round-2 + round-3) review edits to manuscript.md."""
import io, re
P = '/sandbox/workspace/heavymetal/manuscript.md'
s = io.open(P, encoding='utf-8').read()
def rep(o, n, tag):
    global s
    if s.count(o) != 1:
        raise SystemExit('MATCH FAIL (%d) [%s]: %s' % (s.count(o), tag, o[:100]))
    s = s.replace(o, n)
    print('ok', tag)

# ---------- Abstract ----------
i = s.index('**Background.**'); j = s.index('## Introduction')
ABS = '''**Background.** Blood metals are associated with mortality, but the joint (mixture) effect is usually summarized as a single number and opposing directions are rarely separated.

**Methods.** We analyzed 33,104 adults aged \u226520 years from ten NHANES cycles (1999\u20132018) with blood lead, cadmium and mercury and mortality follow-up to 2019; selenium and manganese were available for 2011\u20132018 (n=13,460). Survey-weighted Cox models, quantile g-computation with directional decomposition, and exploratory BKMR were used. Individual-metal and three-metal mixture analyses were primary; per-cause, essential-element and ratio analyses were secondary.

**Results.** There were 4,260 all-cause and 1,332 cardiovascular deaths (median follow-up 9.3 years). Per SD of log-concentration (fully adjusted), cadmium (HR 1.25, 95% CI 1.18\u20131.33) and lead (1.08, 1.02\u20131.13) were positively associated with all-cause mortality and mercury inversely (0.86, 0.82\u20130.89); for cardiovascular mortality cadmium (1.19, 1.08\u20131.31), lead (1.18, 1.10\u20131.27) and mercury (0.90, 0.83\u20130.98) were associated. The net mixture effect was null for all-cause mortality (1.05, 0.98\u20131.13) and borderline for cardiovascular mortality (1.13, 1.00\u20131.28; P=0.049): the positive-direction partial effect (cadmium, lead) was 1.21 (1.14\u20131.28) and the negative-direction effect (mercury) 0.87 (0.84\u20130.90). Mercury was inversely associated with cancer (0.87, 0.79\u20130.95) and cerebrovascular deaths (0.79, 0.67\u20130.93) but not significantly with heart disease (0.92, 0.84\u20131.02); the association persisted with entry delayed to 60 months and among never smokers, whereas the lead estimate was attenuated in sensitivity analyses. Cadmium was also associated with chronic lower respiratory disease (2.12, 1.71\u20132.64), but not among never smokers, so residual confounding by smoking dose is likely. Selenium was approximately monotonically inverse (0.81, 0.72\u20130.92) and manganese U-shaped (nadir \u22488.3 \u00b5g/L; P=0.002).

**Conclusions.** The inverse mercury association was confined to cancer and cerebrovascular deaths rather than heart disease \u2014 a pattern that does not fit a simply cardioprotective, fish-related mechanism \u2014 and was not explained by marine n-3 fatty acids or selenium. Because blood mercury largely reflects recent fish intake, reverse causation cannot be excluded and the association should not be read as evidence that mercury is protective.

**Keywords:** blood metals; cadmium; lead; mercury; selenium; manganese; metal mixtures; quantile g-computation; mortality; NHANES

'''
s = s[:i] + ABS + s[j:]
ab = re.sub(r'\*\*[^*]*\*\*', '', ABS[:ABS.index('**Keywords')])
print('abstract words:', len(re.findall(r"[A-Za-z0-9][A-Za-z0-9\u2013./\-]*", ab)))

# ---------- Methods ----------
rep('Blood-metal subsample weights (WTSH2YR) were not present in the public-use demographic files for all cycles (they are unavailable for 2013\u20132016 in the current release), so MEC weights \u2014 the conventional choice for NHANES mortality analyses of the full examination sample \u2014 were used throughout; estimates are therefore representative of the examined sample.',
 'Blood-metal subsample weights (WTSH2YR) are provided only for the 2013\u20132014 and 2015\u20132016 cycles, and in the laboratory files (PBCD_H and PBCD_I) rather than in the demographic files; because blood metals were measured in a random one-half subsample in those two cycles and no corresponding subsample weight exists for the other cycles, the MEC examination weight (WTMEC2YR) was used throughout so that the weighting is consistent across the pooled period, and estimates are therefore representative of the examined sample. A sensitivity analysis substituting WTSH2YR for 2013\u20132016, and a comparison of the two weightings within that subsample, are reported in Table S9.',
 'M:WTSH2YR')

rep('Non-linearity was assessed by adding a quadratic term and, graphically, with restricted cubic splines (4 knots).',
 'Non-linearity was assessed by adding a quadratic term and, graphically, with restricted cubic splines (4 knots); trend across quartiles was tested by entering each quartile\u2019s median concentration as a continuous variable.',
 'M:trend')

rep('Sensitivity analyses: (i) delayed entry with risk sets beginning 24 months after examination (to examine reverse causation without dropping participants); (ii) exclusion of participants with baseline cardiovascular disease; (iii) omission of hypertension and diabetes (possible mediators); (iv) additional adjustment for serum cotinine; (v) additional adjustment for estimated glomerular filtration rate, and exclusion of participants with eGFR<60 mL/min/1.73 m\u00b2; and (vi) additional adjustment for alcohol intake and physical activity.',
 'Sensitivity analyses: (i) delayed entry with risk sets beginning 24, 48 or 60 months after examination (to examine reverse causation without dropping participants; longer delays are pertinent to mercury, whose biological half-life is only about 50\u201370 days, whereas blood cadmium and lead largely reflect cumulative exposure); (ii) exclusion of participants with baseline cardiovascular disease; (iii) omission of hypertension and diabetes (possible mediators); (iv) additional adjustment for serum cotinine; (v) additional adjustment for estimated glomerular filtration rate, and exclusion of participants with eGFR<60 mL/min/1.73 m\u00b2; (vi) additional adjustment for alcohol intake and physical activity; (vii) substitution of the blood-metal subsample weight for 2013\u20132016 and additional adjustment for survey cycle; and (viii) exclusion of observations at the assay floor (the lowest reportable concentration within each cycle).',
 'M:sens-list')

rep('Period-specific estimates were obtained for 1999\u20132010, 2011\u20132014 and 2015\u20132018.',
 'Period-specific estimates were obtained for 1999\u20132010, 2011\u20132014 and 2015\u20132018 (Table S12). Survey strata and primary sampling units are nested within survey cycle, so cycle was not included in the main models; a sensitivity analysis adding cycle as a covariate was performed. The selenium\u2013mercury interaction was tested with a product term and by selenium tertile, because the Se:Hg molar ratio is an algebraic combination of the two marginal exposures rather than an independent test of antagonism. E-values were calculated for the inverse mercury association.',
 'M:period-cycle-se')

# ---------- Results ----------
rep('Mercury was inversely associated with both all-cause (0.86, 0.82\u20130.89) and cardiovascular mortality (0.90, 0.83\u20130.98); the crude to adjusted attenuation was substantial (crude 0.83 and 0.86, age/sex/race-adjusted 0.75 and 0.78), indicating that much of the crude inverse association reflected sociodemographic and lifestyle factors.',
 'Mercury was inversely associated with both all-cause (0.86, 0.82\u20130.89) and cardiovascular mortality (0.90, 0.83\u20130.98). Adjustment for age, sex and race/ethnicity strengthened the inverse association (crude 0.83 and 0.86; 0.75 and 0.78 after adjustment), whereas further adjustment for socioeconomic position, smoking, adiposity and comorbidity attenuated it (0.75\u21920.86 and 0.78\u21920.90). Demographic factors therefore masked rather than produced the inverse association, which is expected because blood mercury rises with age and with fish intake while age is a strong risk factor for death.',
 'R:mercury-crude')

rep('although its association with all-cause mortality attenuated towards the null in some sensitivity analyses.',
 'although its association with all-cause mortality attenuated towards the null in some sensitivity analyses; its quartile pattern was not monotonic (Q2 0.95, Q3 0.89, Q4 1.07), although the trend test based on quartile medians was significant (P=0.009).',
 'R:lead-quartile')

rep('In the exploratory BKMR analysis the overall mixture association was small and non-monotonic, with both extreme tails differing from the median (probit-scale difference \u22120.48 at the 10th and +0.43 at the 90th percentile), consistent with opposing directions.',
 'In the exploratory BKMR analysis (800 participants, 106 deaths; unweighted) the joint effect function increased across percentiles of the mixture (probit-scale difference from the median \u22120.48, 95% credible interval \u22120.75 to \u22120.21, at the 10th percentile and +0.43, 0.14 to 0.72, at the 90th). This monotonic pattern is consistent with cadmium and lead dominating the joint contrast rather than with a protective effect of mercury, but with 106 deaths the credible intervals are wide and the analysis remains exploratory.',
 'R:bkmr')

rep('Mercury was inversely associated with deaths from cancer (HR 0.87, 0.79\u20130.95) and cerebrovascular disease (0.79, 0.67\u20130.93) but not with heart disease (0.92, 0.84\u20131.02). Cadmium was positively associated with all four causes examined, most strongly with chronic lower respiratory disease (2.12, 1.71\u20132.64), with a gradient across smoking status (current smokers 2.27, 1.58\u20133.28; former 2.12, 1.54\u20132.92; never 1.77, 0.97\u20133.23).',
 'Mercury was inversely associated with deaths from cancer (HR 0.87, 0.79\u20130.95) and cerebrovascular disease (0.79, 0.67\u20130.93); for heart disease the estimate was close to the null and not statistically significant (0.92, 0.84\u20131.02). Cadmium was positively associated with all four causes examined, most strongly with chronic lower respiratory disease (2.12, 1.71\u20132.64). That association was not significant among never smokers (1.77, 0.97\u20133.23), and the apparent gradient across smoking status (current 2.27, 1.58\u20133.28; former 2.12, 1.54\u20132.92) is likely to reflect residual confounding by smoking dose rather than genuine effect modification, since blood cadmium is dominated by tobacco smoke and a single serum cotinine measurement (half-life about 16 hours) cannot capture lifetime smoking dose.',
 'R:cause-specific')

rep('Manganese showed a U-shaped association: the linear term was close to the null (all-cause 1.04, 0.94\u20131.15) but the non-linear component was significant (P for non-linearity 0.002 for all-cause and <0.001 for cardiovascular mortality), with a nadir near 8.3 \u00b5g/L.',
 'Manganese showed a U-shaped association: the linear term was close to the null (all-cause 1.04, 0.94\u20131.15) but the non-linear component was significant (P for non-linearity 0.002 for all-cause and <0.001 for cardiovascular mortality), with a nadir near 8.3 \u00b5g/L, close to the inflection point of about 7.1 \u00b5g/L reported independently in NHANES [22]. The inverse selenium association is consistent with the L-shaped inverse relation reported in a recent NHANES analysis [23].',
 'R:emn-refs')

rep('A higher Se:Hg molar ratio was weakly associated with higher all-cause mortality (per SD 1.14, 1.02\u20131.28); however, models adjusting mercury for the molar ratio are not interpretable because the log-ratio is collinear with log-mercury by construction.',
 'A higher Se:Hg molar ratio was weakly associated with higher all-cause mortality (per SD 1.14, 1.02\u20131.28); this ratio is an algebraic reparameterisation of the two marginal associations under a constrained model and does not test selenium\u2013mercury antagonism, and a product term between selenium and mercury was in fact null (interaction P=0.90 for all-cause and P=0.22 for cardiovascular mortality; Table S11).',
 'R:sehq')

rep('In the smaller subset of never smokers with serum cotinine <10 ng/mL (n=7,821) the point estimates were similar but less precise (Table S3).',
 'In the smaller subset of never smokers with serum cotinine <10 ng/mL (n=7,821) the point estimates were similar but less precise (Table S3). The inverse mercury association was essentially unchanged when entry was delayed to 48 or 60 months (all-cause 0.854 and 0.854), which argues against reverse causation given mercury\u2019s short biological half-life (Table S10). Excluding observations at the assay floor left the all-cause association unchanged (0.86, 0.82\u20130.90) but moved the cardiovascular estimate to the null (0.92, 0.83\u20131.01; P=0.083) (Table S10) \u2014 that particular estimate should therefore be regarded as fragile. Substituting the blood-metal subsample weight for 2013\u20132016 and adding survey cycle as a covariate changed no estimate materially (Table S9), and the E-value for the inverse mercury association was 1.61 (1.49 using the confidence limit) for all-cause and 1.46 (1.16) for cardiovascular mortality. Estimates were consistent across the three periods (cadmium positive and mercury inverse in each period, with imprecise lead estimates in 2015\u20132018; Table S12).',
 'R:sens-extended')

# ---------- Discussion ----------
old_mer = s[s.index('The inverse association of blood mercury with mortality is consistent with some previous reports'):s.index('\n\nThe inverse and approximately monotonic association of selenium')]
new_mer = '''The inverse association of blood mercury with mortality is consistent with some previous reports [4\u20136,20,27]: serum mercury predicted a lower risk of death and myocardial infarction in Swedish women [5], meta-analyses of fish consumption suggest that the cardiovascular benefits of fish outweigh the risks of methylmercury at typical exposure levels [4,6], a large NHANES analysis found blood mercury inversely related to cardiovascular mortality [27], and serum mercury was inversely associated with all-cause and cardiovascular mortality in adults with cardiometabolic multimorbidity [20]. The pattern we observed, however, does not fit a simple cardioprotective, fish-related mechanism. First, and most importantly, the inverse association was confined to cancer and cerebrovascular deaths, whereas the estimate for heart disease \u2014 where a cardioprotective effect of fish would be expected to appear \u2014 was close to the null and not statistically significant (0.92, 0.84\u20131.02). A benefit mediated by fish intake should have shown up most clearly for ischaemic heart disease, so the cause-specific pattern argues either for residual confounding shared with other fish-associated behaviours or for a pathway unrelated to atherothrombosis. Second, the association was not explained by marine n-3 fatty acids or by selenium, although serum n-3 is only a partial proxy for fish intake and is measured with error, and a selenium\u2013mercury product term was null. Third, the association was robust to delaying entry into the risk set by 48\u201360 months, which is a fairly demanding test for a biomarker with a half-life of only about 50\u201370 days; this makes reverse causation less likely than it would be for a cumulative biomarker, although it does not exclude it. Fourth, the E-value of 1.61 indicates that an unmeasured confounder associated with both blood mercury and mortality by a risk ratio of about 1.6 \u2014 for example a healthier overall diet or a higher socioeconomic position among regular fish eaters \u2014 could account for the observed association. Blood mercury concentrations in the US general population are low and predominantly methylmercury from fish; in this range a causal protective effect of mercury itself remains biologically implausible, and we therefore do not interpret the association causally.'''
s = s.replace(old_mer, new_mer); print('ok D:mercury-para')

rep('are consistent with essential trace elements having optimal ranges, with both deficiency and excess being harmful [9,10]. Published evidence on manganese is not uniform: a cohort of dietary manganese intake reported lower cardiovascular mortality at higher intake [10], whereas we observed a U-shaped biomarker association \u2014 a discrepancy plausibly reflecting differences between dietary intake and blood biomarker, dose range and confounding structure.',
 'are consistent with essential trace elements having optimal ranges, with both deficiency and excess being harmful [9,23,26]. For manganese, our U-shaped pattern with a nadir near 8.3 \u00b5g/L is supported by an independent NHANES analysis reporting an inflection point near 7.1 \u00b5g/L [22]; by contrast, a cohort of dietary manganese intake found lower cardiovascular mortality at higher intake [10] \u2014 a discrepancy plausibly reflecting differences between dietary intake and blood biomarker, dose range and confounding structure. One caveat is that manganese enters quantile g-computation under a monotone, dose-weighted assumption that is inconsistent with a U-shaped relation, so its small positive weight in the five-metal mixture (0.05) should not be over-interpreted.',
 'D:mn-se')

rep('Third, we used MEC weights rather than blood-metal subsample weights because the latter were not available for all cycles; estimates are representative of the examined sample.',
 'Third, blood-metal subsample weights exist only for 2013\u20132016, so MEC weights were used for cross-cycle consistency; substituting the subsample weight for those two cycles did not change the estimates, and results are representative of the examined sample. Survey cycle was not included as a covariate because strata and primary sampling units are nested within cycle; adding it had no material effect (Table S9).',
 'D:lim-weights')

rep('In never smokers the number of cardiovascular deaths (581) limited precision, so those estimates are less stable than the all-cause ones.',
 'Ninth, a substantial fraction of cadmium measurements lay at the assay floor (19\u201350% depending on cycle; Table S8); excluding them did not weaken the cadmium association or the inverse mercury association for all-cause mortality, but it moved the mercury\u2013cardiovascular estimate to the null, and we therefore treat that estimate as fragile. Tenth, metals were measured in whole blood without haematocrit adjustment, although most of the metal mass in blood resides in erythrocytes. In never smokers the number of cardiovascular deaths (581) limited precision, so those estimates are less stable than the all-cause ones.',
 'D:lim-floor')

# ---------- Conclusion ----------
rep('In US adults, blood cadmium and lead were associated with higher all-cause and cardiovascular mortality and blood mercury with lower mortality. Because these associations act in opposite directions, the net mixture effect conceals both; separating the positive-direction and negative-direction partial effects is more informative than a single mixture estimate. The inverse mercury association spanned cancer and stroke and was not explained by marine n-3 fatty acids or selenium, and should not be interpreted causally. Blood selenium and manganese showed non-linear associations with mortality.',
 'In US adults, blood cadmium and lead were associated with higher all-cause and cardiovascular mortality and blood mercury with lower mortality. The inverse mercury association was confined to cancer and cerebrovascular deaths rather than heart disease \u2014 a pattern that does not fit a simply cardioprotective, fish-related mechanism \u2014 and was not explained by marine n-3 fatty acids or selenium. Because blood mercury largely reflects recent fish intake, reverse causation cannot be excluded, and the association should not be read as evidence that mercury is protective. Because cadmium and lead act in the opposite direction to mercury, the net mixture effect conceals both, and reporting positive- and negative-direction partial effects is more informative than a single mixture estimate. Blood selenium was inversely associated in an approximately monotonic manner within the observed range, whereas manganese showed a U-shaped association with a nadir near 8.3 \u00b5g/L.',
 'C:conclusion')

# ---------- References ----------
rep('20. Zhang A, Li Y, Wang X, et al. Associations of serum lead, cadmium, and mercury concentrations',
 '20. Zhang A, Wei P, Ding L, et al. Associations of serum lead, cadmium, and mercury concentrations', 'REF20')
rep('24. Wang Y, Zhang Y, Li S, et al. Low-grade systemic inflammation',
 '24. Wang Y, Wang Y, Li R, et al. Low-grade systemic inflammation', 'REF24')
rep('26. Chang CP, Lin YT, Chen YC, et al. Serum selenium and reduced mortality',
 '26. Chang CP, You CH. Serum selenium and reduced mortality', 'REF26')
rep('25. Reyes-Avila AD, Laws EA, Herrmann AD, DeLaune RD, Blanchard TP. Mercury and selenium levels, and Se:Hg molar ratios in freshwater fish from South Louisiana. *J Environ Sci Health A Tox Hazard Subst Environ Eng*. 2019;54(3):238\u2013245. doi:10.1080/10934529.2018.1546492. PMID 30601090.',
 '25. Ralston NVC, Raymond LJ. Dietary selenium\u2019s protective effects against methylmercury toxicity. *Toxicology*. 2010;278(1):112\u2013123. doi:10.1016/j.tox.2010.06.004. PMID 20561558.', 'REF25')
rep('26. Chang CP, You CH. Serum selenium and reduced mortality in middle-aged and older adults with prefrailty or frailty: the mediating role of inflammatory status. *Front Nutr*. 2025;12:1560167. doi:10.3389/fnut.2025.1560167. PMID 40791232.',
 '26. Chang CP, You CH. Serum selenium and reduced mortality in middle-aged and older adults with prefrailty or frailty: the mediating role of inflammatory status. *Front Nutr*. 2025;12:1560167. doi:10.3389/fnut.2025.1560167. PMID 40791232.\n27. Gonuguntla K, Alruwaili W, Thyagaturu H, et al. Impact of environmental heavy metal exposure on cardiovascular disease mortality. *Am J Prev Cardiol*. 2026;29:101488. doi:10.1016/j.ajpc.2026.101488. PMID 42403446.', 'REF27')

# ---------- Figure legends ----------
rep('**Figure S1.** Flow diagram of participant selection.',
 '**Figure S1.** Flow diagram of participant selection.\n\n**Figure S2.** Exploratory BKMR: overall joint effect of the three-metal mixture on all-cause mortality as a function of all metals being set at a given percentile, relative to the median (probit scale), with 95% credible bands (800 participants, 106 deaths; unweighted).\n\n**Figure S3.** Exploratory BKMR: univariate exposure\u2013response functions for blood lead, cadmium and mercury, with the other metals held at their median.', 'FIG-LEG')

# ---------- Supplementary table list ----------
rep('**Table S8.** Assay floor and proportion of the sample at the floor, by NHANES cycle and metal.',
 '''**Table S8.** Assay floor and proportion of the sample at the floor, by NHANES cycle and metal.

**Table S9.** Sensitivity analyses using the blood-metal subsample weight (WTSH2YR) for 2013\u20132016 and adding survey cycle as a covariate.

**Table S10.** Delayed entry into the risk set at 24, 48 and 60 months, and exclusion of observations at the assay floor.

**Table S11.** Selenium\u2013mercury interaction (product term and selenium tertiles) and mercury hazard ratios by selenium tertile.

**Table S12.** Period-specific estimates for 1999\u20132010, 2011\u20132014 and 2015\u20132018, and E-values for the inverse mercury association.''', 'SUPP-LIST')

# ---------- Data / code availability ----------
rep('**Data availability.** NHANES data are publicly available from the US National Center for Health Statistics (https://wwwn.cdc.gov/nchs/nhanes/). Analysis scripts are available from the corresponding author on reasonable request.',
 '**Data availability.** All NHANES data are publicly available from the US National Center for Health Statistics (https://wwwn.cdc.gov/nchs/nhanes/), and the public-use linked mortality files from https://www.cdc.gov/nchs/data-linkage/mortality-public.htm. **Code availability.** The complete analysis code (data download, cleaning, statistical analysis and figure generation) is openly available in a public repository [GitHub URL, to be inserted on submission] and will be archived on Zenodo with a persistent DOI upon acceptance.', 'DATA-AVAIL')

io.open(P, 'w', encoding='utf-8').write(s)
print('WROTE', P)
