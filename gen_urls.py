#!/usr/bin/env python3
import os
BASE = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public"
CYC = [(1999, ""), (2001, "_B"), (2003, "_C"), (2005, "_D"), (2007, "_E"),
       (2009, "_F"), (2011, "_G"), (2013, "_H"), (2015, "_I"), (2017, "_J")]
STEMS = ["DEMO", "BMX", "BPX", "SMQ", "DIQ", "BPQ", "MCQ", "ALQ", "PAQ"]
# 血清可替宁：1999/2001 在 LAB06/L06_B 内(LBXCOT)；2003=L06COT_C；2005+=COT_x
COT = {1999: [], 2001: [], 2003: ["L06COT_C"], 2005: ["COT_D"], 2007: ["COT_E"],
       2009: ["COT_F"], 2011: ["COT_G"], 2013: ["COT_H"], 2015: ["COT_I"], 2017: ["COT_J"]}
SPECIAL = {
    1999: ["LAB18", "LAB13"],
    2001: ["L40_B", "L13_B"],
    2003: ["L40_C", "L13_C"],
    2005: ["BIOPRO_D", "TCHOL_D", "HDL_D"],
    2007: ["BIOPRO_E", "TCHOL_E", "HDL_E"],
    2009: ["BIOPRO_F", "TCHOL_F", "HDL_F"],
    2011: ["BIOPRO_G", "TCHOL_G", "HDL_G"],
    2013: ["BIOPRO_H", "TCHOL_H", "HDL_H"],
    2015: ["BIOPRO_I", "TCHOL_I", "HDL_I"],
    2017: ["BIOPRO_J", "TCHOL_J", "HDL_J"],
}
lines = []
for yr, suf in CYC:
    for s in STEMS:
        lines.append(f"{BASE}/{yr}/DataFiles/{s}{suf}.XPT")
    for s in SPECIAL[yr]:
        lines.append(f"{BASE}/{yr}/DataFiles/{s}.XPT")
    for s in COT[yr]:
        lines.append(f"{BASE}/{yr}/DataFiles/{s}.XPT")
open("/tmp/urls.txt", "w").write("\n".join(sorted(set(lines))) + "\n")
print(f"{len(lines)} urls written")
for u in lines:
    print(" ", u.split('/DataFiles/')[1])
