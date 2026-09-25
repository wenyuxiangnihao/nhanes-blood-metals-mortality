#!/usr/bin/env python3
"""下载血重金属死亡分析的协变量文件（10 周期）+ 完整性校验."""
import os, urllib.request, time

RAW = "/sandbox/workspace/heavymetal/data_raw"
os.makedirs(RAW, exist_ok=True)
CYC = [(1999, ""), (2001, "_B"), (2003, "_C"), (2005, "_D"), (2007, "_E"),
       (2009, "_F"), (2011, "_G"), (2013, "_H"), (2015, "_I"), (2017, "_J")]
STEMS = ["DEMO", "BMX", "BPX", "SMQ", "DIQ", "BPQ", "MCQ", "ALQ", "PAQ", "COTNAL"]
# 生化(肌酐) / 血脂(总胆固醇、HDL) 周期命名不规则
BIOPRO = {1999: "LAB18", 2001: "L40_B", 2003: "L40_C", 2005: "BIOPRO_D", 2007: "BIOPRO_E",
          2009: "BIOPRO_F", 2011: "BIOPRO_G", 2013: "BIOPRO_H", 2015: "BIOPRO_I", 2017: "BIOPRO_J"}
CHOL = {1999: ["LAB13"], 2001: ["L13_B"], 2003: ["L13_C"],
        2005: ["TCHOL_D", "HDL_D"], 2007: ["TCHOL_E", "HDL_E"], 2009: ["TCHOL_F", "HDL_F"],
        2011: ["TCHOL_G", "HDL_G"], 2013: ["TCHOL_H", "HDL_H"], 2015: ["TCHOL_I", "HDL_I"],
        2017: ["TCHOL_J", "HDL_J"]}

def size_ok(p):
    return os.path.exists(p) and os.path.getsize(p) > 20000 and os.path.getsize(p) != 131072

def get(fn, yr):
    dest = f"{RAW}/{fn}.XPT"
    if size_ok(dest):
        return "skip"
    url = f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{yr}/DataFiles/{fn}.XPT"
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=120) as r:
                data = r.read()
            if len(data) < 20000 or len(data) == 131072:
                raise IOError(f"bad size {len(data)}")
            open(dest, "wb").write(data)
            return f"OK {len(data):,}"
        except Exception as e:
            last = str(e)
            time.sleep(2)
    return f"FAIL {last[:60]}"

ok = fail = 0
for yr, suf in CYC:
    for stem in STEMS:
        r = get(f"{stem}{suf}", yr)
        flag = r.split()[0]
        ok += flag in ("OK", "skip"); fail += flag == "FAIL"
        if flag == "FAIL":
            print(f"  !! {stem}{suf} ({yr}): {r}")
    r = get(BIOPRO[yr], yr)
    if r.split()[0] == "FAIL": print(f"  !! BIOPRO {BIOPRO[yr]}: {r}")
    for c in CHOL[yr]:
        r = get(c, yr)
        if r.split()[0] == "FAIL": print(f"  !! CHOL {c}: {r}")
print(f"\nfiles ok/skip={ok}  fail={fail}")
print("total files in data_raw:", len(os.listdir(RAW)))
