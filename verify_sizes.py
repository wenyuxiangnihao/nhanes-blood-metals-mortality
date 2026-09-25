#!/usr/bin/env python3
"""对照服务器 Content-Length 校验本地文件完整性，输出需重下的清单."""
import os, re, urllib.request, json

RAW = "/sandbox/workspace/heavymetal/data_raw"
BASE = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public"
MORT = "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality"

CYC = [(1999, "", "LAB06"), (2001, "_B", "L06_B"), (2003, "_C", "L06C"),
       (2005, "_D", "PBCD_D"), (2007, "_E", "PBCD_E"), (2009, "_F", "PBCD_F"),
       (2011, "_G", "PBCD_G"), (2013, "_H", "PBCD_H"), (2015, "_I", "PBCD_I"),
       (2017, "_J", "PBCD_J")]
STEMS = ["DEMO", "BMX", "BPX", "SMQ", "DIQ", "BPQ", "MCQ", "ALQ", "PAQ"]
SPECIAL = {1999: ["LAB18", "LAB13"], 2001: ["L40_B", "L13_B"], 2003: ["L40_C", "L13_C"],
           2005: ["BIOPRO_D", "TCHOL_D", "HDL_D"], 2007: ["BIOPRO_E", "TCHOL_E", "HDL_E"],
           2009: ["BIOPRO_F", "TCHOL_F", "HDL_F"], 2011: ["BIOPRO_G", "TCHOL_G", "HDL_G"],
           2013: ["BIOPRO_H", "TCHOL_H", "HDL_H"], 2015: ["BIOPRO_I", "TCHOL_I", "HDL_I"],
           2017: ["BIOPRO_J", "TCHOL_J", "HDL_J"]}

urls = []
for yr, suf, _ in CYC:
    if suf == "_C":
        continue  # 2003-2004 has no blood metals; skip its covariates too
    for s in STEMS:
        urls.append(f"{BASE}/{yr}/DataFiles/{s}{suf}.XPT")
    for s in SPECIAL[yr]:
        urls.append(f"{BASE}/{yr}/DataFiles/{s}.XPT")
# blood metals
for yr, suf, pb in CYC:
    urls.append(f"{BASE}/{yr}/DataFiles/{pb}.XPT")
urls = sorted(set(urls))

def remote_len(u):
    try:
        req = urllib.request.Request(u, method="HEAD", headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=60) as r:
            return int(r.headers.get("Content-Length", -1))
    except Exception as e:
        return None

missing, mismatch, ok = [], [], 0
for u in urls:
    fn = u.split("/")[-1]
    p = os.path.join(RAW, fn)
    rl = remote_len(u)
    if rl is None or rl < 0:
        missing.append((u, fn, "HEAD_FAIL", os.path.getsize(p) if os.path.exists(p) else -1)); continue
    ls = os.path.getsize(p) if os.path.exists(p) else -1
    if ls != rl:
        mismatch.append((u, fn, rl, ls))
    else:
        ok += 1
print(f"checked {len(urls)} | OK {ok} | mismatch {len(mismatch)} | head-fail {len(missing)}")
print("\n--- MISMATCH (need re-download) ---")
for u, fn, rl, ls in mismatch:
    print(f"  {fn:<28} remote={rl:>10,}  local={ls:>10,}")
for u, fn, why, ls in missing:
    print(f"  {fn:<28} {why} local={ls}")
with open("/tmp/redownload.txt", "w") as f:
    for u, fn, rl, ls in mismatch:
        f.write(u + "\n")
print(f"\n{len(mismatch)} urls -> /tmp/redownload.txt")
