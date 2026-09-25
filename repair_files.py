#!/usr/bin/env python3
"""重下截断/可疑文件（先下到 .tmp，校验大小后再替换）."""
import os, urllib.request, time

RAW = "/sandbox/workspace/heavymetal/data_raw"
CDC = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public"
MORT = "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality"

TARGETS = [
    (f"{CDC}/1999/DataFiles/DEMO.XPT", "DEMO.XPT"),
    (f"{CDC}/1999/DataFiles/MCQ.XPT", "MCQ.XPT"),
    (f"{CDC}/2007/DataFiles/BPQ_E.XPT", "BPQ_E.XPT"),
    (f"{CDC}/2009/DataFiles/PBCD_F.XPT", "PBCD_F.XPT"),
    (f"{CDC}/2009/DataFiles/SMQ_F.XPT", "SMQ_F.XPT"),
    (f"{CDC}/2015/DataFiles/ALQ_I.XPT", "ALQ_I.XPT"),
] + [(f"{MORT}/NHANES_{c}_MORT_2019_PUBLIC.dat", f"NHANES_{c}_MORT_2019_PUBLIC.dat")
     for c in ["1999_2000","2001_2002","2005_2006","2007_2008","2009_2010",
               "2011_2012","2013_2014","2015_2016","2017_2018"]]

def rlen(u, tries=4):
    for _ in range(tries):
        try:
            req = urllib.request.Request(u, method="HEAD", headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=90) as r:
                return int(r.headers.get("Content-Length", -1))
        except Exception:
            time.sleep(3)
    return None

def fetch(u, dest, tries=6):
    for _ in range(tries):
        try:
            req = urllib.request.Request(u, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=600) as r, open(dest, "wb") as f:
                f.write(r.read())
            return True
        except Exception as e:
            time.sleep(3)
    return False

for u, fn in TARGETS:
    p = os.path.join(RAW, fn)
    rl = rlen(u)
    ls = os.path.getsize(p) if os.path.exists(p) else -1
    if rl is not None and rl == ls:
        print(f"  OK      {fn}  ({ls:,})"); continue
    print(f"  fixing  {fn}  remote={rl}  local={ls} ...", flush=True)
    tmp = p + ".tmp"
    if fetch(u, tmp):
        tl = os.path.getsize(tmp)
        if rl is None or tl == rl:
            os.replace(tmp, p); print(f"  FIXED   {fn}  -> {tl:,}")
        else:
            print(f"  BADSIZE {fn}  got {tl:,} expected {rl}")
            os.remove(tmp)
    else:
        print(f"  FAILED  {fn}")
print("\nrepair done")
