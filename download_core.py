#!/usr/bin/env python3
"""下载血重金属 LAB06/PBCD 全 10 周期 + 2019 公共使用死亡链接."""
import os, urllib.request, time

RAW = "/sandbox/workspace/heavymetal/data_raw"
os.makedirs(RAW, exist_ok=True)

PBCD = [(1999, "LAB06"), (2001, "L06B"), (2003, "L06C"),
        (2005, "PBCD_D"), (2007, "PBCD_E"), (2009, "PBCD_F"),
        (2011, "PBCD_G"), (2013, "PBCD_H"), (2015, "PBCD_I"), (2017, "PBCD_J")]

MORT = ["1999_2000", "2001_2002", "2003_2004", "2005_2006", "2007_2008",
        "2009_2010", "2011_2012", "2013_2014", "2015_2016", "2017_2018"]

def get(url, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        print(f"  skip {os.path.basename(dest)} ({os.path.getsize(dest)} B)")
        return True
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=120) as r, open(dest, "wb") as f:
            f.write(r.read())
        print(f"  OK   {os.path.basename(dest)}  {os.path.getsize(dest):,} B")
        return True
    except Exception as e:
        print(f"  FAIL {url} -> {e}")
        return False

print("== PBCD / blood metals ==")
for yr, fn in PBCD:
    url = f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{yr}/DataFiles/{fn}.XPT"
    get(url, f"{RAW}/{fn}.XPT")
    time.sleep(0.4)

print("== 2019 linked mortality (public use) ==")
for cyc in MORT:
    fn = f"NHANES_{cyc}_MORT_2019_PUBLIC"
    url = f"https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/{fn}.dat"
    get(url, f"{RAW}/{fn}.dat")
    time.sleep(0.4)

print("\n== inventory ==")
for f in sorted(os.listdir(RAW)):
    print(f"  {os.path.getsize(RAW+'/'+f):>12,}  {f}")
