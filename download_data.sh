#!/usr/bin/env bash
# =============================================================================
# download_data.sh
#   Fetch every public NHANES file used by this study plus the NCHS 2019
#   Linked Mortality File (public-use, fixed-width).
#
# Usage:
#   ./download_data.sh                 # writes into ./data_raw
#   RAW_DIR=/mnt/nhanes ./download_data.sh
#
# NHANES file URL pattern:
#   https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{startyear}/DataFiles/{FILE}.XPT
# Linked Mortality files:
#   https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW="${RAW_DIR:-$HERE/data_raw}"
BASE="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public"
MORTBASE="https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality"
mkdir -p "$RAW"

file_size() { stat -c%s "$1" 2>/dev/null || stat -f%z "$1" 2>/dev/null || echo 0; }

dl() {  # $1 = FILE stem, $2 = start year
  local fn="$1" yr="$2" out="$RAW/$1.XPT"
  if [ -s "$out" ] && [ "$(file_size "$out")" -gt 20000 ]; then
    echo "  skip $1.XPT"; return 0
  fi
  if curl -fsSL --retry 4 --retry-delay 2 -o "$out" "$BASE/$yr/DataFiles/$1.XPT"; then
    echo "  ok   $1.XPT"
  else
    echo "  FAIL $1.XPT"; rm -f "$out"; return 1
  fi
}

dl_opt() { dl "$1" "$2" || true; }   # optional file, never abort

dlm() { # $1 = cycle label (e.g. 1999_2000)
  local fn="NHANES_${1}_MORT_2019_PUBLIC"
  local out="$RAW/${fn}.dat"
  if [ -s "$out" ]; then echo "  skip $fn.dat"; return 0; fi
  if curl -fsSL --retry 4 --retry-delay 2 -o "$out" "$MORTBASE/${fn}.dat"; then
    echo "  ok   $fn.dat"
  else
    echo "  FAIL $fn.dat"; rm -f "$out"; return 1
  fi
}

# start-year:suffix  ("" for 1999-2000, "_B" for 2001-2002, ... "_J" for 2017-2018)
CYCLES=("1999:" "2001:_B" "2003:_C" "2005:_D" "2007:_E" \
        "2009:_F" "2011:_G" "2013:_H" "2015:_I" "2017:_J")

for c in "${CYCLES[@]}"; do
  yr="${c%%:*}"; suf="${c##*:}"
  echo "== cycle $yr (suffix '${suf:-none}') =="

  # --- blood metals: Pb (LBXBPB), Cd (LBXBCD), Hg (LBXTHG), Se (LBXBSE), Mn (LBXBMN)
  case "$yr" in
    1999) dl LAB06    "$yr" ;;
    2001) dl L06_B    "$yr" ;;
    2003) dl L06BMT_C "$yr" ;;
    *)    dl "PBCD${suf}" "$yr" ;;
  esac

  # --- serum cotinine: 1999/2001 ship inside LAB06 / L06_B; 2003 = L06COT_C; 2005+ = COT_x
  #     2007-2008 does not publish COT_E; the cycle ships COTNAL_E instead (fetched
  #     as an optional file below). These are tolerated as optional because cotinine
  #     is not part of the complete-case covariate set of the primary analysis.
  case "$yr" in
    1999|2001) : ;;
    2003)      dl L06COT_C "$yr" ;;
    *)         dl_opt "COT${suf}" "$yr" ;;
  esac

  # --- demographics / examination / questionnaire
  dl "DEMO${suf}" "$yr"     # age, sex, race/ethnicity, education, PIR, weights, PSU, strata, pregnancy
  dl "BMX${suf}"  "$yr"     # BMI, waist
  dl "BPX${suf}"  "$yr"     # blood pressure (systolic / diastolic)
  dl "SMQ${suf}"  "$yr"     # smoking
  dl "DIQ${suf}"  "$yr"     # diabetes
  dl "BPQ${suf}"  "$yr"     # hypertension
  dl "MCQ${suf}"  "$yr"     # self-reported CVD / cancer
  dl "ALQ${suf}"  "$yr"     # alcohol
  dl "PAQ${suf}"  "$yr"     # physical activity
  dl_opt "COTNAL${suf}" "$yr"   # (cotinine, questionnaire variant; not used)

  # --- serum chemistry: creatinine (for eGFR)
  case "$yr" in
    1999) dl LAB18  "$yr" ;;
    2001) dl L40_B  "$yr" ;;
    2003) dl L40_C  "$yr" ;;
    *)    dl "BIOPRO${suf}" "$yr" ;;
  esac

  # --- lipids: total cholesterol + HDL
  case "$yr" in
    1999) dl LAB13 "$yr" ;;
    2001) dl L13_B "$yr" ;;
    2003) dl L13_C "$yr" ;;
    *)    dl "TCHOL${suf}" "$yr"; dl "HDL${suf}" "$yr" ;;
  esac
done

echo "== optional: alternative heavy-metal files (2001/2003) & fatty acids (marine n-3) =="
dl_opt L06HM_B 2001     # 2001-2002 heavy metals (alternate release)
dl_opt L06HM_C 2003     # 2003-2004 heavy metals (alternate release)
dl_opt FAS_G   2011     # 2011-2012 fatty acids (EPA/DHA)  -> marine n-3 sensitivity
dl_opt FAS_H   2013     # 2013-2014 fatty acids (EPA/DHA)

echo "== NCHS 2019 Linked Mortality File (public-use, fixed-width) =="
for cyc in 1999_2000 2001_2002 2003_2004 2005_2006 2007_2008 \
           2009_2010 2011_2012 2013_2014 2015_2016 2017_2018; do
  dlm "$cyc"
done

echo
echo "All files downloaded into: $RAW"
echo "For the Linked Mortality File documentation and the official SAS/R read-in"
echo "program, see: https://www.cdc.gov/nchs/data-linkage/mortality-public.htm"
