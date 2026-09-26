#!/bin/bash
# 协变量批量下载（curl 版，带完整性校验）
cd ${NHANES_ROOT:-$(pwd)}/data_raw || exit 1
BASE="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public"
OK=0; FAIL=0; FAILED=""

dl() {
  local fn="$1" yr="$2"
  local f="${fn}.XPT"
  if [ -s "$f" ] && [ "$(stat -c%s "$f")" -gt 20000 ] && [ "$(stat -c%s "$f")" -ne 131072 ]; then
    OK=$((OK+1)); return
  fi
  for i in 1 2 3 4; do
    curl -s --retry 4 --retry-all-errors --retry-delay 2 -m 100 -o "$f" "$BASE/$yr/DataFiles/$f"
    local s=$(stat -c%s "$f" 2>/dev/null || echo 0)
    if [ "$s" -gt 20000 ] && [ "$s" -ne 131072 ]; then
      echo "  OK   $f ($s)"; OK=$((OK+1)); return
    fi
    sleep 1
  done
  echo "  FAIL $f"; FAIL=$((FAIL+1)); FAILED="$FAILED $f"
}

# 周期: "起始年 后缀"
CYCLES=("1999 " "2001 _B" "2003 _C" "2005 _D" "2007 _E" "2009 _F" "2011 _G" "2013 _H" "2015 _I" "2017 _J")
STEMS=(DEMO BMX BPX SMQ DIQ BPQ MCQ ALQ PAQ COTNAL)

for c in "${CYCLES[@]}"; do
  set -- $c; yr=$1; suf=$2
  for s in "${STEMS[@]}"; do dl "${s}${suf}" "$yr"; done
  case $yr in
    1999) dl LAB18 $yr; dl LAB13 $yr;;
    2001) dl L40_B $yr; dl L13_B $yr;;
    2003) dl L40_C $yr; dl L13_C $yr;;
    2005) dl BIOPRO_D $yr; dl TCHOL_D $yr; dl HDL_D $yr;;
    2007) dl BIOPRO_E $yr; dl TCHOL_E $yr; dl HDL_E $yr;;
    2009) dl BIOPRO_F $yr; dl TCHOL_F $yr; dl HDL_F $yr;;
    2011) dl BIOPRO_G $yr; dl TCHOL_G $yr; dl HDL_G $yr;;
    2013) dl BIOPRO_H $yr; dl TCHOL_H $yr; dl HDL_H $yr;;
    2015) dl BIOPRO_I $yr; dl TCHOL_I $yr; dl HDL_I $yr;;
    2017) dl BIOPRO_J $yr; dl TCHOL_J $yr; dl HDL_J $yr;;
  esac
done
echo "=== DONE ok=$OK fail=$FAIL ==="
[ -n "$FAILED" ] && echo "FAILED:$FAILED"
echo "total XPT: $(ls *.XPT | wc -l)"
