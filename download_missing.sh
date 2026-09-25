#!/bin/bash
# 只下载缺失/可疑文件，提高并发
cd /sandbox/workspace/heavymetal/data_raw || exit 1
: > /tmp/missing.txt
while read -r u; do
  f=$(basename "$u")
  if [ -s "$f" ]; then
    s=$(stat -c%s "$f")
    if [ "$s" -gt 20000 ] && [ "$s" -ne 131072 ] && [ "$s" -ne 1323008 ]; then continue; fi
  fi
  echo "$u" >> /tmp/missing.txt
done < /tmp/urls.txt
echo "missing: $(wc -l < /tmp/missing.txt)"
xargs -a /tmp/missing.txt -P 10 -I {} curl -s --retry 8 --retry-all-errors --retry-delay 2 -m 400 -O \
  -w "%{http_code} %{size_download} %{filename_effective}\n" {} 2>&1
echo "MISSING_DONE"
