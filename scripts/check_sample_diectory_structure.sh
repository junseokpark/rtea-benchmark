#!/bin/bash
set -euo pipefail

DATA_HOME="/home/junseokp/workspaces/data/rTea-simul/sims"

cd "$DATA_HOME" || { echo "DATA_HOME directory not found: $DATA_HOME"; exit 1; }

# Get directories containing *.fq.gz (GNU find: -printf; BSD/macOS fallback: dirname)
if find . -type f -name "*.fq.gz" -printf "%h\n" >/dev/null 2>&1; then
  list_cmd='find . -type f -name "*.fq.gz" -printf "%h\n"'
else
  list_cmd='find . -type f -name "*.fq.gz" -exec dirname {} \;'
fi

eval "$list_cmd" \
| awk '
  {
    n=split($0,a,"/");
    p=a[1];
    seen["."]=1;
    for(i=2;i<=n;i++){ p=p"/"a[i]; seen[p]=1 }
  }
  END{ for(p in seen) print p }
' \
| LC_ALL=C sort -u \
| tree --fromfile -d