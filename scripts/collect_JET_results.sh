#!/usr/bin/env bash
set -euo pipefail

BASE="/home/junseokp/workspaces/projects/rtea/results"
DEST_ROOT="${1:-./collected_results}"

# File suffixes (prefix is <sample>)
FILES_SUFFIXES=(
  "_Fusions_chim2.junc2e7.size11.ids.txt"
  "_Fusions_chim2.junc2e7.size11.RData"
  "_Chimeric.out.annotatedJET.txt"
  "_Fusions_chim2.junc2e7.size11.genomic.txt"
  "_Log.progress.out"
  "_Log.out"
)

say() { printf "[%s] %s\n" "$(date '+%F %T')" "$*"; }

# Derive compact rel-path under /output/ keeping either:
#   nonReferenceTE/.../fq
#   referenceTE/.../fq
derive_rel_under_output() {
  local p="$1"
  [[ "$p" == *"/output/"* ]] || return 1
  local after_output="${p#*/output/}"

  local key=""
  if [[ "$after_output" == *"nonReferenceTE/"* ]]; then
    key="nonReferenceTE"
  elif [[ "$after_output" == *"referenceTE/"* ]]; then
    key="referenceTE"
  else
    return 1
  fi

  # start from "<key>/..."
  local rel="${key}/${after_output#*${key}/}"

  # cut after first "/fq"
  [[ "$rel" == *"/fq"* ]] || return 1
  rel="${rel%%/fq*}/fq"

  printf "%s" "$rel"
}

copy_dir_if_exists() {
  local src="$1"
  local dst="$2"
  if [[ -d "$src" ]]; then
    mkdir -p "$(dirname "$dst")"
    if command -v rsync >/dev/null 2>&1; then
      rsync -a --delete "$src"/ "$dst"/
    else
      rm -rf "$dst"
      cp -a "$src" "$dst"
    fi
  fi
}

say "Searching for fq directories under: $BASE"
mapfile -d '' FQ_DIRS < <(
  find "$BASE" -type d -name "fq" -path "*/output/*" -print0 2>/dev/null | sort -z
)

if [[ "${#FQ_DIRS[@]}" -eq 0 ]]; then
  say "ERROR: No fq directories found under $BASE"
  exit 1
fi

say "Found ${#FQ_DIRS[@]} fq directories. Collecting into: $DEST_ROOT"
mkdir -p "$DEST_ROOT"

for fq in "${FQ_DIRS[@]}"; do
  # sample = directory immediately before "output"
  sample="$(sed -n 's#^.*/\([^/]\+\)/output/.*#\1#p' <<<"$fq")"
  if [[ -z "${sample:-}" ]]; then
    say "WARN: Could not parse sample name from: $fq"
    continue
  fi

  rel="$(derive_rel_under_output "$fq" || true)"
  if [[ -z "${rel:-}" ]]; then
    say "WARN: Could not derive compact rel-path for fq dir: $fq"
    continue
  fi

  dest_dir="${DEST_ROOT}/${rel}/${sample}"
  mkdir -p "$dest_dir"

  say "----"
  say "Sample:  $sample"
  say "Run fq:   $fq"
  say "Dest dir: $dest_dir"

  # Copy requested files (search within fq up to depth 3)
  for suf in "${FILES_SUFFIXES[@]}"; do
    fname="${sample}${suf}"
    src_path="$(find "$fq" -maxdepth 3 -type f -name "$fname" -print -quit 2>/dev/null || true)"

    if [[ -n "$src_path" ]]; then
      cp -a "$src_path" "$dest_dir/"
      say "Copied: $fname"
    else
      if [[ "$suf" == "_Fusions_chim2.junc2e7.size11.ids.txt" ]]; then
        : > "${dest_dir}/${fname}"
        say "Missing ids.txt -> created empty: $fname"
      else
        say "WARN: Missing file (skipped): $fname"
      fi
    fi
  done

  # Copy err/log directories (try common nearby locations)
  fq_parent="$(dirname "$fq")"

  copy_dir_if_exists "$fq/log"        "$dest_dir/log"
  copy_dir_if_exists "$fq/err"        "$dest_dir/err"
  copy_dir_if_exists "$fq_parent/log" "$dest_dir/log"
  copy_dir_if_exists "$fq_parent/err" "$dest_dir/err"

  # bounded local search as fallback
  if [[ ! -d "$dest_dir/log" ]]; then
    near_log="$(find "$fq_parent" -maxdepth 3 -type d -name "log" -print -quit 2>/dev/null || true)"
    [[ -n "$near_log" ]] && copy_dir_if_exists "$near_log" "$dest_dir/log"
  fi
  if [[ ! -d "$dest_dir/err" ]]; then
    near_err="$(find "$fq_parent" -maxdepth 3 -type d -name "err" -print -quit 2>/dev/null || true)"
    [[ -n "$near_err" ]] && copy_dir_if_exists "$near_err" "$dest_dir/err"
  fi

  # Ensure log/err exist even if empty
  mkdir -p "$dest_dir/log" "$dest_dir/err"
done

say "Done."
say "Inspect with:"
say "  find \"$DEST_ROOT\" -maxdepth 6 -type f | sort"