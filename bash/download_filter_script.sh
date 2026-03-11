#!/usr/bin/env bash
# Download & subset REMSS CCMP .nc files for a year/month range, preserving source folder structure.
# After successful local subsetting the original raw file is removed (to save space).
set -euo pipefail
shopt -s globstar nullglob

BASE_URL="https://data.remss.com/ccmp/v03.1"   # root on the host
START_YEAR=2018
END_YEAR=2019

# months: use quoted elements
MONTHS=( "01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" )

RAW_ROOT="raw"
SUBSET_ROOT="subset"

# bounding box (edit to your desired lat/lon)
LONMIN=-40
LONMAX=8
LATMIN=8
LATMAX=40

# If CLEAN_RAW=1 the raw downloaded file will be deleted after a successful subset to save disk space.
CLEAN_RAW=1

WAIT=1
RANDOM_WAIT="--random-wait"
LIMIT_RATE=""

VERBOSE=1

mkdir -p "$RAW_ROOT" "$SUBSET_ROOT"

log() { [ "$VERBOSE" -eq 1 ] && printf '%s\n' "$*"; }

# Helper: list .nc URLs for a single remote directory (year/month dir)
list_nc_urls_from_dir() {
  local dirurl="$1"
  wget --spider -r -l1 -np -nv "$dirurl" 2>&1 \
    | grep -Eo 'https?://[^[:space:]]+\.nc' \
    | sort -u
}

for year in $(seq "$START_YEAR" "$END_YEAR"); do
  for month in "${MONTHS[@]}"; do
    dir_url="${BASE_URL}/Y${year}/M${month}/"
    log "Listing: $dir_url"
    urls=$(list_nc_urls_from_dir "$dir_url") || urls=""
    if [ -z "$urls" ]; then
      log "  (no .nc URLs found at $dir_url) — skipping"
      continue
    fi

    while IFS= read -r url; do
      [ -z "$url" ] && continue
      path="$(echo "$url" | sed -E 's#https?://[^/]+##')"
      path="${path#/}"

      out_path="$SUBSET_ROOT/$path"
      raw_path="$RAW_ROOT/$path"
      mkdir -p "$(dirname "$out_path")" "$(dirname "$raw_path")"

      if [ -f "$out_path" ]; then
        log "skip (already exists): $out_path"
        continue
      fi

      log "Processing: $url"
      log "  target subset: $out_path"

      # Try remote cdo first (no raw file created)
      if cdo sellonlatbox,${LONMIN},${LONMAX},${LATMIN},${LATMAX} "$url" "$out_path" >/dev/null 2>&1; then
        log "  remote cdo succeeded -> $out_path"
        continue
      fi

      log "  remote cdo failed; downloading and running local cdo"

      tmpdl="${raw_path}.part"
      mkdir -p "$(dirname "$tmpdl")"
      wget -c -nv ${LIMIT_RATE:+--limit-rate="$LIMIT_RATE"} --wait="$WAIT" ${RANDOM_WAIT:+$RANDOM_WAIT} -O "$tmpdl" "$url"
      mv -f "$tmpdl" "$raw_path"

      if [ ! -f "$raw_path" ]; then
        log "  download failed for $url, skipping"
        continue
      fi

      if cdo sellonlatbox,${LONMIN},${LONMAX},${LATMIN},${LATMAX} "$raw_path" "$out_path" >/dev/null 2>&1; then
        log "  local cdo succeeded -> $out_path"

        # CLEAN-UP: remove the raw downloaded file to save disk space if requested
        if [ "${CLEAN_RAW:-0}" -eq 1 ]; then
          rm -f "$raw_path" || true
          # attempt to remove parent dir if empty (ignore errors)
          rmdir --ignore-fail-on-non-empty "$(dirname "$raw_path")" 2>/dev/null || true
          # also try to remove empty parent trees up to RAW_ROOT (optional)
          parent="$(dirname "$raw_path")"
          while [ "$parent" != "." ] && [ "$parent" != "$RAW_ROOT" ] && [ "$parent" != "/" ]; do
            rmdir --ignore-fail-on-non-empty "$parent" 2>/dev/null || break
            parent="$(dirname "$parent")"
          done
        fi

      else
        log "  local cdo FAILED on $raw_path"
      fi

    done <<< "$urls"
  done
done

log "All done. Subsets are under: $SUBSET_ROOT (raw originals under: $RAW_ROOT)"
