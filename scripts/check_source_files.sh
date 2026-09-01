#!/usr/bin/env bash
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TABLE="$ROOT/docs/SOURCE_PATHS.tsv"

printf 'status\tscreen\tcategory\trole\tpath\n'
tail -n +2 "$TABLE" | while IFS=$'\t' read -r screen category role path; do
  if [[ -e "$path" ]]; then
    status=FOUND
  else
    status=MISSING
  fi
  printf '%s\t%s\t%s\t%s\t%s\n' "$status" "$screen" "$category" "$role" "$path"
done
