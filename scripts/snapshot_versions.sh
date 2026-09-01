#!/usr/bin/env bash
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$ROOT/environment/software_versions_$(date +%Y%m%d_%H%M%S).txt"
{
  echo "Snapshot date: $(date -Is)"
  echo "Host: $(hostname)"
  echo
  for cmd in python python3 R Rscript mageck git; do
    echo "===== $cmd ====="
    if command -v "$cmd" >/dev/null 2>&1; then
      command -v "$cmd"
      case "$cmd" in
        R|Rscript) "$cmd" --version 2>&1 | head -n 4 ;;
        mageck) "$cmd" --version 2>&1 || "$cmd" -v 2>&1 || true ;;
        *) "$cmd" --version 2>&1 | head -n 4 ;;
      esac
    else
      echo "NOT FOUND"
    fi
    echo
  done
  if command -v conda >/dev/null 2>&1; then
    echo "===== conda list ====="
    conda list 2>&1 || true
  elif command -v mamba >/dev/null 2>&1; then
    echo "===== mamba list ====="
    mamba list 2>&1 || true
  fi
} > "$OUT"
echo "$OUT"
