#!/usr/bin/env bash
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

bad=0

# Reject large tracked files.
while IFS= read -r -d '' file; do
  size=$(stat -c %s "$file")
  if (( size > 20 * 1024 * 1024 )); then
    echo "LARGE TRACKED FILE: $file ($size bytes)" >&2
    bad=1
  fi
done < <(git ls-files -z)

# Reject common raw-data extensions if tracked.
if git ls-files | grep -Ei '\.(fastq|fq)(\.gz)?$|\.(bam|cram|sam)$|\.tar(\.gz)?$' >/dev/null; then
  echo "Raw sequencing/archive file appears to be tracked:" >&2
  git ls-files | grep -Ei '\.(fastq|fq)(\.gz)?$|\.(bam|cram|sam)$|\.tar(\.gz)?$' >&2
  bad=1
fi

# Look for common credential assignments. This is intentionally conservative.
if git grep -nEI '(IGVF_SECRET_KEY|ANTHROPIC_API_KEY|GITHUB_TOKEN|ghp_[A-Za-z0-9]{20,}|github_pat_)' -- . \
   ':(exclude)scripts/prepush_safety_check.sh' >/tmp/igvf_repo_secret_scan.$$ 2>/dev/null; then
  echo "Possible credential material detected:" >&2
  cat /tmp/igvf_repo_secret_scan.$$ >&2
  bad=1
fi
rm -f /tmp/igvf_repo_secret_scan.$$

if (( bad )); then
  echo "Safety check FAILED." >&2
  exit 1
fi

echo "Safety check passed."
