# hiPSC-EC genome-wide and validation CRISPRi screen analysis

This repository is a provenance snapshot of the analysis code used for the
hiPSC-to-endothelial-cell genome-wide CRISPRi differentiation screen and its
focused validation screen.

## Biological comparisons

### Genome-wide screen

- Primary differentiation comparison: Day 3 CDH5-positive versus Day 3
  CDH5-negative fractions.
- Survival comparison: Day 3 total (combined CDH5-positive and CDH5-negative)
  versus Day 0.
- Analysis methods documented in the source handoff: MAGeCK and a U-test
  workflow.

### Validation screen

- Primary comparison: Day 3 CDH5-positive versus Day 3 CDH5-negative fractions.
- Analysis methods documented in the source handoff: MAGeCK and a U-test
  workflow, followed by cross-screen summaries and volcano/correlation plots.

## Repository organization

- `genome_wide/`: copied genome-wide analysis notebooks and plotting scripts.
- `validation/`: copied validation-screen guide-counting, MAGeCK, U-test, and
  plotting code.
- `docs/SOURCE_PATHS.tsv`: source inputs, outputs, and code locations recorded
  in the offboarding document.
- `docs/CODE_SNAPSHOT_MANIFEST.tsv`: exact source-to-repository copy manifest,
  file sizes, and SHA-256 checksums.
- `docs/REPRODUCIBILITY_GAPS.md`: items that still need confirmation before the
  repository can be described as a fully reproducible rerun.
- `config/paths.example.env`: configurable Sherlock/OAK paths.
- `scripts/check_source_files.sh`: verifies documented Sherlock paths.
- `scripts/snapshot_versions.sh`: records currently available software versions.
- `scripts/prepush_safety_check.sh`: rejects obvious raw-data or credential
  files before a push.

## Data policy

Raw FASTQs, count matrices, MAGeCK result tables, and other large or controlled
files are not committed. The analysis code expects those files to remain on
Sherlock/OAK or to be supplied through paths in `config/paths.example.env`.
