# Analysis workflow represented by this repository

## Step 1: Guide-read processing and guide quantification

Inputs are guide-sequencing FASTQs, sample metadata, and the corresponding
library design. The output is a guide-by-sample count table used by MAGeCK and
the U-test workflow.

## Step 2: Genome-wide screen analysis

The genome-wide screen evaluates:

1. Day 3 CDH5-positive versus Day 3 CDH5-negative guide frequencies.
2. Day 3 total versus Day 0 guide frequencies as a survival-related comparison.

The documented methods are MAGeCK and a U-test analysis. Native result files
are kept outside Git and are listed in `docs/SOURCE_PATHS.tsv`.

## Step 3: Validation-screen analysis

The validation screen evaluates Day 3 CDH5-positive versus Day 3 CDH5-negative
fractions with MAGeCK and U-test workflows. Downstream code combines validation
and genome-wide effect estimates for summary plots and hit comparisons.

## Step 4: Reporting

R scripts generate volcano plots, correlation plots, rank plots, and additional
screen summaries. Plot outputs are not committed by default.

## Proposed IGVF workflow mapping

For IGVF metadata, the code can be represented by separate AnalysisStep records
for guide quantification, MAGeCK promoter-level enrichment, U-test
promoter-level enrichment, and optional reporting/cross-screen comparison.
The precise SoftwareVersion and AnalysisStepVersion records should be based on
versions and parameters confirmed from this repository snapshot and the
original run logs.
