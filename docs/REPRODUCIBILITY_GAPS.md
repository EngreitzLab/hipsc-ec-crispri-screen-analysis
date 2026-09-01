# Reproducibility gaps to resolve before final IGVF workflow submission

The offboarding document identifies code and result locations but does not fully
specify every executable command. Confirm the following from copied notebooks,
shell histories, run logs, or the original analysts:

1. Exact MAGeCK version.
2. Exact `mageck count` and `mageck test` commands.
3. Treatment/control sample names and replicate grouping for each comparison.
4. Normalization method and any control-guide list.
5. Filtering thresholds for guides or promoters.
6. Guide-to-promoter aggregation rules, including P1/P2/P1P2 handling.
7. Exact U-test implementation, grouping rule, sidedness, and multiple-testing
   correction.
8. The normalization factor used in the validation summary and how it was
   derived.
9. R and Python package versions used for plots and table assembly.
10. Whether the additional sequencing directories were merged or one was a
    redundant rerun.
11. Resolution of the documented `NEDD8` duplicate-name issue.
12. A definitive primary output for the validation-screen per-promoter MAGeCK
    analysis rather than only a gene-level combined summary.

Do not fill these gaps by guessing. Record confirmed values in a release or an
IGVF AnalysisStepVersion description.
