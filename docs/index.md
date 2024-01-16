# HarmoniSV

A toolkit to harmonize and filter structural variations across methods and samples.

*** Last updated: 2024-01-16 ***

## Features
- Harmonize SVs discovered by different SV calling methods
- Filter high-confidence SVs with a random forest classifier
- Fast VCF manipulation, annotation, and conversion

## Important
The documentation for below commands are ready now:

- [harmonize-header]
- [harmonize]
- [represent]
- [genotype]

For other commands, please type `harmonisv <command> -h` to get help. I am simplifying the input format and command line options for other commands. The documentation will be updated after that.

----
[harmonize]: VCF_manipulation/harmonize.md
[harmonize-header]: VCF_manipulation/harmonize_header.md
[sample2pop]: VCF_manipulation/sample2pop.md
[intersect]: VCF_manipulation/intersect.md
[represent]: SV_analysis/represent.md
[genotype]: SV_analysis/genotype.md
[filter]: SV_analysis/filter.md
[concordance]: SV_analysis/concordance.md