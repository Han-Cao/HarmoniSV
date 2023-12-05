# HarmoniSV

A toolkit to harmonize and filter structural variations across methods and samples.

*** Last updated: 2023-12-05 ***

## Features
- Harmonize SVs discovered by different SV calling methods
- Filter high-confidence SVs with a random forest classifier
- Fast VCF manipulation, annotation, and conversion

## Important
The documentation for below commands are ready now:

- [harmonize-header]
- [harmonize]
- [represent]

For other commands, please type `harmonisv <command> -h` to get help. I am simplifying the input format and command line options for other commands. The documentation will be updated after that.

----
[harmonize]: docs/VCF_manipulation/harmonize.md
[harmonize-header]: docs/VCF_manipulation/harmonize_header.md
[sample2pop]: docs/VCF_manipulation/sample2pop.md
[intersect]: docs/VCF_manipulation/intersect.md
[represent]: docs/SV_analysis/represent.md
[genotype]: docs/SV_analysis/genotype.md
[filter]: docs/SV_analysis/filter.md
[concordance]: docs/SV_analysis/concordance.md