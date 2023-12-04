# HarmoniSV
A toolkit to harmonize and filter structural variations across methods and samples.

**Important**: The document is under development. We have tested HarmoniSV to perform population-scale SV calling using SV calls from `Sniffles2`, `cuteSV`, and `SVIM`. It should be able to work with any SV callers whose output follows VCF specification. Please open an issue for questions or bug reports.

## Features
- Harmonize SVs discovered by different SV calling methods
- Filter high-confidence SVs with a random forest classifier
- Fast VCF manipulation, annotation, and conversion

## Installation
``` bash
git clone https://github.com/Han-Cao/HarmoniSV.git
```
## Dependencies
HarmoniSV is written in python3.8. The following python modules are required:
```
pysam
pandas
numpy
matplotlib
scikit-learn
pyranges
```

## Quick start
``` bash
cd harmoniSV
./harmonisv

HarmoniSV: A toolkit to harmonize and filter structural variantions across methods and samples
Version: 0.1.0

Usage: harmonisv <command> [options]

Commands:

 -- VCF manipulation
    harmonize          Harmonize SV VCFs from across samples and SV calling methods
    harmonize-header   Harmonize VCF headers
    sample2pop         Convert single-sample VCF to multi-sample VCF
    intersect          Intersect SVs with genomic features

 -- Analysis on SV callset
    represent          Select the representative SV from merged SVs
    genotype           Genotype SVs across SV genotyping methods
    filter             Random forest filter for SVs
    concordance        Calculate genotype concordance between two VCFs


Note:
    1. All input VCFs MUST follow the VCF specification
    2. Some commands assume specific variant ID format to index SVs from different methods and samples, 
       please check the required ID format before you use
    3. The input/output VCF format (i.e., vcf, vcf.gz, bcf) will be automatically detected. However, a 
       temporary uncompressed VCF file will be generated if the output is vcf.gz or bcf.
     

For help on a specific command, run:
    harmonisv <command> -h

```
