# Quick start

*** Last updated: 2023-03-19 ***

## Installation
``` bash
git clone https://github.com/Han-Cao/HarmoniSV.git
cd HarmoniSV/src/harmoniSV
chmod +x harmonisv
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

## Usage
Go to the src/harmoniSV directory and run `harmonisv`:
``` bash
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
    2. HarmoniSV assume specific variant ID format to index SVs from different methods and samples, 
       please check the required ID format of different commands before you use
    3. The input/output VCF format (i.e., vcf, vcf.gz, bcf) will be automatically detected

For help on a specific command, run:
    harmonisv <command> -h

```