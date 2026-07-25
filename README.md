# HarmoniSV
A toolkit to harmonize and filter structural variations across methods and samples.

- [Features](#features)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Documentation](#documentation)
- [License](#license)

## Features
- Harmonize SVs discovered by different SV calling methods
- Filter high-confidence SVs with a random forest classifier
- Fast VCF manipulation, annotation, and conversion for SVs

## Installation

Clone the repository:
``` bash
git clone https://github.com/Han-Cao/HarmoniSV.git
```

Make sure you have Python 3.8+ and pip installed. Then, install the tool and its dependencies:
``` bash
pip install .
```

## Quick start
``` bash
harmonisv

HarmoniSV: A toolkit to harmonize and filter structural variantions across methods and samples
Version: 0.1.0

Usage: harmonisv <command> [options]

Commands:

 -- VCF manipulation
    harmonize          Harmonize SV VCFs across samples and SV calling methods
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
    2. Some commands assume unique SV IDs to index SVs
    3. The input/output VCF format (i.e., vcf, vcf.gz, bcf) will be automatically detected. However, a 
       temporary uncompressed VCF file will be generated if the output is vcf.gz or bcf
              

For help on a specific command, run:
    harmonisv <command> -h

```

## Documentation

For more details and tutorials, please refer to the [documentation].

## License

[MIT License](LICENSE)

[documentation]: https://han-cao.github.io/HarmoniSV/