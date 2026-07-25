# *concordance*

Calculate genotype concordance between two VCFs

***Last updated: 2026-07-25***

## Overview

The `concordance` command compares genotype calls between two VCF files (`--invcf` and `--refvcf`) and computes per-variant concordance metrics. For each variant in the input VCF, a matching variant is located in the reference VCF by genomic position and either alleles or variant ID. This is useful to benchmark short-read SV genotyping against long-read SV calls using the same samples.

## Input requirements

- **Input VCF** (`--invcf`):
    - **Format**: bi-allelic, multi-sample VCF/BCF following the VCF specification.
    - **ID**: unique within the VCF. Used for matching if `--compare id` is specified.
    - **GT**: genotype field required for all samples to be compared.

- **Reference VCF** (`--refvcf`):
    - **Format**: bi-allelic, multi-sample VCF/BCF following the VCF specification. Must be sorted by position.
    - **GT**: genotype field required for all samples to be compared.

- **Sample map** (`--map`, optional):
    - **Format**: tab-separated file without header.
    - **Columns**: `input_sample`, `reference_sample`.
    - If not provided, samples with matching names in both VCFs are compared.

## Output

- **Output VCF**: `--outvcf`
    - A copy of the input VCF with three additional INFO fields:
        - `CON_GT`: proportion of compared samples with identical genotypes (0/0, 0/1, 1/1).
        - `CON_EXIST`: proportion of compared samples where the variant is present or absent in both VCFs.
        - `CON_WEIGHTED`: weighted genotype concordance — the average of concordance rates computed separately for each genotype class (homozygous reference, heterozygous, homozygous alternate). See [PanGenie](https://www.nature.com/articles/s41588-022-01043-w) paper for details on weighted concordance.

## Usage

``` bash
harmonisv concordance [options] -i <input_vcf> -r <reference_vcf> -o <output_vcf>
```

## Examples

##### 1. Basic genotype concordance between a test callset and a truth set

Samples with identical names in both VCFs are compared. Variants are matched by position and alleles. Only variants with `FILTER == PASS` in the reference VCF are considered.

``` bash
harmonisv concordance \
  -i test.vcf.gz \
  -r truth.vcf.gz \
  -o concordance.vcf.gz \
  --pass-only
```

##### 2. Concordance with sample name mapping and variant ID matching

When sample names differ between VCFs, use `--map` to specify the correspondence. Use `--compare id` if variant IDs are consistent but alleles may differ.

``` bash
harmonisv concordance \
  -i test.vcf.gz \
  -r truth.vcf.gz \
  -o concordance.vcf.gz \
  --map sample_map.txt \
  --compare id
```

```
HG002_test    HG002
HG003_test    HG003
HG004_test    HG004
```


## Arguments

#### Input/Output arguments:

-i, --invcf VCF
:   input VCF/BCF file (the callset to evaluate)

-r, --refvcf VCF
:   reference VCF/BCF file (the truth set to compare against)

-o, --outvcf VCF
:   output VCF/BCF file. A copy of `--invcf` with `CON_GT`, `CON_EXIST`, and `CON_WEIGHTED` added to INFO.

#### Comparison arguments:

--compare STR
:   method to match variants between the two VCFs. `alleles` (default): match by position and REF/ALT alleles. `id`: match by position and variant ID.

--map FILE
:   tab-separated file mapping input VCF sample names (column 1) to reference VCF sample names (column 2). If not provided, samples with identical names in both VCFs are compared.

--region STR
:   genomic region to compare, in bcftools format (e.g., `chr1` or `chr1:1000-2000`). Requires indexed VCFs.

--include-missing
:   include samples where the input VCF genotype is missing (`./.`) in the comparison. By default, only samples with non-missing genotypes in both VCFs are compared.

--pass-only
:   only consider variants with `FILTER == PASS` in the reference VCF. Variants at matching positions that do not pass filters are ignored.

#### Other arguments:

-h, --help
:   show help message and exit
