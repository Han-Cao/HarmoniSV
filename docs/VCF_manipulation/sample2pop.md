# *sample2pop*

Convert per-sample VCFs to a multi-sample population VCF

***Last updated: 2026-07-25***

## Overview

The `sample2pop` command merges multiple per-sample VCFs into a single multi-sample (population) VCF. Variants are merged across samples by their ID, and INFO and FORMAT fields are merged according to user-defined rules.

## Input requirements

- **Manifest**:
    - **Format**: tab-separated file with header.
    - **Columns**:
        - `file`: path to per-sample VCF/BCF file
        - `sample`: sample name to be used in the output VCF. Sample names must be unique.

- **Per-sample VCF**:
    - **Format**: bi-allelic VCF/BCF files following the VCF specification.
    - **ID**: unique within each per-sample VCF. Variants are merged across samples by ID.

## Output

- **Multi-sample VCF**: a single VCF/BCF file containing all samples.
    - **ID**: same as the per-sample VCF ID.
    - **INFO**: merged according to user-specified rules (see below). `AC` and `AN` are automatically recalculated as the sum across samples.
    - **FORMAT**: per-sample fields. `GT` is always included; `FT_SAMPLE` records the per-sample FILTER status. Additional FORMAT tags are preserved via `--keep-format`, and INFO tags can be moved to FORMAT via `--info-to-format`.

## Merge rules

For each INFO field and variant-level attribute, the user specifies a merge rule that determines how values from multiple samples are combined. Fields not mentioned in any rule are dropped. Rules accept wildcards (`*`) to match groups of INFO tags.

| Rule | Option | Description |
|------|--------|-------------|
| First | `--info-first` | Use the first non-missing value across samples |
| Sum | `--info-sum` | Sum of values across samples (missing treated as 0) |
| Average | `--info-avg` | Mean of values across samples (missing excluded) |
| Minimum | `--info-min` | Minimum value across samples (missing excluded) |
| Maximum | `--info-max` | Maximum value across samples (missing excluded) |
| Concatenate | `--info-concat` | Comma-separated string of all non-missing values |

## INFO to FORMAT

The `--info-to-format` option moves specified INFO tags into the per-sample FORMAT field. This is useful for retaining useful information of individual SV calls (e.g., `RF_SCORE`).

The per-sample FILTER status is always recorded in `FORMAT/FT_SAMPLE` as a string. If `--filter-GT` is used, genotypes of variants that did not pass filters are set to `./.`

## Usage

``` bash
harmonisv sample2pop [options] --manifest <manifest.txt> -o <output.vcf.gz>
```

## Examples

##### 1. Basic merge of two per-sample filtered VCFs

The following command merges two per-sample VCFs into a multi-sample VCF, with the `GT` FORMAT field and the `SVTYPE`, `SVLEN`, and `END` INFO fields unchanged. All depth-related INFO fields (`DP*` and `RE*`) are summed across samples, and all support-related INFO fields (`SUPP_*`) are averaged across samples. The maximum `RF_SCORE` value per variant is used as site-level `RF_SCORE` in `INFO`, and `RF_SCORE` of each SV call are moved to the per-sample FORMAT field.

``` bash
harmonisv sample2pop \
--manifest manifest_sample2pop.txt \
--outvcf population.vcf.gz \
--filter-GT \
--info-first SVTYPE,SVLEN,END \
--info-sum "DP*,RE*" \
--info-avg "SUPP_*" \
--info-max RF_SCORE \
--info-to-format RF_SCORE \
--keep-format GT
```

## Arguments

#### Input/Output arguments:

--manifest FILE
:   tab-delimited manifest file with columns: `file` (path to per-sample VCF), `sample` (sample name). Sample names must be unique.

-o, --outvcf FILE
:   output VCF/BCF file. Supports `.vcf`, `.vcf.gz` (bgzip-compressed), and `.bcf` extensions.

#### Merge rule arguments (all accept wildcards `*`):

--info-first TAGS
:   comma-separated INFO fields merged by taking the first non-missing value

--info-sum TAGS
:   comma-separated INFO fields merged by summing across samples

--info-avg TAGS
:   comma-separated INFO fields merged by averaging across samples

--info-min TAGS
:   comma-separated INFO fields merged by taking the minimum across samples

--info-max TAGS
:   comma-separated INFO fields merged by taking the maximum across samples

--info-concat TAGS
:   comma-separated INFO fields merged by concatenating non-missing values (comma-separated)

#### FORMAT arguments:

--keep-format TAGS
:   comma-separated FORMAT fields to preserve in the output (default: `GT`). Accepts wildcards.

--info-to-format TAGS
:   comma-separated INFO fields to move to per-sample FORMAT. Each tag is removed from merged INFO and written per sample.

--filter-GT
:   set genotype to `./.` if the per-sample FILTER is not missing and not `PASS`. The original FILTER value is recorded in `FORMAT/FT_SAMPLE`.

#### Other arguments:

-r, --region STR
:   genomic region to process, in bcftools format (e.g., `chr1` or `chr1:1000-2000`). Requires indexed VCFs.
