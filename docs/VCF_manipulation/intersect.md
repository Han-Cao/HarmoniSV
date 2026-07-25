# *intersect*

Intersect structural variants with genomic features

***Last updated: 2026-07-25***

## Overview

The `intersect` command annotates structural variants with overlap information against user-provided genomic feature BED files. For each SV, the tool computes either a binary flag or the number of overlapping base pairs. Results can be written as annotated VCF or as a tab-separated TSV file.

## Input requirements

- **VCF**:
    - **Format**: bi-allelic VCF/BCF following the VCF specification.
    - **INFO**: `SVTYPE` and `SVLEN` are required. For insertions (`SVTYPE=INS`), the SV is treated as a single base-pair interval at the start position.
    - **ID**: unique within the VCF.

- **Feature BED files**:
    - **Format**: standard 3-column BED (CHR, START, END). 0-based half-open coordinates.
    - Multiple features can be specified as a comma-separated list to `--bed`.

## Output

- **VCF output** (`--out-type vcf`, default): a copy of the input VCF with one additional INFO field per feature. For binary mode, the value is 0 (no overlap) or 1 (overlap). For amount mode, the value is the number of overlapping base pairs.

- **TSV output** (`--out-type tsv`): a tab-separated file with columns `ID`, `CHR`, `POS`, `SVTYPE`, `SVLEN`, and feature overlap values.

## Overlap modes

Each feature is configured independently via the `--mode` and `--overlap` options:

- **Binary mode** (`b`, default): reports 1 if the SV overlaps the feature by at least `--overlap` base pairs, 0 otherwise.

- **Amount mode** (`a`): reports the number of overlapping base pairs. Only overlaps meeting the `--overlap` threshold are reported; others are set to 0.

- **Overlap threshold**: the `--overlap` value specifies the minimum required overlap in base pairs. A value of `-1` requires the SV to be fully contained within the feature (i.e., overlap equals SV length).

## Usage

``` bash
harmonisv intersect [options] -i <input_vcf> --bed <feature.bed> --name <feature_name> -o <output>
```

## Examples

##### 1. Binary annotation of a single feature

Flag each SV that overlaps a segmental duplication with at least 1 bp:

``` bash
harmonisv intersect \
  -i sv_calls.vcf.gz \
  --bed segdup.bed \
  --name SEGDUP \
  -o sv_calls.segdup.vcf.gz
```

The output VCF will contain `INFO/SEGDUP=1` for SVs overlapping a segmental duplication, and `0` otherwise.

##### 2. Amount annotation with a minimum overlap threshold

Compute the number of base pairs each SV overlaps with a repetitive element, requiring at least 50 bp of overlap:

``` bash
harmonisv intersect \
  -i sv_calls.vcf.gz \
  --bed repeats.bed \
  --name REPEAT_OVERLAP \
  --mode a \
  --overlap 50 \
  -o sv_calls.repeats.vcf.gz
```

The output VCF will contain `INFO/REPEAT_OVERLAP` with the overlap amount in base pairs.

##### 3. Multiple features with different modes

Annotate SVs with three features simultaneously — two binary and one amount-based:

``` bash
harmonisv intersect \
  -i sv_calls.vcf.gz \
  --bed exons.bed,repeats.bed,centromeres.bed \
  --name EXON,REPEAT,CENTROMERE \
  --mode b,a,b \
  --overlap 1,50,-1 \
  -o sv_calls.annotated.vcf.gz
```

- `EXON`: binary, ≥1 bp overlap
- `REPEAT`: amount, ≥50 bp overlap
- `CENTROMERE`: binary, full overlap

# Arguments

#### Input/Output arguments:

-i, --invcf VCF
:   input VCF/BCF file

--bed FILE
:   comma-separated list of genomic feature BED files (3 columns: CHR, START, END)

--name NAME
:   comma-separated list of feature names added to VCF INFO. Must be in the same order as `--bed`.

-o, --output FILE
:   output file path. Extension determines format: `.vcf`, `.vcf.gz`, `.bcf`, or any extension with `--out-type tsv`.

--out-type STR
:   output file type: `vcf` (default) or `tsv`.

#### Overlap arguments:

--mode STR
:   comma-separated list of intersection modes for each feature, in the same order as `--bed`. `b` (default): binary flag (0/1). `a`: amount of overlapping base pairs.

--overlap STR
:   comma-separated list of minimum required overlap in base pairs for each feature, in the same order as `--bed`. `1` (default): at least 1 bp overlap. `-1`: SV must be fully contained within the feature. Values ≥1 specify the minimum overlap.
