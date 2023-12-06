# *represent*

Select the representative SV from merged SVs

*** Last updated: 2023-12-05 ***

## Input
- **VCF**: 
    - **Format**: bi-allelic VCF/BCF files following the VCF specification
    - **Required INFO**: `SVTYPE`, `SVLEN`, `AC`
    - **Required ID**: unique ID across all input VCFs
- **VCF file list**:
    - **Format**: one VCF per line
- **SV merging results**:
    - **Format**: each line is comma-separated list of IDs of merged SVs, or the ID of unique SVs. Any SV merging method can be used. Particularly, we have tested [Jasmine] in the workflow.

## Output
- **VCF**: harmonized VCF/BCF files.
- **ID format**: `{PREFIX}.{SVTYPE}.{NUMBER}` if `--rename-id --id-prefix PREFIX` is specified, otherwise the same as input.

## Usage

``` bash
harmonisv represent [options] -i <input_vcf> -o <output_vcf> --merge <merge_file>
harmonisv represent [options] -f <file_list> -o <output_vcf> --merge <merge_file>
```

## Representative SV
SV calling results can vary significantly across different methods or samples. Therefore, merging SVs across samples and methods is a common practice to identify non-redundant SVs. Once the SVs are merged, the position, length, and alleles of the merged SV can be determined based on the representative SV. The `represent` command is used to select the representative SV from the merged SVs using one of the following criteria:

- `--by-max TAG`: This option selects the SV with the maximum value of `INFO/TAG`.
- `--by-freq`: This option selects the SV with the maximum frequency of `POS` and `SVLEN`. If more than one SV has the same frequency, it selects the one closest to the average `POS` and `SVLEN`.

## Example

##### 1. Select representative SVs from SVs merged across samples and methods
In this example, we will first use [Jasmine] to merge SVs across samples and methods, and then select the representative SV based on the prevalence of SVs' positions and lengths.

``` bash

# prepare VCF file list
for vcf in SampleA.Method1.vcf SampleA.Method2.vcf SampleB.Method1.vcf SampleB.Method2.vcf; do
    echo ${vcf}
done > vcf_list.txt

# SV merging across samples and methods
jasmine \
file_list=vcf_list.txt \
out_file=merged.vcf \
--keep_var_ids \
--ignore_strand 

# extract SV merging results
bcftools query -f '%IDLIST\n' merged.vcf > merge.txt

# Select representative SV based on the frequency of POS and SVLEN
harmonisv represent \
-f vcf_list.txt \
-o representative.vcf \
--merge merge.txt \
--by-freq \
--id-prefix Merged \            # rename output SV ID as Merged.SVTYPE.NUMBER
--save-id  \                    # save merged SVs to INFO/ID_LIST
--min-len-input 35 \            # only SVs with SVLEN < 35 are considered when selecting representative SV
--min-len-output 50             # drop representative SVs with SVLEN < 50
```


##### 2. Remove duplicated SVs within the same sample and method

SV calling methods may produce duplicated SVs within the same sequencing data. This example will first use [Jasmine] to merge SVs within the same sample and method, and then select the representative SV based on the maximum value of reads supporting the SV.


``` bash
# SV merging within the same sample and method
jasmine \
file_list=input.vcf \
out_file=jasmine.vcf \
--comma_filelist \
max_dist=200 \
--allow_intrasample \
--nonlinear_dist \
--ignore_strand \
--keep_var_ids

# extract SV merging results
bcftools query -f '%INTRASAMPLE_IDLIST\n' jasmine.vcf > dup_sv.txt

# Select representative SV based on the maximum value of reads supporting the SV (INFO/RE)
harmonisv represent \
-i input.vcf \
-o dedup.vcf \
--merge dup_sv.txt \
--by-max RE \
--min-len-input 30 \
--min-len-output 30
```

## Arguments

#### Input/Output arguments:
-i, --invcf VCF
:   input VCF

-f, --file-list TSV
:   list of VCF files, one VCF per line

-o, --outvcf VCF
:   output VCF

--merge FILE
:   SV merging results, each line is a comma-separated list of IDs of merged SVs, or the ID of unique SVs

-r, --region CHR
:   genomic region to work on (requires VCF with index)

#### Representative SV selection arguments:
--id-prefix PREFIX
:   Rename output SV ID as PREFIX.SVTYPE.NUMBER

--by-max TAG
:   Select representative SV by maximum INFO/TAG value

--by-freq
:   Select representative SV by frequency of POS, SVLEN, if >1 SVs have same frequency, select the one closest to average POS, SVLEN

--save-id
:   Save original IDs to INFO/ID_LIST

--keep-hom-ref
:   Keep SVs with AC=0 before selecting representative SV (Default: False)

--min-len-input N
:   Remove SVLEN < min_len_input before selecting representative SV (Default: 35)

--min-len-output N
:   Remove SVLEN < min_len_output from output (Default: 50)

[Jasmine]: https://github.com/mkirsche/Jasmine