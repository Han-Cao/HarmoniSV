# *harmonize*

Harmonize SV VCFs across samples and SV calling methods

*** Last updated: 2023-03-19 ***

## Input
- **VCF**: any VCF/BCF files following the VCF specification.
- **ID format**: no requirements

## Output
- **VCF**: VCF/BCF files following the VCF specification.
- **ID format**: `{PREFIX}.{SVTYPE}.{NUMBER}` if `--rename-id --id-prefix PREFIX` is specified, otherwise the same as input.
- 

## Usage

``` bash
harmonisv harmonize [options] -i <input_vcf> -o <output_vcf> 
```

## Examples

###### Example 1: Harmonize sniffles2 (ver 2.0.6) output

Output VCF of `Sniffles2` store SV information (e.g., SVTYPE, SVLEN) in the INFO field, and the read depth of reference (DR) and variant (DV) alleles in the FORMAT field. The following command keep basic SV information and extract read depths to INFO fields.

``` bash
harmonisv harmonize \
-i input.vcf.gz \                           # input VCF
-o output.vcf.gz \                          # output VCF
--info SVTYPE,SVLEN,END,STRANDS=STRAND \    # INFO fields to be kept, rename STRAND to STRANDS
--format-to-info DP=DR,DP=DV,RE=DV \        # Extract FORMAT fields to INFO fields
--sum \                                     # Sum all values assigned to the same key (DP = DR + DV)
--header header.txt \                       # New VCF header, please define new keys here (i.e., DP, RE)
--rename-id \                               # rename ID to {PREFIX}.{SVTYPE}.{NUMBER}
--id-prefix HG002.minimap2.sniffles2        # {SAMPLE}.{ALIGNER}.{CALLER}
```

## Arguments

##### Required arguments

-i, --invcf *vcf*   
:   input vcf

-o, --outvcf *vcf*  
:   output vcf

##### Optional arguments
--info TAG
:   Comma separated INFO tags to extract, can rename tag by NEW=OLD. Can give multiple candidate old tags for 1 new tag, priority from high to low. E.g., 'NEW=OLD1,NEW=OLD2' means if OLD1 is present, use OLD1, otherwise use OLD2.

--info-to-alt TAG
:   Comma separated INFO tags to fill in ALT, from high prioirty to low priority. This is useful for insertion sequence stored in INFO.

--format-to-info TAG
:   Comma separated FORMAT tags to be sum and add in INFO, from high prioirty to low priority. New headers must given by --header.

--sum
:   Change merge logic to sum, all tags must exist for all records. E.g., '--format-to-info DP=DR,DP=DV --sum' means 'INFO/DP = sum(FORMAT/DR + FORMAT/DV)'

--header FILE
:   New vcf header file to replace the header of invcf

--header-str string
:   Semicolon separated INFO header string added to new header (metadata separated by comma). e.g., 'DP,1,Integar,Sequencing depth;AF,1,Float,Allele frequency'

--id-prefix PREFIX
:   Rename SV ID to PREFIX.raw_ID. Final ID should be Sample.Aligner.Caller.unique_id for downstream analysis

--rename-id
:   Rename SV ID to PREFIX.SVTYPE.No., must use with --id-prefix

--keep-old
:   Keep raw INFO when rename, useful when merging (default: False)

--keep-all
:   Whether keep all INFO fields. If specified, only rename work in --info (default: False)

--no-AC
:   Disable automatically add AC to INFO (default: False)

--no-check
:   Disable vcf check and filter on ID, SVTYPE for downstream analysis (default: False)
