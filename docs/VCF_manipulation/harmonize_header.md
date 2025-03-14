# *harmonize-header*

Harmonize VCF headers

***Last updated: 2025-03-14***

## Input requirmenets
- **VCF**: VCF/BCF files following the VCF specification.


## Output
- **VCF header**: harmonized VCF header containing all header records of input VCFs


## Usage

``` bash
harmonisv harmonize-header [options] -i <input_vcf> -o <output_header> 
```

## Examples

In this example, the VCF headers from cuteSV, sniffles, and svim are harmonized. If one header occurs in more than one VCF, the priority is: sniffles > cuteSV > svim (based on the input order).

``` bash
harmonisv harmonize-header \
-i HG002.minimap2.cuteSV.vcf,HG002.minimap2.sniffles.vcf,HG002.minimap2.svim.vcf \
-o harmonized_header.txt \
-r HG002.minimap2.sniffles.vcf # reference VCF has the highest priority
```

Input:
``` bash
HG002.minimap2.cuteSV.vcf:
##fileformat=VCFv4.2
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variant">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variant">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NULL

HG002.minimap2.sniffles.vcf:
##fileformat=VCFv4.2
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Structural variation with precise breakpoints">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Structural variation with imprecise breakpoints">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variation">
##INFO=<ID=STDEV_POS,Number=1,Type=Float,Description="Standard deviation of structural variation start position">
##INFO=<ID=STDEV_LEN,Number=1,Type=Float,Description="Standard deviation of structural variation length">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE

HG002.minimap2.svim.vcf:
##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=STD_SPAN,Number=1,Type=Float,Description="Standard deviation in span of merged SV signatures">
##INFO=<ID=STD_POS,Number=1,Type=Float,Description="Standard deviation in position of merged SV signatures">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample
```

Output:

`PRECISE`, `IMPRECISE`, `SVLEN`, `SVTYPE` are defined in multiple VCFs, and the output use the definition in sniffles as it has the highest priority. `STDEV_POS`, `STDEV_LEN`, `STD_SPAN`, `STD_POS` are defined in single VCF and are appended to the output header.

``` bash
##fileformat=VCFv4.2
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Structural variation with precise breakpoints">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Structural variation with imprecise breakpoints">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variation">
##INFO=<ID=STDEV_POS,Number=1,Type=Float,Description="Standard deviation of structural variation start position">
##INFO=<ID=STDEV_LEN,Number=1,Type=Float,Description="Standard deviation of structural variation length">
##INFO=<ID=STD_SPAN,Number=1,Type=Float,Description="Standard deviation in span of merged SV signatures">
##INFO=<ID=STD_POS,Number=1,Type=Float,Description="Standard deviation in position of merged SV signatures">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
```

## Arguments

##### Input/Output arguments:
-i, --invcf *VCF*
:   Comma-separated list of input VCF files. Duplicate headers will use the first one, including SAMPLE header. For multi-sample VCF, please make sure all input VCFs have the same SAMPLE order.
  
-f, --file-list *FILE_LIST*
:   File containing a list of input VCF files, one VCF per line. VCFs from both -i and -f will be used.

-o, --output *OUTPUT*
:   Output VCF header file

##### optional arguments:
-r, --ref-vcf *VCF*
:   Reference VCF file, highest priority for duplicated headers
