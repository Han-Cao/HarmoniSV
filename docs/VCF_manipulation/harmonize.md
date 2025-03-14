# *harmonize*

Harmonize SV VCFs across samples and SV calling methods

***Last updated: 2025-03-14***

## Input requirmenets
- **VCF**: 
    - **Format**: bi-allelic VCF/BCF files following the VCF specification
    - **INFO**: `SVTYPE`, if other name is used, specify by `--svtype`

## Output
- **VCF**: harmonized VCF/BCF files.
    - **ID**: `{PREFIX}.{SVTYPE}.{NUMBER}` if `--rename-id --id-prefix PREFIX` is specified, otherwise the same as input.

## Usage

``` bash
harmonisv harmonize [options] -i <input_vcf> -o <output_vcf> 
```


## Input filtering

By default, the following records will be removed:

- Records with missing ID or duplicated ID if `--rename-id` is not specified
- Records without `INFO/SVTYPE`
- Records' `INFO/SVTYPE` not in "INS, DEL, DUP, INV, CNV, BND"

If the types of SVs are stored in `INFO` field with different names, use `--svtype` to specify the tag name. 

If the SV types are not in the above list, use `--INS`, `--DEL`, `--DUP`, `--INV`, `--CNV`, `--BND` to normalize the SV types. For example, `--info SVTYPE_raw=SVTYPE --DUP DUP,DUP:TANDEM` will normalize `SVTYPE=DUP:TANDEM` to `SVTYPE=DUP;SVTYPE_raw=DUP:TANDEM`.

To disable the check, use `--no-check`.

## VCF headers
The VCF headers produced by different SV calling methods are usually not compatable. It is necessary to have a unified VCF header for downstream analysis. If only INFO tags from the original VCFs will be used, one simple method is to use `harmonisv harmonize-header` to combine the headers from all VCFs (see [harmonize-header]). However, if any new INFO tags will be added / renamed, please use `--header-str` to specify them.

For example, the following command will combine the headers of method A, B, and C, and add `INFO/STRANDS` header.

``` bash
harmonisv harmonize-header \
-i A.vcf,B.vcf,C.vcf \
-o harmonized_header.txt \

for vcf in A.vcf B.vcf C.vcf; do
    harmonisv harmonize \
    -i $vcf \
    -o harmonized.$vcf \
    --header harmonized_header.txt \
    --header-str 'STRANDS,1,String,Strand orientation of supporting reads'
done
```


## Examples

In this example, we will rename variant ID, normalize SVTYPE, and extract the basic SV calling information:

- SVTYPE
- SVLEN
- END
- STRANDS (if exists)
- DP (sequencing depth)
- RE (number of reads supporting the SV)

Please note that SV-calling methods may store this information under different names. For instance, `INFO/SUPPORT` might be used for `RE`, and `FORMAT/DP` for `DP`. The sequencing depth might also be provided separately for REF and ALT alleles (e.g., `FORMAT/DR` and `FORMAT/DV`). Additionally, the same type of SV might have different names in different methods (e.g., `DUP:TANDEM` vs `DUP`). All the harmonized information will be stored in the `INFO` field.

** Sniffles2 (ver 2.0.6) **

``` bash
harmonisv harmonize \
-i HG002.minimap2.sniffles.vcf \             # input VCF
-o HG002.minimap2.sniffles.harmonized.vcf \  # output VCF
--info SVTYPE,SVLEN,END,STRANDS=STRAND \     # INFO fields to be kept, rename STRAND to STRANDS
--format-to-info RE=DV \                     # Extract FORMAT/DV to INFO/RE
--format-to-info-sum DP=DR,DP=DV \           # Calculate INFO/DP = FORMAT/DR + FORMAT/DV
--header harmonized_header.txt \             # Replace VCF header
--id-prefix HG002.minimap2.sniffles \        # {SAMPLE}.{ALIGNER}.{CALLER}
--rename-id                                  # Rename all variant ID
```

** SVIM (ver 2.0.0) **

``` bash
harmonisv harmonize \                         
-i HG002.minimap2.svim.vcf \                 # input VCF
-o HG002.minimap2.svim.harmonized.vcf \      # output VCF
--info SVTYPE,SVLEN,END,RE=SUPPORT \         # INFO fields to be kept, rename SUPPORT to RE
--format-to-info DP=DP \                     # Extract FORMAT/DP to INFO/DP
--DUP DUP,DUP:TANDEM,DUP:INT \               # Normalize the name of duplications
--header harmonized_header.txt \             # Harmonized VCF header
--id-prefix HG002.minimap2.svim \            # {SAMPLE}.{ALIGNER}.{CALLER}
--rename-id                                  # Rename all variant ID
```

** cuteSV (ver 2.0.3) **

``` bash
harmonisv harmonize \
-i HG002.minimap2.cuteSV.vcf \               # input VCF
-o HG002.minimap2.cuteSV.harmonized.vcf \    # output VCF
--info SVTYPE,SVLEN,END,RE \                 # INFO fields to be kept
--format-to-info-sum DP=DR,DP=DV \           # Calculate INFO/DP = FORMAT/DR + FORMAT/DV
--header harmonized_header.txt \             # Harmonized VCF header
--id-prefix HG002.minimap2.cuteSV \          # {SAMPLE}.{ALIGNER}.{CALLER}
--rename-id                                  # Rename all variant ID
```

By default, all INFO tags not specified in `--info` will be removed. To keep all original INFO tags, use `--keep-all`. To keep the original INFO tags before renaming, use `--keep-old`. Meanwhile, allele count (AC) and allele number (AN) will also be automatically computed from genotypes. To disable this feature, use `--no-AC`.

Particularly, `harmonisv harmonize -i input.vcf -o output.vcf --keep-all --no-AC --no-check` should produce the same variant records as the input (some VCF headers are still added).


## Arguments

#### Input/Output arguments:
-i, --invcf VCF
:   input VCF

-o, --outvcf VCF
:   output VCF

#### VCF INFO manipulation:
--info TAG
:   Comma separated INFO tags to extract or rename. INFO tags can be renamed by NEW=OLD, from high prioirty to low priority, e.g., 'NEW=OLD1,NEW=OLD2' means if OLD1 is present, use OLD1, otherwise use OLD2.

--info-sum TAG
:   Comma separated INFO tags to extract and sum. Old tags with the same new tag name will be summed up, e.g., 'NEW=OLD1,NEW=OLD2' means 'INFO/NEW = INFO/OLD1 + INFO/OLD2'. Please define the header of new tags in --header or --header-str

--format-to-info TAG
:   Comma separated FORMAT tags to sum across samples and add to INFO, from high prioirty to low priority, e.g., DP=DP means 'INFO/DP = sum(FORMAT/DP)'.

--format-to-info-sum TAG
:   Comma separated FORMAT tags to sum across samples and tags and add to INFO, e.g., 'DP=DR,DP=DV' means 'INFO/DP = sum(FORMAT/DR) + sum(FORMAT/DV)'. Please define the header of new tags in --header or --header-str

--info-to-alt TAG
:   Comma separated INFO tags to fill in ALT, from high prioirty to low priority. This is useful if insertion sequence stored in INFO.

--keep-old
:   Keep original INFO tags after renaming or sum (default: False)

--keep-all
:   Keep all original INFO tags (default: False)

--no-AC
:   Disable automatically add AC and AN to INFO (default: False)

#### VCF header manipulation:
--header FILE
:   New VCF header file to replace the header of input VCF.

--header-str string
:   Semicolon separated INFO header string added to new header (metadata separated by comma), e.g., 'DP,1,Integar,Sequencing depth;AF,1,Float,Allele frequency'

#### Structural variation format:
--id-prefix PREFIX
:   Rename SV ID to PREFIX.raw_ID. Final ID should be Sample.Aligner.Caller.unique_id for downstream analysis
  
--rename-id
:   Rename SV ID to PREFIX.SVTYPE.No., must use with --id-prefix
  
--svtype SVTYPE
:   INFO tag stores the structural variation type (default: SVTYPE), will rename it to SVTYPE if not.
  
--INS INS
:   Comma separated SVTYPE string for insertions, will be nomalized as INS (default: INS)
  
--DEL DEL
:   Comma separated SVTYPE string for deletions, will be nomalized as DEL (default: DEL)
  
--DUP DUP
:   Comma separated SVTYPE string for duplications, will be nomalized as DUP (default: DUP)

--INV INV
:   Comma separated SVTYPE string for inversions, will be nomalized as INV (default: INV)
  
--CNV CNV
:   Comma separated SVTYPE string for copy number variations, will be nomalized as CNV (default: CNV)

--BND BND
:   Comma separated SVTYPE string for breakends, will be nomalized as BND (default: BND)

--no-check
:   Disable check and filter on variant ID and SVTYPE (default: False)



[harmonize-header]: harmonize_header.md