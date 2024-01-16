# Tutorial

*** Last updated: 2024-01-16 ***

This tutorial provides a step-by-step guide on how to utilize `harmonisv` for post-processing the outcomes of various SV calling methods and conducting joint SV calling across multiple samples and methods. The necessary input data and scripts for this tutorial can be found in the [test] folder on GitHub.

## 1. SV Discovery Using Various Methods

Please note that `harmonisv` itself does not directly perform SV calling. Instead, it is designed to harmonize and integrate the results from any SV calling methods that generate a VCF file following the VCF specification and contains the essential information:

- Type of SV
- Length of SV
- Sequencing depth for both reference and alternative alleles

In this tutorial, we employ two aligners, minimap2 and NGMLR, and three SV callers, cuteSV, sniffles, and SVIM, to call SVs in the HG002 sample.

Align long-read sequencing data to reference:

```bash
ref="hs37d5.fa"
input="HG002.fasta"

# minimap2
minimap2 -L -t 36 --MD -a -x map-pb $ref $input | samtools sort -@ 4 -o $bam

# NGMLR
ngmlr --bam-fix -t 40 -x pacbio \
-r $ref \
-q $input \
-o $bam
```

Call SVs from alignments:

```bash
# sniffles (ver 2.0.6)
sniffles \
-i $bam \
-v $vcf \
--reference $ref \
--tandem-repeats "human_hs37d5.trf.bed"

# SVIM (ver 2.0.0)
svim alignment $outpath $input $ref \
--max_sv_size 1000000

# cuteSV (ver 2.0.3)
cuteSV \
--max_cluster_bias_INS 100 \
--diff_ratio_merging_INS 0.3 \
--max_cluster_bias_DEL 200 \
--diff_ratio_merging_DEL 0.5 \
--genotype \
$input $ref $vcf $workdir

```

After running the above commands, we get the following VCF files:

- HG002.minimap2.sniffles.vcf
- HG002.minimap2.svim.vcf
- HG002.minimap2.cuteSV.vcf
- HG002.NGMLR.sniffles.vcf
- HG002.NGMLR.svim.vcf
- HG002.NGMLR.cuteSV.vcf

## 2. Harmonize VCFs

The output VCF files from different SV calling methods usually have different formats in their `INFO` and `FORMAT` fields. To integrate the results from different methods, we first need to harmonize the VCFs to a standard format. This is because when we rename or move VCF tags (e.g., FORMAT/tag to INFO/tag), the newly added tags should also have definition in the header according to the VCF specification. Here, we first use `harmonize-header` to combine the headers from all input VCFs.

```bash
# If one tag is defined in multiple VCFs
# the one from the reference VCF or the first VCF in the list will be used
harmonisv harmonize-header \
-f vcf_list.txt \                # read VCF list from file
-o harmonized_header.txt \       # output file
-r HG002.minimap2.sniffles.vcf   # reference VCF file (optional)
```

Then, we can use `harmonize` to standardize the VCFs:

1. reheader with `harmonized_header.txt`
2. for the same information stored in different tags, rename them to the same tag
3. normalize SV types to the same format (e.g., `DUP:TANDEM` to `DUP`)
4. extract all required tags to `INFO` field
5. rename SV ID to make them unique across samples and methods

```bash
# sniffles
harmonisv harmonize \
-i HG002.minimap2.sniffles.vcf \
-o HG002.minimap2.sniffles.harmonized.vcf \
--info SVTYPE,SVLEN,END,STRANDS=STRAND \        # INFO tags to be kept
--format-to-info RE=DV \                        # extract FORMAT/DV to INFO/RE
--format-to-info-sum DP=DR,DP=DV \              # INFO/DP = FORMAT/DR + FORMAT/DV
--header harmonized_header.txt \                # reheader
--id-prefix HG002.minimap2.sniffles \           # add sample ID and method as prefix to SV ID
--rename-id                                     # rename SV ID to prefix.chr.svtype.number

# SVIM
harmonisv harmonize \
-i HG002.minimap2.svim.vcf \
-o HG002.minimap2.svim.harmonized.vcf \
--info SVTYPE,SVLEN,END,RE=SUPPORT \            # INFO tags to be kept
--format-to-info DP=DP \                        # extract FORMAT/DP to INFO/DP
--DUP DUP,DUP:TANDEM,DUP:INT \                  # normalize SV types
--header harmonized_header.txt \                # reheader
--id-prefix HG002.minimap2.svim \               # add sample ID and method as prefix to SV ID
--rename-id                                     # rename SV ID to prefix.chr.svtype.number

# cuteSV
harmonisv harmonize \
-i HG002.minimap2.cuteSV.vcf \
-o HG002.minimap2.cuteSV.harmonized.vcf \
--info SVTYPE,SVLEN,END,RE \                   # INFO tags to be kept
--format-to-info-sum DP=DR,DP=DV \             # INFO/DP = FORMAT/DR + FORMAT/DV
--header harmonized_header.txt \               # reheader    
--id-prefix HG002.minimap2.cuteSV \            # add sample ID and method as prefix to SV ID
--rename-id                                    # rename SV ID to prefix.chr.svtype.number
```

## 3. Merge VCFs
Ongoing...

[test]: https://github.com/Han-Cao/HarmoniSV/tree/master/test