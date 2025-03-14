# Tutorial

***Last updated: 2025-03-14***

This tutorial provides a step-by-step guide on how to utilize `harmonisv` for post-processing the results of various SV calling methods and conducting joint SV calling across multiple samples and methods. The necessary input data and scripts for this tutorial can be found in the [test] folder on GitHub.

In this tutorial, we use the following tools to showcase the usage of `harmonisv`:

- aligners: `minimap2` and `NGMLR`
- SV callers: `cuteSV`, `sniffles`, and `SVIM`
- SV merging: `jasmine`

## Table of Contents

- [1. SV discovery using multiple methods](#1-sv-discovery-using-mutliple-methods)
- [2. Harmonize VCFs](#2-harmonize-vcfs)
- [3. Remove Duplicated SVs](#3-remove-duplicated-svs)
- [4. SV merging](#4-sv-merging)
- [5. SV re-genotyping](#5-sv-re-genotyping)
- [6. Combine per-sample SV genotyping results](#6-combine-per-sample-sv-genotyping-results)

## 1. SV discovery using mutliple methods

`harmonisv` is designed to harmonize and integrate the results from different SV calling methods. It requires the input VCF files include the following information:

- Type of SV
- Length of SV
- Sequencing depth for both reference and alternative alleles

We have generated the example VCF files of HG002 chr22 in the [test/raw] folder:

- HG002.minimap2.sniffles.vcf
- HG002.minimap2.svim.vcf
- HG002.minimap2.cuteSV.vcf
- HG002.NGMLR.sniffles.vcf
- HG002.NGMLR.svim.vcf
- HG002.NGMLR.cuteSV.vcf

These VCF files were generated using the following commands:

```bash
ref="hs37d5.fa"
input="HG002.fasta"

# align long-read to reference
# minimap2
minimap2 -L -t 36 --MD -a -x map-pb $ref $input | samtools sort -@ 4 -o $bam

# NGMLR
ngmlr --bam-fix -t 40 -x pacbio \
-r $ref \
-q $input \
-o $bam

# SV calling
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

## 2. Harmonize VCFs

The output VCF files from different SV calling methods usually have different formats in their `INFO` and `FORMAT` fields. To integrate the results from different methods, we need to harmonize the VCFs to a standard format. Here, we first use `harmonize-header` to combine the headers from all input VCFs.

```bash
# If one tag is defined in multiple VCFs
# the one from the reference VCF or the first VCF in the list will be used
dir_header="output/harmonize_header/"
[[ ! -d $dir_header ]] && mkdir -p $dir_header

ls raw/*NGMLR* > $dir_header/NGMLR_vcf_list.txt

harmonisv harmonize-header \
-i raw/HG002.minimap2.cuteSV.vcf,raw/HG002.minimap2.sniffles.vcf,raw/HG002.minimap2.svim.vcf \
-f $dir_header/NGMLR_vcf_list.txt \
-o $dir_header/harmonized_header.txt \
-r raw/HG002.minimap2.sniffles.vcf
```

Then, we can use `harmonize` to standardize the VCFs:

1. standardize VCF headers and tag names
2. normalize SV types (e.g., `DUP:TANDEM` to `DUP`)
3. extract depth and number of supporting reads to `INFO/DP` and `INFO/RE`, respectively
4. rename SV ID to make them unique across samples and methods. By using `--id-prefix`` and `--rename-id`, a unique ID in the format of `prefix.chr.svtype.number` will be assigned to each SV.

```bash
dir_harmonize="output/harmonize/"
[[ ! -d $dir_harmonize ]] && mkdir -p $dir_harmonize

# sniffles
# note: INFO/DP = FORMAT/DR + FORMAT/DV
ls raw/*sniffles* | while read vcf; do
    name=$(basename $vcf)
    name=${name%.vcf}

    harmonisv harmonize \
    -i $vcf \
    -o ${dir_harmonize}/${name}.harmonized.vcf \
    --info SVTYPE,SVLEN,END,STRANDS=STRAND \
    --format-to-info RE=DV \
    --format-to-info-sum DP=DR,DP=DV \
    --header $dir_header/harmonized_header.txt \
    --header-str 'STRANDS,1,String,Strand orientation of supporting reads' \
    --id-prefix $name \
    --rename-id

done

# SVIM
# note: convert SVTYPE "DUP:TANDEM" and "DUP:INT" to "DUP"
ls raw/*svim* | while read vcf; do
    name=$(basename $vcf)
    name=${name%.vcf}

    harmonisv harmonize \
    -i $vcf \
    -o ${dir_harmonize}/${name}.harmonized.vcf \
    --info SVTYPE,SVLEN,END,RE=SUPPORT \
    --format-to-info DP=DP \
    --DUP DUP,DUP:TANDEM,DUP:INT \
    --header $dir_header/harmonized_header.txt \
    --header-str 'STRANDS,1,String,Strand orientation of supporting reads' \
    --id-prefix $name \
    --rename-id
done

# cuteSV
ls raw/*cuteSV* | while read vcf; do
    name=$(basename $vcf)
    name=${name%.vcf}

    harmonisv harmonize \
    -i $vcf \
    -o ${dir_harmonize}/${name}.harmonized.vcf \
    --info SVTYPE,SVLEN,END,RE \
    --format-to-info-sum DP=DR,DP=DV \
    --header $dir_header/harmonized_header.txt \
    --header-str 'STRANDS,1,String,Strand orientation of supporting reads' \
    --id-prefix $name \
    --rename-id
done
```

## 3. Remove duplicated SVs
SV calling methods may redundantly call the same SV. We can use `jasmine` to identify intra-sample duplicated SVs. The results can be found in [test/dup_call]

```bash
dir_dupcall="dup_call/"
[[ ! -d $dir_dupcall ]] && mkdir -p $dir_dupcall

ls output/harmonize/*harmonized.vcf | while read vcf; do
    name=$(basename $vcf)
    name=${name%.vcf}

    jasmine file_list=$vcf \
    out_file=${dir_dupcall}/${name}.dup_call.vcf \
    genome_file=hs37d5.fa \
    --comma_filelist \
    max_dist=200 \
    --allow_intrasample \
    --nonlinear_dist \
    --ignore_strand \
    --keep_var_ids 

    # write duplicated SV list
    bcftools query -f '%INTRASAMPLE_IDLIST\n' ${dir_dupcall}/${name}.dup_call.vcf \
    > ${dir_dupcall}/${name}.dup_call.txt
done
```

Then, we use `represent` to remove duplicated SVs by keeping the one with more supporting reads.

```bash
dir_dedup="output/dedup/"
[[ ! -d dir_dedup ]] && mkdir -p $dir_dedup

ls ${dir_harmonize}/*harmonized.vcf | while read vcf; do
    name=$(basename $vcf)
    name=${name%.vcf}

    harmonisv represent \
    -i $vcf \
    -o ${dir_dedup}/${name}.dedup.vcf \
    --merge ${dir_dupcall}/${name}.dup_call.txt \
    --by-max RE \
    --min-len-input 30 \
    --min-len-output 30
done
```

## 4. SV merging

After generating individual SV discovery call set, we next identify non-redundant SVs across samples and methods by using `jasmine`. The results can be found in the [test/sv_merge] folder.
```bash
dir_merge="sv_merge/"
[[ ! -d $dir_merge ]] && mkdir -p $dir_merge

ls ${dir_harmonize}/*harmonized.vcf > ${dir_merge}/All_method.merge_vcf_list.txt

jasmine \
file_list=${dir_merge}/All_method.merge_vcf_list.txt \
out_file=${dir_merge}/All_method.merge.vcf \
--keep_var_ids \
--ignore_strand 

bcftools query -f '%IDLIST\n' ${dir_merge}/All_method.merge.vcf > ${dir_merge}/All_method.merge.txt
```

To select representative SVs among redundant SV calls, we use the `represent` command to keep the one with most frequent `POS` and `SVLEN`.
```bash
dir_represent="output/represent/"
[[ ! -d $dir_represent ]] && mkdir -p $dir_represent

harmonisv represent \
-f ${dir_merge}/All_method.merge_vcf_list.txt \
-o ${dir_represent}/All_method.representative.vcf \
--merge ${dir_merge}/All_method.merge.txt \
--by-freq \
--id-prefix All_method \
--save-id
```

## 5. SV re-genotyping
The VCF `All_method.representative.vcf` include all non-redundant SVs discovered from different methods and samples (if more than one). To get a fully genotyped VCF for each sample, we re-genotype those SVs using `sniffles` and `cuteSV`. The results can be found in the [test/force_call] folder.

```bash
dir_force_call="force_call/"
[[ ! -d $dir_force_call ]] && mkdir -p $dir_force_call

# 1. re-genotyping
# sniffles
sniffles \
-i $bam \
-v "${dir_force_call}/HG002.minimap2.sniffles.vcf" \
--reference "hs37d5.fa" \
-t 40 \
--tandem-repeats "human_hs37d5.trf.bed" \
--genotype-vcf ${dir_represent}/All_method.representative.vcf

# cuteSV
cuteSV \
--max_cluster_bias_INS 100 \
--diff_ratio_merging_INS 0.3 \
--max_cluster_bias_DEL 200 \
--diff_ratio_merging_DEL 0.5 \
--genotype \
-t 40 \
-L -1 \
-Ivcf ${dir_represent}/All_method.representative.vcf \
$bam hs37d5.fa "${dir_force_call}/HG002.minimap2.cuteSV.vcf" $workdir
```

After re-genotyping, we need to harmonize them:

```bash
# 2. harmonize 
dir_force_call_harmonize="output/harmonize_force_call/"
# sniffles
ls ${dir_force_call}/*sniffles* | while read vcf; do
    name=$(basename $vcf)
    name=${name%.vcf}

    harmonisv harmonize \
    -i $vcf \
    -o ${dir_force_call_harmonize}/${name}.harmonized.vcf \
    --info SVTYPE,SVLEN,END \
    --format-to-info RE=DV \
    --format-to-info-sum DP=DR,DP=DV \
    --header $dir_header/harmonized_header.txt \
    --header-str 'STRANDS,1,String,Strand orientation of supporting reads'
done

# cuteSV
ls ${dir_force_call}/*cuteSV* | while read vcf; do
    name=$(basename $vcf)
    name=${name%.vcf}

    harmonisv harmonize \
    -i $vcf \
    -o ${dir_force_call_harmonize}/${name}.harmonized.vcf \
    --format-to-info-sum DP=DR,DP=DV \
    --header $dir_header/harmonized_header.txt \
    --header-str 'STRANDS,1,String,Strand orientation of supporting reads' \
    --keep-all
done
```

## 6. Combine per-sample SV genotyping results

Now, we get 10 VCFs per sample, including 4 discovery VCFs and 6 re-genotyping VCFs. We can use the `genotype` command to merge the results among all VCFs. This will generate the final per-sample VCF for downstream analysis.

```bash
dir_genotype="output/genotype/"
[[ ! -d $dir_genotype ]] && mkdir -p $dir_genotype

$harmonisv genotype \
-i ${dir_represent}/All_method.representative.vcf \
-f manifest.txt \
-o ${dir_genotype}/HG002.representative.genotyped.vcf \
--sample HG002
```

Particularly, the `manifest.txt` is a tab-separated file with the following columns:

- `file`: path to SV discovery and force calling VCF/BCF file
- `sample`: sample ID
- `aligner`: aligner used to generate the SV call set
- `caller`: SV caller used to generate the SV call set
- `if_force_call`: whether the SV call set is generated by force-calling (1: True, 0: False)
- `info` (optional): INFO tags in additional to `DP`, `RE` to be added to the output VCF

| file                                                                | sample | aligner  | caller   | is_force_call | info                |
|---------------------------------------------------------------------|--------|----------|----------|---------------|---------------------|
| output/harmonize_force_call/HG002.NGMLR.cuteSV.harmonized.vcf       | HG002  | NGMLR    | cuteSV   | 1             | PRECISE,CIPOS,CILEN |
| output/harmonize_force_call/HG002.NGMLR.sniffles.harmonized.vcf     | HG002  | NGMLR    | sniffles | 1             |                     |
| output/harmonize_force_call/HG002.minimap2.cuteSV.harmonized.vcf    | HG002  | minimap2 | cuteSV   | 1             | PRECISE,CIPOS,CILEN |
| output/harmonize_force_call/HG002.minimap2.sniffles.harmonized.vcf  | HG002  | minimap2 | sniffles | 1             |                     |
| output/harmonize/HG002.NGMLR.cuteSV.harmonized.vcf                  | HG002  | NGMLR    | cuteSV   | 0             |                     |
| output/harmonize/HG002.NGMLR.sniffles.harmonized.vcf                | HG002  | NGMLR    | sniffles | 0             |                     |
| output/harmonize/HG002.NGMLR.svim.harmonized.vcf                    | HG002  | NGMLR    | svim     | 0             |                     |
| output/harmonize/HG002.minimap2.cuteSV.harmonized.vcf               | HG002  | minimap2 | cuteSV   | 0             |                     |
| output/harmonize/HG002.minimap2.sniffles.harmonized.vcf             | HG002  | minimap2 | sniffles | 0             |                     |
| output/harmonize/HG002.minimap2.svim.harmonized.vcf                 | HG002  | minimap2 | svim     | 0             |                     |



[test]: https://github.com/Han-Cao/HarmoniSV/tree/master/test
[test/raw]: https://github.com/Han-Cao/HarmoniSV/tree/master/test/raw
[test/dup_call]: https://github.com/Han-Cao/HarmoniSV/tree/master/test/dup_call
[test/sv_merge]: https://github.com/Han-Cao/HarmoniSV/tree/master/test/sv_merge
[test/force_call]: https://github.com/Han-Cao/HarmoniSV/tree/master/test/force_call
