#!/bin/bash

# test run of harmoniSV
harmonisv="../src/harmoniSV/harmonisv"
chmod +x $harmonisv

# harmonize headers from list of input vcf files
dir_header="output/harmonize_header/"
[[ ! -d $dir_header ]] && mkdir -p $dir_header

ls raw/*NGMLR* > $dir_header/NGMLR_vcf_list.txt

$harmonisv harmonize-header \
-i raw/HG002.minimap2.cuteSV.vcf,raw/HG002.minimap2.sniffles.vcf,raw/HG002.minimap2.svim.vcf \
-f $dir_header/NGMLR_vcf_list.txt \
-o $dir_header/harmonized_header.txt \
-r raw/HG002.minimap2.sniffles.vcf 2> $dir_header/harmonize_header.log

# harmonize VCFs
dir_harmonize="output/harmonize/"
[[ ! -d $dir_harmonize ]] && mkdir -p $dir_harmonize

# sniffles
ls raw/*sniffles* | while read vcf; do
    name=$(basename $vcf)
    name=${name%.vcf}

    $harmonisv harmonize \
    -i $vcf \
    -o ${dir_harmonize}/${name}.harmonized.vcf \
    --info SVTYPE,SVLEN,END,STRANDS=STRAND \
    --format-to-info RE=DV \
    --format-to-info-sum DP=DR,DP=DV \
    --header $dir_header/harmonized_header.txt \
    --id-prefix $name \
    --rename-id 2> ${dir_harmonize}/${name}.harmonized.vcf.log

done

