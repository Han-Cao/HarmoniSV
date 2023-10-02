#!/bin/bash

# test run of harmoniSV
harmonisv="../harmoniSV/harmonisv"

# harmonize headers from list of input vcf files
dir_header="output/harmonize_header"
[[ ! -d $dir_header ]] && mkdir -p $dir_header

ls raw/*NGMLR* > $dir_header/NGMLR_vcf_list.txt

$harmonisv harmonize-header \
-i raw/HG002.minimap2.cuteSV.vcf,raw/HG002.minimap2.sniffles.vcf,raw/HG002.minimap2.svim.vcf \
-f $dir_header/NGMLR_vcf_list.txt \
-o $dir_header/harmonized_header.txt \
-r raw/HG002.minimap2.sniffles.vcf 2> $dir_header/harmonize_header.log