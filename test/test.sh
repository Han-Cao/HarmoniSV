#!/bin/bash

# test run of harmoniSV
harmonisv="../src/harmoniSV/harmonisv"
chmod +x $harmonisv

# harmonize headers from list of input vcf files
echo "################## 1. harmonize-header ##################"
dir_header="output/harmonize_header/"
[[ ! -d $dir_header ]] && mkdir -p $dir_header

ls raw/*NGMLR* > $dir_header/NGMLR_vcf_list.txt

$harmonisv harmonize-header \
-i raw/HG002.minimap2.cuteSV.vcf,raw/HG002.minimap2.sniffles.vcf,raw/HG002.minimap2.svim.vcf \
-f $dir_header/NGMLR_vcf_list.txt \
-o $dir_header/harmonized_header.txt \
-r raw/HG002.minimap2.sniffles.vcf

# harmonize VCFs
echo "################## 2. harmonize VCFs ##################"
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
    --header-str 'STRANDS,1,String,Strand orientation of supporting reads' \
    --id-prefix $name \
    --rename-id

done

# SVIM
ls raw/*svim* | while read vcf; do
    name=$(basename $vcf)
    name=${name%.vcf}

    $harmonisv harmonize \
    -i $vcf \
    -o ${dir_harmonize}/${name}.harmonized.vcf \
    --info SVTYPE,SVLEN,END,RE=SUPPORT \
    --format-to-info DP=DP \
    --DUP DUP,DUP:TANDEM,DUP:INT \
    --header $dir_header/harmonized_header.txt \
    --header-str 'STRANDS,1,String,Strand orientation of supporting reads' \
    --id-prefix $name \
    --rename-id 2> ${dir_harmonize}/${name}.harmonized.vcf.log
done

# cuteSV
ls raw/*cuteSV* | while read vcf; do
    name=$(basename $vcf)
    name=${name%.vcf}

    $harmonisv harmonize \
    -i $vcf \
    -o ${dir_harmonize}/${name}.harmonized.vcf \
    --info SVTYPE,SVLEN,END,RE \
    --format-to-info-sum DP=DR,DP=DV \
    --header $dir_header/harmonized_header.txt \
    --header-str 'STRANDS,1,String,Strand orientation of supporting reads' \
    --id-prefix $name \
    --rename-id
done

# Jasmine code to identify duplicated SV calls
dir_dupcall="dup_call/"
[[ ! -d $dir_dupcall ]] && mkdir -p $dir_dupcall

# ls output/harmonize/*harmonized.vcf | while read vcf; do
#     name=$(basename $vcf)
#     name=${name%.vcf}

#     jasmine file_list=$vcf \
#     out_file=${dir_dupcall}/${name}.dup_call.vcf \
#     genome_file=hs37d5.fa \
#     --comma_filelist \
#     max_dist=200 \
#     --allow_intrasample \
#     --nonlinear_dist \
#     --ignore_strand \
#     --keep_var_ids 

#     # write merged SV list
#     bcftools query -f '%INTRASAMPLE_IDLIST\n' ${dir_dupcall}/${name}.dup_call.vcf \
#     > ${dir_dupcall}/${name}.dup_call.txt
# done

# remove duplicated SVs
# i.e., find representative SVs among duplicated SVs based on RE
echo "################## 3. Remove duplicated SVs within the same method ##################"
dir_dedup="output/dedup/"
[[ ! -d dir_dedup ]] && mkdir -p $dir_dedup

ls ${dir_harmonize}/*harmonized.vcf | while read vcf; do
    name=$(basename $vcf)
    name=${name%.vcf}

    $harmonisv represent \
    -i $vcf \
    -o ${dir_dedup}/${name}.dedup.vcf \
    --merge ${dir_dupcall}/${name}.dup_call.txt \
    --by-max RE \
    --min-len-input 30 \
    --min-len-output 30
done

# Jasmine code to merge SVs from different methods
dir_merge="sv_merge/"
[[ ! -d $dir_merge ]] && mkdir -p $dir_merge

# ls ${dir_harmonize}/*harmonized.vcf > ${dir_merge}/All_method.merge_vcf_list.txt

# jasmine \
# file_list=${dir_merge}/All_method.merge_vcf_list.txt \
# out_file=${dir_merge}/All_method.merge.vcf \
# --keep_var_ids \
# --ignore_strand 

# bcftools query -f '%IDLIST\n' ${dir_merge}/All_method.merge.vcf > ${dir_merge}/All_method.merge.txt

# Find representative SVs among merged SVs based on frequency / distance
echo "################## 4. Find representative SVs among different methods ##################"
dir_represent="output/represent/"
[[ ! -d $dir_represent ]] && mkdir -p $dir_represent

$harmonisv represent \
-f ${dir_merge}/All_method.merge_vcf_list.txt \
-o ${dir_represent}/All_method.representative.vcf \
--merge ${dir_merge}/All_method.merge.txt \
--by-freq \
--id-prefix All_method \
--save-id
