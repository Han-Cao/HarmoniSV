#!/usr/bin/env python3

# Intersect SVs with genomic features
# Last update: 25-Nov-2024
# Author: Han Cao


import argparse
import logging
import pysam
import pysam.bcftools
import pandas as pd
import pyranges as pr


from utils import read_vcf, vcf_to_df, parse_cmdargs

# parse arguments
parser = argparse.ArgumentParser(prog="harmonisv intersect",
                                 description="Intersect SVs with genomic features",
                                 add_help=False)
parser.add_argument("-i", "--invcf", metavar="vcf", type=str, required=True,
                    help="input vcf file")
parser.add_argument("-o", "--output", metavar="vcf", type=str, required=True,
                    help="ouput file")
parser.add_argument("--out-type", metavar="vcf", type=str, required=False, default="vcf",
                    help="output file type, vcf or bed")
parser.add_argument("--bed", metavar="bed", type=str, required=True,
                    help="comma-separated genomic feature bed files with 3 columns: CHR, START, END")
parser.add_argument("--name", metavar="name", type=str, required=True,
                    help="comma-separated feature names added to VCF INFO, should be the same order as feature files")
parser.add_argument("--mode", metavar="b", type=str, required=False, default='b',
                    help="comma-separated intersecting modes for each feature (a: amount of overlapped base pair, b: binary), should be the same order as feature files (Default: b)")
parser.add_argument("--overlap", metavar="1", type=str, required=False, default='1',
                    help="comma-separated minimum required overlap bp for each feature (-1 indicate SV reside within feature), should be the same order as feature files (Default: 1)")

optional_arg = parser.add_argument_group('optional arguments')
optional_arg.add_argument("-h", "--help", action='help', help='show this help message and exit')


def vcf_to_pyranges(vcf: pysam.VariantFile) -> pr.PyRanges:
    """
    Convert vcf to pyranges

    Input vcf must have SVTYPE and SVLEN INFO tags

    Output pyranges will contain columns: ID, Chromosome, Start, End
    """
    rows = []
    for variant in vcf.fetch():
        new_row = {}
        new_row['ID'] = variant.id
        new_row['Chromosome'] = variant.chrom
        new_row['Start'] = variant.start
        if variant.info['SVTYPE'] == 'INS':
            new_row['End'] = new_row['Start'] + 1
        else:
            new_row['End'] = new_row['Start'] + abs(variant.info['SVLEN'])
        rows.append(new_row)
    
    df = pd.DataFrame(rows)
    prg = pr.PyRanges(df)
    
    return prg

class Feature_range():
    """Parse and intersect feature bed file"""
    def __init__(self, bed: str, name: str, mode: str, overlap: int) -> None:
        self.bed = pr.read_bed(bed)
        self.name = name
        self.mode = mode
        self.overlap = overlap

        self.logger = logging.getLogger(__name__)


    def intersect(self, ranges: pr.PyRanges) -> pd.DataFrame:
        """
        Intersect SV with features

        Return: dataframe of overlapped ranges
        """
        result = ranges.join(self.bed, strandedness=False, report_overlap=True, suffix=f"_{self.name}", apply_strand_suffix=False).as_df()
        if self.overlap > 0:
            result = result[result['Overlap'] >= self.overlap]
        elif self.overlap == -1:
            result = result[result['Overlap'] == (result['End'] - result['Start'])]
        else:
            self.logger.error(f"Invalid overlap value {self.overlap} of feature {self.name}")
            raise SystemExit()
        
        result = result.rename(columns={'Overlap': f"Overlap_{self.name}"})

        return result


def parse_feature(beds: str, names: str, modes: str, overlaps: str) -> list:
    """Parse input features and settings"""
    logger = logging.getLogger(__name__)

    bed_list = beds.split(',')
    name_list = names.split(',')
    mode_list = modes.split(',')
    overlap_list = [int(x) for x in overlaps.split(',')]


    # check input
    if len(bed_list) != len(name_list) != len(mode_list) != len(overlap_list):
        logger.error("Number of features, names, modes and overlaps are not the same")
        raise SystemExit()
    for x in overlap_list:
        if x < 1 and x != -1:
            logger.error("Overlap value should be >= 1 or -1")
            raise SystemExit()
    
    features = []
    for i in range(len(bed_list)):
        features.append(Feature_range(bed_list[i], name_list[i], mode_list[i], overlap_list[i]))

    return features


def add_feature_header(vcf: pysam.VariantFile, feature_list: list) -> pysam.VariantHeader:
    """Add header of features"""
    header = vcf.header
    for feature in feature_list:
        if feature.mode == 'a':
            header.info.add(feature.name, "1", "Integer", f"Amount of overlapped bp between SV and {feature.name}")
        elif feature.mode == 'b':
            header.info.add(feature.name, "1", "Integer", f"Whether SV overlap with {feature.name} (1: yes, 0: no)")
    return header


def intersect_features(invcf: pysam.VariantFile, feature_list: list, df_invcf: pd.DataFrame=None) -> pd.DataFrame:
    """Intersect SV with features and return dataframe"""
    if df_invcf is None:
        df_out = vcf_to_df(invcf, info=['SVTYPE', 'SVLEN'])
        df_out = df_out.set_index('ID')
    else:
        df_out = df_invcf.copy()
    pr_invcf = vcf_to_pyranges(invcf)

    for feature in feature_list:
        df_intersect = feature.intersect(pr_invcf)
        if feature.mode == 'b':
            df_intersect.loc[df_intersect[f'Overlap_{feature.name}'] > 0, f'Overlap_{feature.name}'] = 1
        df_intersect = df_intersect[['ID', f'Overlap_{feature.name}']].set_index('ID')
        df_out = df_out.join(df_intersect, how='left', on='ID')
        df_out[f'Overlap_{feature.name}'] = df_out[f'Overlap_{feature.name}'].fillna(0).astype(int)
    
    return df_out


def intersect_vcf(invcf: pysam.VariantFile, out: str, feature_list: list) -> None:
    """Intersect SV with features and add to vcf INFO"""
    df_invcf = intersect_features(invcf, feature_list)
    
    # write to outvcf
    header = add_feature_header(invcf, feature_list)
    outvcf = pysam.VariantFile(out, 'w', header=header)

    for variant in invcf.fetch():
        new_variant = variant.copy()
        for feature in feature_list:
            new_variant.info[feature.name] = df_invcf.at[variant.id, f'Overlap_{feature.name}'].item()

        outvcf.write(new_variant)
    
    outvcf.close()
    
    if out.endswith("vcf.gz"):
        pysam.bcftools.index("-t", out)


def intersectSV_main(cmdargs) -> None:
    """Main function"""
    logger = logging.getLogger(__name__)

    args = parse_cmdargs(parser, cmdargs)

    # parse input
    feature_list = parse_feature(args.bed, args.name, args.mode, args.overlap)

    # intersect
    invcf = read_vcf(args.invcf)
    if args.out_type == 'vcf':
        intersect_vcf(invcf, args.output, feature_list)
    elif args.out_type == 'bed':
        df_intersect = intersect_features(invcf, feature_list)
        df_intersect.to_csv(args.output, sep='\t')
    else:
        logger.error("Invalid output type")
        raise SystemExit()

    logger.info(f"Write output to {args.output}")

