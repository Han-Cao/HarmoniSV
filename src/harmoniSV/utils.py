#!/usr/bin/env python3

# Utils of harmoniSV
# Created: 19/9/2022
# Author: Han Cao

import logging
import sys
import argparse

import pysam
import pandas as pd
import pyranges as pr


SVTYPE_ALLOW = {'INS', 'DEL', 'DUP', 'INV'}

def read_vcf(file_vcf: str, check_id: bool=False, check_svtype: bool=False, check_sample: bool=False, check_n: int=10000):
    """Read vcf file"""
    logger = logging.getLogger(__name__)

    logger.info(f"Read vcf file: {file_vcf}")
    vcf = pysam.VariantFile(file_vcf, 'r')

    if check_sample:
        vcf_samples = list(vcf.header.samples)
        if len(vcf_samples) > 1:
            logger.error(f"More than one sample in {file_vcf}")
            raise SystemExit()
    
    if (check_id or check_svtype) and check_n > 0:
        id_set = set()
        for i, variant in enumerate(vcf.fetch(), start=1):
            if check_id:
                if variant.id in id_set:
                    logger.error(f"Duplicate ID: {variant.id} in {file_vcf}")
                    raise SystemExit()
                else:
                    id_set.add(variant.id)
            
            if check_svtype:
                if variant.info['SVTYPE'] not in SVTYPE_ALLOW:
                    logger.error(f"SVTYPE {variant.info['SVTYPE']} of {variant.id} is not allowed in {file_vcf}")
                    raise SystemExit()
            
            if i >= check_n:
                break
    
    vcf.reset()
    return vcf


class ID_generator:
    """ Generate new SV IDs """
    def __init__(self, prefix: str, chr: str=None) -> None:
        self.prefix = prefix
        self.chr = chr
        self.id = {}

    def get(self, sv_type: str):
        """ Generate new ID and update ID No. """
        new_id = self.id.get(sv_type, 1)
        self.id[sv_type] = new_id + 1
        if self.chr is not None:
            return f"{self.prefix}.{sv_type}.{self.chr}_{new_id}"
        else:
            return f"{self.prefix}.{sv_type}.{new_id}"


def parse_region(region: str) -> tuple:
    """Parse region string"""
    logger = logging.getLogger(__name__)
    if ':' in region:
        chrom, pos = region.split(':')
        if '-' in pos:
            start, end = pos.split('-')
            return chrom, int(start), int(end)
        else:
            logger.error(f"Invalid region: {region}, must be in format of chr:start-end")
            raise SystemExit()
    else:
        return region, None, None


def parse_info(variant: pysam.VariantRecord, info: list='all') -> dict:
    """Parse INFO tags of a variant"""
    new_info = {}
    if info == 'all':
        tag_exist = variant.info.keys()
    else:
        tag_exist = [x for x in info if x in variant.info]

    for tag in tag_exist:
        tag_value = variant.info[tag]
        if isinstance(tag_value, tuple):
            if len(tag_value) == 1:
                tag_value = tag_value[0]
            else:
                tag_value = ','.join(tag_value)
        
        new_info[tag] = tag_value
    
    return new_info


def vcf_to_df(vcf: pysam.VariantFile, info: list='all', region: str=None) -> pd.DataFrame:
    """
    Convert vcf to dataframe, ID, CHR, POS and all INFO tags are included by default
    
    info: list of INFO tags to be extracted. If 'all', all INFO tags will be extracted (default)
    region: bcftools region string (default: None)

    For INFO tags with multiple value (i.e., comma in value), store as string
    """
    rows = []
    if region is None:
        vcf_iter = vcf.fetch()
    else:
        vcf_iter = vcf.fetch(region=region)

    for variant in vcf_iter:
        new_row = {}
        new_row['ID'] = variant.id
        new_row['CHR'] = variant.chrom
        new_row['POS'] = variant.pos
        new_row.update(parse_info(variant, info))
        rows.append(new_row)
    
    df = pd.DataFrame(rows)
    
    return df


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


def parse_cmdargs(parser: argparse.ArgumentParser, cmdargs: list) -> argparse.Namespace:
    """Parse arguments"""
    logger = logging.getLogger(__name__)
    if len(cmdargs) == 0:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args(cmdargs)
    # print arguments
    logger.info(f'Command: {parser.prog} {" ".join(cmdargs)}' )
    return args