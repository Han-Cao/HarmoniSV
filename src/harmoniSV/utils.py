#!/usr/bin/env python3

# Utils of harmoniSV
# Created: 19/9/2022
# Author: Han Cao

import logging
import sys
import argparse
import os

import pysam
import pysam.bcftools
import pandas as pd
import pyranges as pr


SVTYPE_ALLOW = {'INS', 'DEL', 'DUP', 'INV', "CNV", "BND"}

def read_vcf(file_vcf: str, 
             check_id: bool=False, 
             check_svtype: bool=False, 
             check_sample: bool=False, 
             check_n: int=10000,
             slient: bool=False) -> pysam.VariantFile:
    """Read vcf file"""
    logger = logging.getLogger(__name__)

    if not slient:
        logger.info(f"Read vcf file: {file_vcf}")
    vcf = pysam.VariantFile(file_vcf, 'r')

    if check_sample:
        vcf_samples = list(vcf.header.samples)
        if len(vcf_samples) > 1:
            logger.error(f"More than one sample in {file_vcf}")
            raise SystemExit()
    
    if (check_id or check_svtype) and check_n > 0:
        logger.info(f"Check format of first {check_n} variants")
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


class OutVcf:
    """ Class to write vcf file """
    def __init__(self, file_vcf: str, header: pysam.VariantHeader) -> None:
        self.file_vcf = file_vcf
        self.header = header
        self.logger = logging.getLogger(__name__)

        # set output mode
        if (file_vcf.endswith('vcf.gz') or file_vcf.endswith('vcf.bgz')):
            self.flag = 'bgz'
            self.file_out = f'{file_vcf}.tmp'
        elif file_vcf.endswith('bcf'):
            self.flag = 'bcf'
            self.file_out = f'{file_vcf}.tmp'
        else:
            self.flag = 'vcf'
            self.file_out = file_vcf
        
        # open file
        self.f = open(self.file_out, 'w')

        # write header
        self.f.write(str(self.header))
    

    def write(self, record: pysam.VariantRecord):
        """ Write a variant record """
        self.f.write(str(record))


    def close(self, slient: bool=False):
        """ Close file """
        self.f.close()

        if self.flag == 'bgz':
            pysam.bcftools.view('--no-version', '-Oz', '-o', self.file_vcf, self.file_out, catch_stdout=False)
            os.remove(self.file_out)
        elif self.flag == 'bcf':
            pysam.bcftools.view('--no-version', '-Ob', '-o', self.file_vcf, self.file_out, catch_stdout=False)
            os.remove(self.file_out)

        if not slient:
            self.logger.info(f"Write output to {self.file_vcf}")        


def read_manifest(file_manifest: str, header, default_info: list=[]) -> pd.DataFrame:
    """
    Read manifest file (tab-delimited without header)

    Format:
    Column 1: <file_vcf>, one file per line
    Column 2 (optional): <info_tag>, comma separated list of additional INFO tags to be extracted, all if not specified
    """

    logger = logging.getLogger(__name__)

    # header = None for file list
    # header = 0 for manifest file
    df_manifest = pd.read_csv(file_manifest, sep='\t', header=header, dtype=str)

    # add header if not exist
    if header is None:
        if df_manifest.shape[1] == 1:
            df_manifest.columns = ['file']
        else:
            logger.error(f"Invalid file list, please make sure there is only one column without header")
            raise SystemExit()
        
    # add info column if not exist
    if 'info' not in df_manifest.columns:
        if len(default_info) > 0:
            df_manifest['info'] = ','.join(default_info)
            df_manifest['info'] = df_manifest['info'].str.split(',')
        else:
            df_manifest['info'] = 'all'
    # add default info to info
    elif len(default_info) > 0:
        df_manifest['info'].fillna('', inplace=True)
        df_manifest['info'] = df_manifest['info'].str.split(',')
        df_manifest['info'] = df_manifest['info'].apply(lambda x: list(set(default_info + x) - set([''])))
    
    return df_manifest
        

def manifest_to_df(df_manifest: pd.DataFrame, region: str=None) -> pd.DataFrame:
    """ Read VCFs in the manifest file and concat to dataframe """

    lst_df = []
    for i, row in df_manifest.iterrows():
        file_vcf = row['file']
        info = row['info']
        vcf = read_vcf(file_vcf)
        df_vcf = vcf_to_df(vcf, info, region)
        df_vcf['file_idx'] = i
        lst_df.append(df_vcf)
        vcf.close()

    df_all = pd.concat(lst_df, ignore_index=True)
    return df_all


class IdGenerator:
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
                tag_value = tag_value
        
        new_info[tag] = tag_value
    
    return new_info


def vcf_to_df(vcf: pysam.VariantFile, info: list='all', region: str=None, fill_tuple_na: bool=False) -> pd.DataFrame:
    """
    Convert vcf to dataframe, ID, CHR, POS and all INFO tags are included by default
    
    info: list of INFO tags to be extracted. If 'all', all INFO tags will be extracted (default)
    region: bcftools region string (default: None)
    fill_tuple_na: whether to replace NaN if INFO tag is tuple (default: False)

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

    # handle NA for tuple types
    if fill_tuple_na:
        for col in df.columns:
            if df[col].dtype != 'object':
                continue
            n_na = df[col].isna().sum()
            if n_na == 0:
                continue
            # if tuple and include NaN, expand NaN to tuple of None
            if df[col].apply(lambda x: isinstance(x, tuple)).any():
                # count number of elements in tuple
                n_elem = df[col].apply(lambda x: len(x) if isinstance(x, tuple) else -1)
                n_elem_uniq = n_elem.unique()
                n_elem_uniq = n_elem_uniq[n_elem_uniq >= 0]
                # return (None,) if not fixed length
                if len(n_elem_uniq) > 1:
                    df.loc[df[col].isna(), col] = [tuple([None])] * n_na
                else:
                    df.loc[df[col].isna(), col] = df.loc[df[col].isna(), col].apply(lambda x: tuple([None] * n_elem_uniq[0]))

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

    if len(cmdargs) == 0:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args(cmdargs)
    # print arguments
    print(f'Command: {parser.prog} {" ".join(cmdargs)}', file=sys.stderr)
    return args


class ProgressLogger:
    """Log progress"""
    def __init__(self, item: str, verbosity: int=1000) -> None:
        self.counter = 0
        self.verbosity = verbosity
        self.item = item
        self.logger = logging.getLogger(__name__)
    
    def log(self):
        self.counter += 1
        if self.counter % self.verbosity == 0:
            self.logger.info(f"{self.counter} {self.item} processed")

    def finish(self):
        self.logger.info(f"{self.counter} {self.item} processed")