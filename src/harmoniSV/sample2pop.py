#!/usr/bin/env python3

# Convert single-sample VCF to multi-sample VCF
# Created: 23/12/2022
# Author: Han Cao

import argparse
import logging
from fnmatch import fnmatchcase
import os

import pysam
import pandas as pd
import numpy as np

from utils import read_vcf, parse_info, parse_cmdargs

# parse arguments
parser = argparse.ArgumentParser(prog="harmonisv sample2pop",
                                 description="Convert single-sample VCF to multi-sample VCF",
                                 add_help=False)
required = parser.add_argument_group('required arguments')
required.add_argument("-i", "--invcf", metavar="vcf", type=str, required=True,
                      help="input vcf file")
required.add_argument("-o", "--outvcf", metavar="PREFIX", type=str, required=True,
                      help="ouput vcf file")


rule_args = parser.add_argument_group('Rules to convert or merge variants, accept wildcards "*"')
rule_args.add_argument("--filter-GT", action="store_true", required=False, default=False,
                       help="set GT to ./. if FILTER is not missing or PASS")
rule_args.add_argument("--info-first", metavar="INFO", type=str, required=False,
                       help="comma separated INFO field where the first value is used")
rule_args.add_argument("--info-sum", metavar="INFO", type=str, required=False,
                       help="comma separated INFO field where the sum of values is used")
rule_args.add_argument("--info-avg", metavar="INFO", type=str, required=False,
                       help="comma separated INFO field where the average of values is used")
rule_args.add_argument("--info-min", metavar="INFO", type=str, required=False,
                       help="comma separated INFO field where the minimum of values is used")
rule_args.add_argument("--info-max", metavar="INFO", type=str, required=False,
                       help="comma separated INFO field where the maximum of values is used")
rule_args.add_argument("--info-concat", metavar="INFO", type=str, required=False, 
                       help="comma separated INFO field where the non-missing values are concatenated by comma")
rule_args.add_argument("--info-to-format", metavar="INFO", type=str, required=False,
                       help="comma separated INFO field where the values are moved to FORMAT field, FILTER will always set as FORMAT/FT_SAMPLE")
rule_args.add_argument("--keep-format", metavar="FORMAT", type=str, required=False, default="GT",
                       help="comma separated FORMAT field to keep. (Default: GT)")

optional_arg = parser.add_argument_group('optional arguments')
optional_arg.add_argument("-r", "--region", metavar="chr[:start-end]", type=str, required=False,
                          help="region string accecpted by bcftools")
optional_arg.add_argument("-d", "--delimiter", metavar=".", type=str, required=False, default=".",
                          help="delimiter for between sample ID and unique variant ID (default: .)")
optional_arg.add_argument("--seed", metavar="N", type=int, required=False, default=42,
                          help="random seed (default: 42)")
optional_arg.add_argument("-h", "--help", action='help', help='show this help message and exit')



def rule_first(x: pd.Series):
    """Return first value of pd.Series"""
    return x.iloc[0]


def rule_sum(x: pd.Series):
    """Return sum of pd.Series"""
    if x.isna().sum() == len(x):
        return np.nan
    else:
        return np.nansum(x)


def rule_avg(x: pd.Series):
    """Return average of pd.Series"""
    if x.isna().sum() == len(x):
        return np.nan
    else:
        return np.nanmean(x)


def rule_min(x: pd.Series):
    """Return minimum of pd.Series"""
    if x.isna().sum() == len(x):
        return np.nan
    else:
        return np.nanmin(x)


def rule_max(x: pd.Series):
    """Return maximum of pd.Series"""
    if x.isna().sum() == len(x):
        return np.nan
    else:
        return np.nanmax(x)

def rule_concat(x: pd.Series):
    """Return concatenated string of pd.Series"""
    x_clean = x.dropna()
    if len(x_clean) == 0:
        return '.'
    else:
        return ','.join(x_clean.astype(str))

RULE_DICT = {'first': rule_first,
             'sum': rule_sum,
             'avg': rule_avg,
             'min': rule_min,
             'max': rule_max,
             'concat': rule_concat}

GT_DICT = {0: '0',
           1: '1',
           None: '.'}

VCF_BASIC_FIELDS = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']



def match_field(field_str: str, all_fields: list) -> list:
    """Match field string"""
    logger = logging.getLogger(__name__)
    if field_str is None:
        return []
    field_list = field_str.split(',')

    result_list = []
    for field in field_list:
        if '*' not in field:
            if field not in all_fields:
                logger.warning(f"Field {field} not in vcf, ignored")
            else:
                result_list.append(field)
        else:
            matched_fields = [x for x in all_fields if fnmatchcase(x, field)]
            if len(matched_fields) == 0:
                logger.warning(f"No field matched {field}, ignored")
                continue
            else:
                result_list.extend(matched_fields)

    return result_list

def parse_rule(rule: str, field_str: str, all_fields: list, rule_dict: dict=RULE_DICT) -> dict:
    """Parse rule string to dict"""
    matched_fields = match_field(field_str, all_fields)

    result_dict = {}
    for matched_field in matched_fields:
        result_dict[matched_field] = rule_dict[rule]

    return result_dict


def parse_all_rules(args, all_fields:list) -> dict:
    """Parse all rules"""
    # basic variant information
    merge_rule = {'CHROM': RULE_DICT['first'],
                  'POS': RULE_DICT['first'],
                  'ID': RULE_DICT['first'],
                  'REF': RULE_DICT['first'],
                  'ALT': RULE_DICT['first'],
                  'QUAL': RULE_DICT['first'],
                  'FILTER': RULE_DICT['first'],
                  'AC': RULE_DICT['sum'],
                  'AN': RULE_DICT['sum']}
    
    # rule first
    if args.info_first is not None:
        merge_rule.update(parse_rule('first', args.info_first, all_fields=all_fields))
    if args.info_sum is not None:
        merge_rule.update(parse_rule('sum', args.info_sum, all_fields=all_fields))
    if args.info_avg is not None:
        merge_rule.update(parse_rule('avg', args.info_avg, all_fields=all_fields))
    if args.info_min is not None:
        merge_rule.update(parse_rule('min', args.info_min, all_fields=all_fields))
    if args.info_max is not None:
        merge_rule.update(parse_rule('max', args.info_max, all_fields=all_fields))
    if args.info_concat is not None:
        merge_rule.update(parse_rule('concat', args.info_concat, all_fields=all_fields))
    
    return merge_rule


def sample_vcf_to_df(vcf: pysam.VariantFile, 
                    info:list, 
                    format: list, 
                    info_to_format: list, 
                    filter: bool=True, 
                    region: str=None):
    """Convert single-sample vcf to dataframe"""
    rows_info = []
    rows_geno = []
    
    # pre-process info and format
    sample_id = vcf.header.samples[0]
    if format != 'all':
        format = [x for x in format if x != 'GT']
    # remove AC and AN from info
    info = [x for x in info if x not in ['AC', 'AN']]
    # if info_to_format is not None:
    #     info = [x for x in info if x not in info_to_format]

    # generate vcf iterator
    if region is None:
        vcf_iter = vcf.fetch()
    else:
        vcf_iter = vcf.fetch(region=region)

    # parse variants
    for variant in vcf_iter:
        # parse INFO
        new_row_info = {}
        new_row_info['CHROM'] = variant.chrom
        new_row_info['POS'] = variant.pos
        new_row_info['uniq_id'] = variant.id
        new_row_info['REF'] = variant.ref
        new_row_info['ALT'] = ','.join(variant.alts)
        new_row_info['QUAL'] = '.' # do not use QUAL for now
        new_row_info['FILTER'] = '.' # add filter after analyze all samples
        new_row_info['END'] = variant.stop


        # parse genotype
        new_row_geno = {}
        new_row_geno['uniq_id'] = variant.id
        
        if format == 'all':
            tag_exist = variant.format.keys()
            tag_exist.remove('GT')
        else:
            tag_exist = [x for x in format if x in variant.format]

        # handle GT
        var_filter = variant.filter.keys()
        filter_flag = False
        if len(var_filter) == 0:
            var_filter = '.'
        else:
            var_filter = ','.join(var_filter)
            if var_filter != 'PASS' and filter:
                filter_flag = True
            
        gt_value = variant.samples[sample_id]['GT']
        gt_an = 2 - gt_value.count(None)
        gt_ac = gt_an - gt_value.count(0)
        if filter_flag:
            gt_value = './.'
            gt_ac = 0
            gt_an = 0
        elif variant.samples[sample_id].phased:
            gt_value = '|'.join([GT_DICT[x] for x in gt_value])
        else:
            gt_value = '/'.join([GT_DICT[x] for x in gt_value])
        new_row_geno['GT'] = gt_value
        new_row_geno['FT_SAMPLE'] = var_filter

        # add AC and AN
        new_row_info['AC'] = gt_ac
        new_row_info['AN'] = gt_an

        # parse remaining info
        new_row_info.update(parse_info(variant, info))
        rows_info.append(new_row_info)
        

        # parse format tags
        for tag in tag_exist:
            tag_value = variant.samples[sample_id][tag]

            if isinstance(tag_value, tuple):
                if len(tag_value) == 1:
                    tag_value = tag_value[0]
                else:
                    tag_value = [str(x) for x in tag_value]
                    tag_value = ','.join(tag_value)
            
            new_row_geno[tag] = tag_value
        
        # move info tags to format
        for tag in info_to_format:
            if tag in variant.info:
                tag_value = variant.info[tag]
                if isinstance(tag_value, tuple):
                    if len(tag_value) == 1:
                        tag_value = tag_value[0]
                    else:
                        tag_value = [str(x) for x in tag_value]
                        tag_value = ','.join(tag_value)
                        
                new_row_geno[tag] = tag_value

             
        rows_geno.append(new_row_geno)
    
    # convert to dataframe, fill None with np.nan
    df_info = pd.DataFrame(rows_info).fillna(value=np.nan)
    df_geno = pd.DataFrame(rows_geno).fillna(value=np.nan)
    
    return df_info, df_geno


def generate_info(df: pd.DataFrame) -> pd.DataFrame:
    """Generate INFO fields from dataframe"""
    df_new = df.copy()

    info_col = [x for x in df.columns if x not in VCF_BASIC_FIELDS]
    float_col = [x for x in info_col if df[x].dtype == 'float64']
    for info in info_col:
        if info in float_col:
            df_new[info] = df_new[info].round(6)
        df_new[info] = df_new[info].fillna('.').astype(str)
        df_new[info] = info + '=' + df_new[info]
    
    return df_new


def generate_fomart(df: pd.DataFrame) -> pd.DataFrame:
    """Generate FORMAT fields from dataframe"""
    df_new = df.copy()

    format_col = [x for x in df.columns if x not in VCF_BASIC_FIELDS]
    float_col = [x for x in format_col if df[x].dtype == 'float64']
    for fmt in format_col:
        if fmt in float_col:
            df_new[fmt] = df_new[fmt].round(6)
        df_new[fmt] = df_new[fmt].fillna('.').astype(str)
    
    return df_new


def merge_vcf(df_info: pd.DataFrame, 
              df_geno: pd.DataFrame, 
              merge_rule: dict,
              keep_info: list,
              keep_format: list,
              sample_list: list) -> pd.DataFrame:
    """Merge vcf dataframe"""
    df_info = df_info.copy()
    df_geno = df_geno.copy()

    # merge info
    df_info = df_info.groupby('ID', sort=False).agg(merge_rule)
    df_info = generate_info(df_info)
    df_info['INFO'] = df_info[keep_info].apply(lambda x: ';'.join(x), axis=1)
    df_info['FORMAT'] = ':'.join(keep_format)
    df_info = df_info.reset_index()[VCF_BASIC_FIELDS + ['INFO', 'FORMAT']]

    # merge geno
    df_geno = df_geno.drop(columns=['uniq_id'])
    df_geno = generate_fomart(df_geno)
    df_geno['value'] = df_geno[keep_format].apply(lambda x: ':'.join(x), axis=1)

    df_geno = df_geno.pivot(index='ID', columns='Sample', values='value')
    df_geno = df_geno.reindex(sample_list, axis=1)

    # combine vcf
    df_vcf = df_info.merge(df_geno, left_on='ID', right_index=True, how='left')
    df_vcf = df_vcf.rename(columns={'CHROM': '#CHROM'})

    return df_vcf



def generate_merge_header(header: pysam.VariantHeader, merge_rule: dict, info_to_format: list) -> list:
    """Create new header for merged vcf"""

    new_header = []

    # add original header record, change INFO type after merge
    for record in header.records:
        if record.key == 'INFO':
            name = record['ID']
            number = record['Number']
            description = record['Description']
            merge_f = merge_rule.get(name, None)
            if merge_f is RULE_DICT['avg']:
                new_record = f'##INFO=<ID={name},Number={number},Type=Float,Description={description}>\n'
            elif merge_f is RULE_DICT['concat']:
                new_record = f'##INFO=<ID={name},Number={number},Type=String,Description={description}>\n'
            else:
                new_record = str(record)
            new_header.append(new_record)
        else:
            new_header.append(str(record))
    
    # add AC and AN if not exist
    if 'AC' not in header.info:
        new_header.append('##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">\n')
    if 'AN' not in header.info:
        new_header.append('##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">\n')
    
    # add FORMAT tags
    new_header.append('##FORMAT=<ID=FT_SAMPLE,Number=1,Type=String,Description="Genotype filter per sample">\n')
    for name in info_to_format:
        number = header.info[name].number
        type_ = header.info[name].type
        description = header.info[name].description
        new_record = f'##FORMAT=<ID={name},Number={number},Type={type_},Description={description}>\n'
        new_header.append(new_record)
    
    return new_header

def sample2pop_main(cmdargs) -> None:
    """Main function of sample2pop"""
    
    logger = logging.getLogger(__name__)

    args = parse_cmdargs(parser, cmdargs)

    # read input
    invcf = read_vcf(args.invcf)

    # extract all info and format tags
    all_info = invcf.header.info.keys()
    all_format = invcf.header.formats.keys()

    # match rules
    merge_rule = parse_all_rules(args, all_info)

    # extract tags to extact
    keep_info = [x for x in merge_rule.keys() if x not in VCF_BASIC_FIELDS]
    keep_format = match_field(args.keep_format, all_format)
    keep_info2format = match_field(args.info_to_format, all_info)

    df_info, df_geno = sample_vcf_to_df(vcf=invcf, 
                                        info=keep_info,
                                        format=keep_format, 
                                        info_to_format=keep_info2format, 
                                        filter=args.filter_GT, 
                                        region=args.region)
    
    # extract still remain tags as some only exist in header
    keep_info = [x for x in keep_info if x in df_info.columns]
    keep_format = [x for x in keep_format if x in df_geno.columns]
    keep_info2format = [x for x in keep_info2format if x in df_geno.columns]

    logger.info(f'INFO to merge: {keep_info}')
    logger.info(f'FORMAT to extract: {keep_format}')
    logger.info(f'INFO extract to FORMAT: {keep_info2format}')

    merge_rule = {k: merge_rule[k] for k in merge_rule.keys() if k in df_info.columns}
    keep_format = keep_format + ['FT_SAMPLE'] + keep_info2format
    
    # extract sample and SV ID
    df_info['ID'] = df_info['uniq_id'].str.split('.', n=1, expand=True)[1]
    df_geno[['Sample', 'ID']] = df_geno['uniq_id'].str.split('.', n=1, expand=True)
    sample_list = df_geno['Sample'].unique().tolist()
    sample_list.sort()

    # merge vcf
    logger.info('Start merging vcf...')
    df_vcf = merge_vcf(df_info=df_info,
                       df_geno=df_geno,
                       merge_rule=merge_rule,
                       keep_info=keep_info,
                       keep_format=keep_format,
                       sample_list=sample_list)
    logger.info('Merging vcf done.')
    
    # generate new header
    new_header = generate_merge_header(header=invcf.header,
                                       merge_rule=merge_rule,
                                       info_to_format=keep_info2format)
    logger.info('Generate new header')

    # write output vcf
    outfile = args.outvcf
    bgz_flag = False
    if outfile.endswith('.gz'):
        outfile = outfile[:-3]
        bgz_flag = True
    
    logger.info(f'Writing output to {outfile}')
    with open(outfile, 'w') as f:
        for line in new_header:
            f.write(line)
        df_vcf.to_csv(f, sep='\t', index=False)
    
    if bgz_flag:
        logger.info(f'Bgzip output to {outfile}.gz')
        pysam.tabix_compress(outfile, outfile + '.gz', force=True)
        os.remove(outfile)

    return None





