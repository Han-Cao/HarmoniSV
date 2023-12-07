#!/usr/bin/env python3

# Genotype SVs across SV genotyping methods
# Created: 25/8/2022
# Author: Han Cao

import logging
import argparse
from collections import defaultdict

import pysam
import pysam.bcftools
import pandas as pd
import numpy as np

from utils import read_vcf, OutVcf, read_manifest, vcf_to_df, ProgressLogger, parse_cmdargs

# parse arguments
parser = argparse.ArgumentParser(prog="harmonisv genotype",
                                 description="Genotype SVs across SV genotyping methods",
                                 add_help=False)
io_arg = parser.add_argument_group('Input/Output arguments')
io_arg.add_argument("-i", "--invcf", metavar="VCF", type=str, required=True,
                    help="input representative SV call set vcf, INFO/ID_LIST store merged SV calls")
io_arg.add_argument("-f", "--manifest", metavar="TSV", type=str, required=True,
                    help="tab separated manifest file of SV calling results. Column headers should be: file, sample, aligner, caller, info, is_force_call")
io_arg.add_argument("--sample", metavar="SAMPLE_ID", type=str, required=True,
                    help="sample ID to be extracted from manifest")
io_arg.add_argument("-o", "--outvcf", metavar="VCF", type=str, required=True,
                    help="output vcf")
io_arg.add_argument("-r", "--region", metavar="chr", type=str, default=None,
                    help="genomic region to genotype (require indexed VCF)")

genotyping_arg = parser.add_argument_group('Genotyping arguments')
genotyping_arg.add_argument("--min-dp", metavar="10", type=int, default=10,
                            help="minimum depth to to genotype SVs as 0/0 if no method has RE>0, otherwise ./.")
genotyping_arg.add_argument("--homozygous-freq", metavar="0.8", type=float, default=0.8,
                            help="minimum average allele frequency to genotype as homozygous")
genotyping_arg.add_argument("--heterozygous-freq", metavar="0.2", type=float, default=0.2,
                            help="minimum average allele frequency to genotype as heterozygous")
genotyping_arg.add_argument("--include", metavar="all", type=str, default="all",
                            help="methods to be used to determine genotype. Options: 'all', 'force_call', or comma-separated list of methods ('ALIGNER_CALLER'). Default: 'all'")

optional_arg = parser.add_argument_group('optional arguments')
optional_arg.add_argument("-h", "--help", action='help', help='show this help message and exit')


class GenotypeManifest:
    """ Parse manifest table and store metadata """
    def __init__(self, df_manifest: pd.DataFrame, sample: str, region: str=None) -> None:
        self.df = df_manifest
        self.sample = sample
        self.methods = set()
        self.force_methods = set()
        self.discovery_methods = set()
        self.info_type = {}
        self._check_manifest()

        # only keep rows of working sample
        self.df = self.df[self.df["sample"] == sample].copy()

        # save manifest to dict
        self.vcf_dict = defaultdict(dict)
        self.method_dict = defaultdict(dict)
        for _, row in self.df.iterrows():
            self._parse_row(row, region)
    
    def _check_manifest(self):
        """ Check the format of manifest table """
        logger = logging.getLogger(__name__)
        # check header
        col_required = set(["file", "sample", "aligner", "caller", "info", "is_force_call"])
        col_exist = set(self.df.columns)
        col_missing = col_required - col_exist
        if len(col_missing) > 0:
            logger.error(f"Columns {col_missing} are not in the manifest file")
            raise SystemExit()

        # check duplicated method
        df_count = self.df.groupby(["sample", "aligner", "caller", "is_force_call"]).size().reset_index(name='count')
        df_dup = df_count[df_count["count"] > 1]
        if df_dup.shape[0] > 0:
            logger.error(f"Duplicate method in manifest file: \n {df_dup.to_string()}")
            raise SystemExit()

        # check duplicated file
        file_count = self.df['file'].value_counts()
        file_dup = file_count[file_count > 1]
        if file_dup.shape[0] > 0:
            logger.error(f"Duplicate file in manifest file: \n {file_dup.tolist()}")
            raise SystemExit()
    
    def _parse_vcf(self, vcf: pysam.VariantFile, info: str, region: str=None) -> pd.DataFrame:
        """ Read VCF and conver to df """

        df_vcf = vcf_to_df(vcf, info, region)
        df_vcf.set_index('ID', inplace=True)
        df_vcf['AF'] = df_vcf['RE'] / df_vcf['DP']
        # AF = NaN if DP = 0 is correct for single-method genotyping
        # However, if any other method call this SV (i.e., DP and AF > 0) 
        # setting AF = 0 means this method does not call this SV
        # We will always set genotype as ./. if no method call this SV or DP < min_depth
        # Therefore, setting AF = 0 will not call missing genotype as 0/0
        df_vcf['AF'] = df_vcf['AF'].fillna(0)

        return df_vcf
    
    def _parse_row(self, row: pd.Series, region: str=None) -> None:
        """ Parse one row of manifest table """
        aligner = row["aligner"]
        caller = row["caller"]
        file = row["file"]
        info = row["info"]
        vcf = read_vcf(file, check_sample=True)
        df_vcf = self._parse_vcf(vcf, info, region)
        method = f"{aligner}_{caller}".upper()
        # if AF is not in INFO, add it
        if "AF" not in info:
            info.append("AF")

        # append method
        if method not in self.methods:
            self.methods.add(method)
            self.method_dict[method]["aligner"] = aligner
            self.method_dict[method]["caller"] = caller

        # append file
        self.vcf_dict[file]["path"] = file
        self.vcf_dict[file]["method"] = method
        self.vcf_dict[file]["info"] = info
        self.vcf_dict[file]["df"] = df_vcf
        self.vcf_dict[file]["header"] = vcf.header

        # group by is_force_call
        if row["is_force_call"]:
            self.force_methods.add(method)
            self.method_dict[method]["has_force"] = True
            self.method_dict[method]["force"] = self.vcf_dict[file]
        else:
            self.discovery_methods.add(method)
            self.method_dict[method]["has_discovery"] = True
            self.method_dict[method]["discovery"] = self.vcf_dict[file]
            # add SV ID to discovery SV set
            self.vcf_dict[file]['sv_id'] = set(df_vcf.index)
        
        vcf.close()

    def get_aligner(self, method: str) -> str:
        """ Get aligner of method """
        return self.method_dict[method]["aligner"]
    
    def get_caller(self, method: str) -> str:
        """ Get caller of method """
        return self.method_dict[method]["caller"]
    
    def has_force(self, method: str) -> bool:
        """ Check if method has force calling """
        return self.method_dict[method]["has_force"]

    def has_discovery(self, method: str) -> bool:
        """ Check if method has discovery calling """
        return self.method_dict[method]["has_discovery"]
    
    def map_discovery(self, id_list: list) -> dict:
        """ Map discovery SVs to method """
        res = {}
        for method in self.discovery_methods:
            overlap_id = self.method_dict[method]["discovery"]["sv_id"].intersection(id_list)
            for id in overlap_id:
                res[id] = method
        return res
    
    def get_info(self, method: str, is_force: bool) -> str:
        """ Get info of method """
        if is_force:
            return self.method_dict[method]["force"]["info"]
        else:
            return self.method_dict[method]["discovery"]["info"]
    
    def query(self, method: str, id: str, is_force: bool) -> pd.Series:
        """ Query SV by ID """
        try:
            if is_force:
                res = self.method_dict[method]["force"]["df"].loc[id]
            else:
                res = self.method_dict[method]["discovery"]["df"].loc[id]
        except KeyError:
            res = None
        
        return res

    def add_info_type(self, tag: str, type: str) -> None:
        """ Add info type """
        self.info_type[tag] = type



def parse_gt_method(method: str):
    """ Parse genotyping method """
    if method == "all":
        return "all"
    elif method == "force_call":
        return "force_call"
    elif "," in method:
        return [x.upper() for x in method.split(",")]
    else:
        logger = logging.getLogger(__name__)
        logger.error(f"Unknown genotyping method: {method}")
        raise SystemExit()


def add_header(header: pysam.VariantHeader, manifest: GenotypeManifest) -> pysam.VariantHeader:
    """ Add VCF header for methods """
    # mandatory INFO tags
    header.info.add("SUPP_METHOD", "1", "Integer", "Number of methods support this SV")
    header.info.add("SUPP_METHOD_FORCE", "1", "Integer", "Number of force calling methods support this SV")
    header.info.add("SUPP_ALIGNER", "1", "Integer", "Number of aligners support this SV")
    header.info.add("SUPP_CALLER", "1", "Integer", "Maximum number of callers support this SV given the same aligner")
    header.info.add("SUPP_CALLER_FORCE", "1", "Integer", "Maximum number of force calling callers support this SV given the same aligner")
    header.info.add("MEAN_AF", "1", "Float", "Average allele frequency of all methods")
    header.info.add("STD_AF", "1", "Float", "Standard deviation of allele frequency of all methods")
    header.info.add("MEAN_AF_CALL", "1", "Float", "Average allele frequency of methods support this SV")
    header.info.add("STD_AF_CALL", "1", "Float", "Standard deviation of allele frequency of methods support this SV")
    header.info.add("MAX_RE", "1", "Integer", "Maximum number of reads support this SV")

    # additional INFO tags
    for _, vcf_dict in manifest.vcf_dict.items():
        method = vcf_dict["method"]
        info = vcf_dict["info"]
        old_header = vcf_dict["header"]

        for tag in info:
            new_tag = f"{tag}_{method}"
            if new_tag in header.info:
                continue
            # extract from from the original header
            old_header_line = old_header.info[tag]
            header.info.add(new_tag, old_header_line.number, old_header_line.type, old_header_line.description)
            # add info type to manifest
            manifest.add_info_type(new_tag, old_header_line.type)

    return header


# TODO: round float
def np2value(x, type: str, na_value=None):
    """ Convert the type required by VCF, return specific value when na """
    # np.isnan(str) will cause error
    if not isinstance(x, str) and np.isnan(x):
        return na_value

    if type == 'Integer':
        return int(x)
    elif type == 'Float':
        return float(x)
    elif type == 'String':
        return str(x)
    elif type == 'Flag':
        return bool(x)
    

def extract_sv_call(variant: pysam.VariantRecord, 
                    manifest: GenotypeManifest,
                    hetero_freq: float) -> tuple:
    """ Extract info from discovery and force calling results """
    # store INFO of each methods
    info_dict = defaultdict(dict)
    exist_methods = set()
    exist_methods_force = set()

    # store number of support methods and callers
    support_methods_set = set()
    support_methods_force_set = set()
    support_callers_set = defaultdict(set)
    support_callers_force_set = defaultdict(set)
    support_callers_set["no_caller"] = set()
    support_callers_force_set["no_caller"] = set()

    # extract force calling INFO if any
    for method in manifest.force_methods:
        sv_info = manifest.query(method, variant.id, is_force=True)
        # skip if missing
        if sv_info is None:
            continue
        # add INFO to dict
        exist_methods.add(method)
        exist_methods_force.add(method)
        for tag in manifest.get_info(method, is_force=True):
            info_dict[tag][method] = sv_info[tag]
        # add method to support dict
        if info_dict["AF"][method] >= hetero_freq:
            support_methods_set.add(method)
            support_methods_force_set.add(method)
            aligner = manifest.get_aligner(method)
            support_callers_set[aligner].add(method)
            support_callers_force_set[aligner].add(method)   
    
    # extract discovery method INFO for method without force calling
    # Note: force calling might have missing site
    remain_methods = manifest.methods - exist_methods

    # extract discovery SV calls of working sample
    # TODO: if additional INFO tags in discovery, always add them
    if len(remain_methods) > 0:
        discovery_call = variant.info["ID_LIST"]
        if isinstance(discovery_call, str):
            discovery_call = [discovery_call]
        discovery_call_map = manifest.map_discovery(discovery_call)

        if len(discovery_call_map) > 0:
            for id, method in discovery_call_map.items():
                if method in exist_methods:
                    continue
                sv_info = manifest.query(method, id, is_force=False)
                exist_methods.add(method)
                # add INFO to dict
                for tag in manifest.get_info(method, is_force=False):
                    info_dict[tag][method] = sv_info[tag]
                # add method to support dict
                if info_dict["AF"][method] >= hetero_freq:
                    support_methods_set.add(method)
                    aligner = manifest.get_aligner(method)
                    support_callers_set[aligner].add(method)
        
    # set still missing methods with DP=., RE=0, AF=0
    remain_methods = manifest.methods - exist_methods
    if len(remain_methods) > 0:
        for method in remain_methods:
            info_dict["DP"][method] = np.nan
            info_dict["RE"][method] = 0
            info_dict["AF"][method] = 0
    
    # support dict
    support_dict = {"methods": support_methods_set,
                    "methods_force": support_methods_force_set,
                    "callers": support_callers_set,
                    "callers_force": support_callers_force_set,
                    "exist_force": exist_methods_force} # exist force is used when genotyping

    return info_dict, support_dict


def add_info(new_var: pysam.VariantRecord,
             old_var: pysam.VariantRecord,
             info_dict: dict,
             manifest: GenotypeManifest,
             keep_info: list) -> pysam.VariantRecord:
    """ Add INFO to variant """

    # keep original INFO if any
    for tag in keep_info:
        if tag in old_var.info:
            new_var.info[tag] = old_var.info[tag]
    
    # add new INFO tags
    for tag, method_values in info_dict.items():
        for method, value in method_values.items():
            new_tag = f"{tag}_{method}"
            new_var.info[new_tag] = np2value(value, manifest.info_type[new_tag])
    
    return new_var


def genotype_variant(variant: pysam.VariantRecord, 
                     manifest: GenotypeManifest,
                     min_depth: int, 
                     homo_freq: float, 
                     hetero_freq: float,
                     gt_method: str,
                     sample_col: str="Sample") -> pysam.VariantRecord:
    """ Genotype variant from multiple calling pipelines """

    # new VariantRecord object for output
    new_variant = variant.copy()

    # clean old info
    new_variant.qual = None
    new_variant.filter.clear()
    new_variant.info.clear()
    new_variant.samples[sample_col].clear() # TODO: clear format

    # store INFO of each methods
    info_dict, support_dict = extract_sv_call(variant, manifest, hetero_freq)

    # add INFO
    new_variant = add_info(new_variant, variant, info_dict, manifest, 
                           keep_info=["SVTYPE", "SVLEN", "END"])
    
    # add statistics
    new_variant.info["SUPP_METHOD"] = len(support_dict["methods"])
    new_variant.info["SUPP_METHOD_FORCE"] = len(support_dict["methods_force"])
    new_variant.info["SUPP_CALLER"] = max([len(x) for x in support_dict["callers"].values()])
    new_variant.info["SUPP_CALLER_FORCE"] = max([len(x) for x in support_dict["callers_force"].values()])

    # maximum RE
    new_variant.info["MAX_RE"] = int(np.nanmax(list(info_dict["RE"].values())))

    # calculate AF and estimate genotype
    if gt_method == "all":
        geno_dp = list(info_dict["DP"].values())
        geno_af = list(info_dict["AF"].values())
    elif gt_method == "force_call":
        geno_dp = [info_dict["DP"][x] for x in support_dict["exist_force"]]
        geno_af = [info_dict["AF"][x] for x in support_dict["exist_force"]]
    elif isinstance(gt_method, list):
        geno_dp = [info_dict["DP"][x] for x in gt_method]
        geno_af = [info_dict["AF"][x] for x in gt_method]
    # TODO: determine to use x > 0 or x >= hetero_freq
    geno_af_call = [x for x in geno_af if x >= hetero_freq]

    new_variant.info["MEAN_AF"] = np.nanmean(geno_af)
    new_variant.info["STD_AF"] = np.nanstd(geno_af)
    if len(geno_af_call) > 0:
        mean_af_call = np.nanmean(geno_af_call)
        new_variant.info["MEAN_AF_CALL"] = mean_af_call
        new_variant.info["STD_AF_CALL"] = np.nanstd(geno_af_call)
        # genotype by average AF if any method call this SV
        if mean_af_call < hetero_freq:
            new_variant.samples[sample_col]["GT"] = (0, 0)
            new_variant.info["AC"] = 0
        elif mean_af_call >= hetero_freq and mean_af_call < homo_freq:
            new_variant.samples[sample_col]["GT"] = (0, 1)
            new_variant.info["AC"] = 1
        else:
            new_variant.samples[sample_col]["GT"] = (1, 1)
            new_variant.info["AC"] = 2
    else:
        # genotype 0/0 or ./. by DP if no method call this SV
        # TODO: determine if use average DP
        new_variant.info["MEAN_AF_CALL"] = 0
        new_variant.info["STD_AF_CALL"] = None
        if np.nanmean(geno_dp) >= min_depth:
            new_variant.samples[sample_col]["GT"] = (0, 0)
        else:
            new_variant.samples[sample_col]["GT"] = (None, None)
        new_variant.info["AC"] = 0

    return new_variant



def genotype_vcf(invcf: pysam.VariantFile, 
                 outvcf: pysam.VariantFile,
                 manifest: GenotypeManifest,
                 region: str,
                 min_depth: int, 
                 homo_freq: float, 
                 hetero_freq: float,
                 gt_method: str,
                 sample_col: str="Sample", 
                 verbosity: int=1000) -> None:
    """ Genotyping all variants in the input callset for one sample """

    progess = ProgressLogger(item="variants", verbosity=verbosity)
    
    if region is None:
        vcf_iter = invcf.fetch()
    else:
        vcf_iter = invcf.fetch(region=region)

    for variant in vcf_iter:
        new_variant = genotype_variant(variant=variant, 
                                       manifest=manifest,
                                       min_depth=min_depth, 
                                       homo_freq=homo_freq, 
                                       hetero_freq=hetero_freq, 
                                       gt_method=gt_method,
                                       sample_col=sample_col)
        progess.log()
        outvcf.write(new_variant)
    
    progess.finish()


def genotypeSV_main(cmdargs) -> None:
    """ Genotype sample from input files """
    logger = logging.getLogger(__name__)

    args = parse_cmdargs(parser, cmdargs)

    # read invcf
    invcf = read_vcf(args.invcf, check_sample=True)
    sample_col = list(invcf.header.samples)[0]

    # read manifest
    df_manifest = read_manifest(args.manifest, header=0, default_info=['DP', 'RE'])
    df_manifest['is_force_call'] = df_manifest['is_force_call'].astype(int).astype(bool)
    manifest = GenotypeManifest(df_manifest, args.sample, args.region)

    # parse genotyping method
    gt_method = parse_gt_method(args.include)

    # prepare output vcf
    new_header = add_header(invcf.header, manifest)
    outvcf = OutVcf(args.outvcf, new_header)

    # genotype sample
    genotype_vcf(invcf=invcf, 
                 outvcf=outvcf, 
                 manifest=manifest,
                 region = args.region,
                 min_depth=args.min_dp, 
                 homo_freq=args.homozygous_freq, 
                 hetero_freq=args.heterozygous_freq,
                 sample_col=sample_col,
                 gt_method=gt_method)

    invcf.close()
    outvcf.close()
    
    logger.info('Done')