#!/usr/bin/env python3

# Genotype SVs across SV genotyping methods
# Created: 25/8/2022
# Author: Han Cao

import pysam
import pysam.bcftools
import argparse
import pandas as pd
import numpy as np
import logging

from utils import read_vcf, parse_cmdargs

# parse arguments
parser = argparse.ArgumentParser(prog="harmonisv genotype",
                                 description="Genotype SVs across SV genotyping methods",
                                 add_help=False)
parser.add_argument("-i", "--invcf", metavar="vcf", type=str, required=True,
                    help="input representative SV call set vcf, INFO/ID_LIST store merged SV calls")
parser.add_argument("-o", "--outvcf", metavar="txt", type=str, required=True,
                    help="output vcf")
parser.add_argument("--sample", metavar="SAMPLE_ID", type=str, required=True,
                    help="sample ID to extract sample record in --sv-info and --method-table")
parser.add_argument("--method-table", metavar="FILE", type=str, required=True,
                    help="tab-separated file of bam files for depth quantification. Column headers should be: sample, aligner, caller, bam, min_qual, min_read_len")
parser.add_argument("--sv-info", metavar="tsv", type=str, required=True,
                    help="tab-separated file for SV discovery calling INFO. Mandatory fields: ID, CHR, AC, DP, RE")
parser.add_argument("--force-call-info", metavar="tsv", type=str, required=False,
                    help="tab-separated file for SV force calling INFO. Mandatory fields: ID, CHR, AC, DP, RE")
parser.add_argument("--min-dp", metavar="10", type=int, default=10,
                    help="minimum depth to to genotype SVs as 0/0 if no method has RE>0, otherwise ./.")
parser.add_argument("--count-missing-dp", action="store_true",
                    help="count depth of missing SV from bam files")
parser.add_argument("--homozygous-freq", metavar="0.8", type=float, default=0.8,
                    help="minimum average allele frequency to genotype as homozygous")
parser.add_argument("--heterozygous-freq", metavar="0.2", type=float, default=0.2,
                    help="minimum average allele frequency to genotype as heterozygous")
parser.add_argument("--genotyping-method", metavar="all", type=str, default="all",
                    help="methods to be used to determine genotype. Options: 'all', 'force_call', or comma-separated list of methods ('ALIGNER_CALLER'). Default: 'all'")

optional_arg = parser.add_argument_group('optional arguments')
optional_arg.add_argument("-h", "--help", action='help', help='show this help message and exit')



def read_filter(min_qual: int, min_read_len: int=0):
    """ Genderate read filter funtion for pysam read_callback filter """

    if(min_read_len == 0):
        def fn(read):
            return (read.mapping_quality >= min_qual and 
                    not (read.flag & (0x4 | 0x100 | 0x200 | 0x400)))
        
        return fn

    elif(min_read_len > 0):
        def fn(read):
            return (read.mapping_quality >= min_qual and 
                    not (read.flag & (0x4 | 0x100 | 0x200 | 0x400)) and 
                    read.infer_query_length() >= min_read_len)
        return fn
    
    else:
        logger = logging.getLogger(__name__)
        logger.error("Minimum requred read length must >= 0")
        raise SystemExit()


class method_table:
    """ Store methods for one sample and corresponding bam file and read filter DP quantification """
    def __init__(self, file: pd.DataFrame, sample: str) -> None:
        logger = logging.getLogger(__name__)
        self.table = pd.read_csv(file, sep="\t")
        self.table = self.table[self.table["sample"] == sample]
        if len(self.table.index) == 0:
            logger.error(f"{sample} is not in the bam table")
            raise SystemExit()

        # upper case method used in VCF INFO
        self.table["method"] = self.table["aligner"] + "_" + self.table["caller"]
        self.table["method"] = self.table["method"].str.upper()

        # sample.aligner.caller prefix used in force calling table
        self.table["prefix"] = self.table["sample"] + "." + self.table["aligner"] + "." + self.table["caller"]
        
        self.bam = {}     # dict to store bams of each method
        self.filter = {}  # dict to store read_filter of each method
        self.prefix = {}  # dict to store SV ID prefix of each method to search force calling table
        self.aligner = {} # dict to store aligner of each method
        self.caller = {}  # dict to store caller of each method
        for _, row in self.table.iterrows():
            method = row["method"]
            self.bam[method] = pysam.AlignmentFile(row["bam"], "rb")
            self.filter[method] = read_filter(row["min_qual"], row["min_read_len"])
            self.prefix[method] = row["prefix"]
            self.aligner[method] = row["aligner"]
            self.caller[method] = row["caller"]
        
        self.methods = list(self.bam.keys())

    def get_bam(self, method) -> pysam.AlignmentFile:
        return self.bam[method]
    
    def get_filter(self, method):
        return self.filter[method]

    def get_prefix(self, method) -> str:
        return self.prefix[method]
    
    def get_aligner(self, method) -> str:
        return self.aligner[method]
    
    def get_caller(self, method) -> str:
        return self.caller[method]
    
    def close(self):
        for bam in self.bam.values():
            bam.close()


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


def parse_sv_info(file: str) -> pd.DataFrame:
    """ Extract sample ID and pipeline """

    df_info = pd.read_csv(file, sep="\t", index_col="ID", na_values=[".", "NA"],
                          dtype={'ID': str, 'CHR': str, 'SVTYPE': str})
    split_ID = df_info.index.str.split(".")
    df_info["Sample"] = [x[0] for x in split_ID]
    df_info["Method"] = [f"{x[1]}_{x[2]}".upper() for x in split_ID] # change method to upper case to match VCF header

    return df_info


def add_header(header: pysam.VariantHeader, methods: list) -> pysam.VariantHeader:
    """ Add VCF header for methods """
    for x in methods:
        header.info.add(f"DP_{x}", "1", "Integer", f"Depth of method {x}")
        header.info.add(f"RE_{x}", "1", "Integer", f"Number of support reads of method {x}")
        header.info.add(f"AF_{x}", "1", "Float", f"Allele frequency of method {x}")
    
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

    return header


def np2value(x, type: str, na_value=None):
    """ Convert numpy type to value, return specific value when na """
    if np.isnan(x):
        return na_value
    elif type == "int":
        return int(x)
    elif type == "float":
        return float(x)
    elif type == "str":
        return str(x)


def genotype_variant(variant: pysam.VariantRecord, 
                     sample: str, 
                     sample_method: method_table,
                     df_info: pd.DataFrame, 
                     df_force_call: pd.DataFrame, 
                     force_call_methods: list,
                     min_depth: int, 
                     count_dp: bool,
                     homo_freq: float, 
                     hetero_freq: float,
                     gt_method: str,
                     sample_col: str="Sample") -> pysam.VariantRecord:
    """ Genotype variant from multiple calling pipelines """

    # new VariantRecord object for output
    new_variant = variant.copy()
    new_variant.id = f"{sample}.{variant.id}"

    # clean old info
    new_variant.qual = None
    new_variant.filter.clear()
    new_variant.info.clear()
    new_variant.samples[sample_col].clear() # TODO: clear format

    # copy necessary info from original variant
    new_variant.info["SVTYPE"] = variant.info["SVTYPE"]
    new_variant.info["SVLEN"] = variant.info["SVLEN"]
    if "END" in variant.info:
        new_variant.info["END"] = variant.info["END"]

    # extract discovery SV calls of working sample
    all_sv_call = variant.info["ID_LIST"]
    if isinstance(all_sv_call, str):
        all_sv_call = [all_sv_call]
    sample_sv_call = [x for x in all_sv_call if sample in x]

    # store DP, RE, AF of methods
    method_dict = {}
    method_dict["DP"] = {}
    method_dict["RE"] = {}
    method_dict["AF"] = {}

    # extract force calling INFO if any
    if len(force_call_methods)>0:
        for method in force_call_methods:
            # extract force calling INFO if any
            force_call_id = sample_method.get_prefix(method) + "." + variant.id
            if force_call_id in df_force_call.index:
                method_dict["DP"][method] = df_force_call.at[force_call_id, 'DP']
                method_dict["RE"][method] = df_force_call.at[force_call_id, 'RE']
                if method_dict["DP"][method] > 0:
                    method_dict["AF"][method] = method_dict["RE"][method] / method_dict["DP"][method]
                else:
                    method_dict["AF"][method] = 0

    # extract discovery method INFO for method without force calling
    # Note: force calling might have missing site, ALWAYS use method_dict["DP"].keys() to get existing samples
    remain_methods = [x for x in sample_method.methods if x not in method_dict["DP"].keys()]
    if len(sample_sv_call) > 0 and len(remain_methods) > 0:
        for sv_call in sample_sv_call:
            method = df_info.at[sv_call, "Method"]
            if method in remain_methods:
                method_dict["DP"][method] = df_info.at[sv_call, 'DP']
                method_dict["RE"][method] = df_info.at[sv_call, 'RE']
                method_dict["AF"][method] = method_dict["RE"][method] / method_dict["DP"][method]

    # set still missing methods with RE=0, AF=0, get DP of missing call from bam file if specified
    # TODO: determine if only quantify depth for INS and DEL
    missing_methods = [x for x in remain_methods if x not in method_dict["DP"].keys()]
    if len(missing_methods) > 0:
        for method in missing_methods:
            # get 0-based END of SV
            # in pysam, variant.pos is 1-based position shown in vcf, variant.start is 0-based
            if count_dp:
                if variant.info["SVTYPE"] == "INS":
                    variant_stop = variant.pos
                else:
                    variant_stop = variant.start + abs(variant.info["SVLEN"])
                # TODO: if multiple method with the same filter, only count depth once
                method_dict["DP"][method] = sample_method.get_bam(method).count(contig=variant.chrom,
                                                                                start=variant.start,
                                                                                stop=variant_stop,
                                                                                read_callback=sample_method.get_filter(method))
            else:
                method_dict["DP"][method] = np.nan
            method_dict["RE"][method] = 0
            method_dict["AF"][method] = 0

    # add DP, RE, AF per method to INFO
    # convert from numpy dtype to python built-in type
    for method in sample_method.methods:
        new_variant.info[f"DP_{method}"] = np2value(method_dict["DP"][method], 'int')
        new_variant.info[f"RE_{method}"] = np2value(method_dict["RE"][method], 'int')
        new_variant.info[f"AF_{method}"] = np2value(method_dict["AF"][method], 'float')
    
    # count number of support methods and callers
    support_methods_set = set()
    support_methods_force_set = set()
    support_callers = {'no_caller': 0}
    support_callers_force = {'no_caller': 0}
    for method, freq in method_dict["AF"].items():
        if freq >= hetero_freq:
            support_methods_set = support_methods_set.union([method])
            aligner = sample_method.get_aligner(method)
            support_callers[aligner] = support_callers.get(aligner, 0) + 1
            if method in force_call_methods:
                support_methods_force_set = support_methods_force_set.union([method])
                support_callers_force[aligner] = support_callers_force.get(aligner, 0) + 1
    

    new_variant.info["SUPP_METHOD"] = len(support_methods_set)
    new_variant.info["SUPP_METHOD_FORCE"] = len(support_methods_force_set)
    new_variant.info["SUPP_CALLER"] = max(support_callers.values())
    new_variant.info["SUPP_CALLER_FORCE"] = max(support_callers_force.values())

    # extract maximum RE
    method_RE = list(method_dict["RE"].values())
    new_variant.info["MAX_RE"] = int(np.nanmax(method_RE))

    # calculate AF and estimate genotype
    if gt_method == "all":
        method_AF = list(method_dict["AF"].values())
    elif gt_method == "force_call":
        method_AF = [method_dict["AF"][x] for x in force_call_methods]
    elif isinstance(gt_method, list):
        method_AF = [method_dict["AF"][x] for x in gt_method]
    # TODO: determine to use x > 0 or x >= hetero_freq
    method_AF_call = [x for x in method_AF if x >= hetero_freq]

    new_variant.info["MEAN_AF"] = np.nanmean(method_AF)
    new_variant.info["STD_AF"] = np.nanstd(method_AF)
    if len(method_AF_call) > 0:
        mean_AF_call = np.nanmean(method_AF_call)
        new_variant.info["MEAN_AF_CALL"] = mean_AF_call
        new_variant.info["STD_AF_CALL"] = np.nanstd(method_AF_call)
        # genotype by average AF if any method call this SV
        if mean_AF_call < hetero_freq:
            new_variant.samples[sample_col]["GT"] = (0, 0)
            new_variant.info["AC"] = 0
        elif mean_AF_call >= hetero_freq and mean_AF_call < homo_freq:
            new_variant.samples[sample_col]["GT"] = (0, 1)
            new_variant.info["AC"] = 1
        else:
            new_variant.samples[sample_col]["GT"] = (1, 1)
            new_variant.info["AC"] = 2
    else:
        # genotype 0/0 or ./. by DP if no method call this SV
        new_variant.info["MEAN_AF_CALL"] = 0
        if np.nanmean(list(method_dict["DP"].values())) >= min_depth:
            new_variant.samples[sample_col]["GT"] = (0, 0)
        else:
            new_variant.samples[sample_col]["GT"] = (None, None)
        new_variant.info["AC"] = 0

    return new_variant
        


def genotype_vcf(invcf: pysam.VariantFile, 
                 outvcf: pysam.VariantFile, 
                 sample: str, 
                 sample_method: method_table,
                 df_info: pd.DataFrame, 
                 df_force_call: pd.DataFrame,
                 min_depth: int, 
                 count_dp: bool,
                 homo_freq: float, 
                 hetero_freq: float,
                 gt_method: str,
                 sample_col: str="Sample", 
                 verbosity: int=1000) -> None:
    """ Genotyping all variants in the input callset for one sample """

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    force_call_methods = df_force_call["Method"].unique()

    for i, variant in enumerate(invcf.fetch()):
        if i % verbosity == 0:
            logger.info(f"Finish {i} SVs. Genotyping SV at {variant.chrom}:{variant.pos} ({variant.id})")
        new_variant = genotype_variant(variant=variant, 
                                       sample=sample,
                                       sample_method=sample_method,
                                       df_info=df_info, 
                                       df_force_call=df_force_call, 
                                       force_call_methods=force_call_methods,
                                       min_depth=min_depth, 
                                       count_dp=count_dp,
                                       homo_freq=homo_freq, 
                                       hetero_freq=hetero_freq, 
                                       gt_method=gt_method,
                                       sample_col=sample_col)
        outvcf.write(new_variant)
    
    logger.info(f"Finish genotyping of all {i+1} SVs")


def genotypeSV_main(cmdargs) -> None:
    """ Genotype sample from input files """
    logger = logging.getLogger(__name__)

    args = parse_cmdargs(parser, cmdargs)

    # read invcf
    invcf = read_vcf(args.invcf, check_sample=True)
    sample_col = list(invcf.header.samples)[0]

    sample = args.sample

    # read SV info table
    logger.info(f"Read SV calling INFO of sample {sample} from {args.sv_info}")
    df_info = parse_sv_info(args.sv_info)
    df_info = df_info.loc[df_info["Sample"] == sample]
    if len(df_info.index) == 0:
        logger.error(f"{sample} is not in the SV calling INFO table")

    # read force calling table
    logger.info(f"Read force calling INFO of sample {sample} from {args.force_call_info}")
    df_force_call = parse_sv_info(args.force_call_info)
    df_force_call = df_force_call.loc[df_force_call["Sample"] == sample]
    if len(df_force_call.index) == 0:
        logger.error(f"{sample} is not in the force calling INFO table")

    # read method table
    logger.info(f"Read method table from {args.method_table}")
    sample_method = method_table(args.method_table, sample)

    # parse genotyping method
    gt_method = parse_gt_method(args.genotyping_method)

    # prepare output vcf
    new_header = add_header(invcf.header, methods=sample_method.methods)
    outvcf = pysam.VariantFile(args.outvcf, "w", header=new_header)

    # genotype sample
    genotype_vcf(invcf=invcf, 
                 outvcf=outvcf, 
                 sample=sample, 
                 sample_method=sample_method,
                 df_info=df_info, 
                 df_force_call=df_force_call,
                 min_depth=args.min_dp, 
                 count_dp=args.count_missing_dp,
                 homo_freq=args.homozygous_freq, 
                 hetero_freq=args.heterozygous_freq,
                 sample_col=sample_col,
                 gt_method=gt_method)

    sample_method.close()
    invcf.close()
    outvcf.close()
    
    logger.info(f"Write output to {args.outvcf}")
    if args.outvcf.endswith("vcf.gz"):
        pysam.bcftools.index("-t", args.outvcf)

