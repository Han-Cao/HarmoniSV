#!/usr/bin/env python3

# Select representative SV after merging
# Created: 18/8/2022
# Author: Han Cao

import pysam
import pysam.bcftools
import argparse
import pandas as pd
import logging

from utils import ID_generator, read_vcf, parse_cmdargs

# parse arguments
parser = argparse.ArgumentParser(prog="harmonisv represent",
                                 description="Select the representative SV from merged SVs",
                                 add_help=False)
parser.add_argument("-i", "--invcf", metavar="vcf", type=str, required=True,
                    help="input vcf")
parser.add_argument("-o", "--outvcf", metavar="vcf", type=str, required=True,
                    help="output vcf")
parser.add_argument("-r", "--region", metavar="chr", type=str, required=False,
                    help="chromosome to work on, if not specified, work on all chromosomes")
parser.add_argument("--merge", metavar="FILE", type=str, required=True,
                    help="SV merge results, each line represent a non-overlapping SV, merged SVs separated by comma in 1 line")
parser.add_argument("--sv-info", metavar="tsv", type=str, required=True,
                    help="Tab-separated file for SV INFO used for selecting representataive SV. Mandatory fields: ID(unique as index), CHR, POS, SVTYPE, SVLEN, AC (NOT MAC, as AC==0 will skip by --remove-hom-ref)")
parser.add_argument("--id-prefix", metavar="ID", type=str, required=False,
                    help="Rename output SV ID as PREFIX.SVTYPE.No")
parser.add_argument("--by-max", metavar="TAG", type=str, required=False,
                    help="Select representative SV by maximum TAG in sv-info file")
parser.add_argument("--by-freq", action="store_true", default=False,
                    help="Select representative SV by frequency of POS, SVLEN, if >1 SVs have same frequency, select the one closest to average POS, SVLEN")
parser.add_argument("--save-id", action="store_true", default=False, required=False,
                    help="Save original IDs to INFO/ID_LIST")
parser.add_argument("--keep-hom-ref", action="store_true", default=False, required=False,
                    help="Keep hom_ref SV before selecting representative SV (Default: False)")
parser.add_argument("--min-len-input", metavar="N", type=int, required=False, default=35,
                    help="Remove SVLEN < min_len_input before selecting representative SV (Default: 35)")
parser.add_argument("--min-len-output", metavar="N", type=int, required=False, default=50,
                    help="Remove SVLEN < min_len_output after selecting representative SV (Default: 50)")

optional_arg = parser.add_argument_group('optional arguments')
optional_arg.add_argument("-h", "--help", action='help', help='show this help message and exit')



class SV_cleaner:
    """ Store and check representative SVs """
    def __init__(self, merge_list: list, df_info: pd.DataFrame, keep_hom_ref: bool=False, min_len_output: int=0) -> None:
        self.merge_list = merge_list
        self.logger = logging.getLogger(__name__)

        self.merge_record = [x.split(",") for x in merge_list if "," in x]             # store merged SVs into list
        self.unique_call = set([x for x in merge_list if "," not in x])                # SVs only detected by 1 method
        self.merged_sv = [x for record in self.merge_record for x in record]           # all SVs with mering event
        self.all_sv = self.unique_call.union(self.merged_sv)                           # all SVs
        self.sv_map = {}                                                               # map representative SV to original SV IDs
        self.ready = False                                                             # state to indicate whether finish initialization

        self.representative_sv = self.unique_call.copy()                               # representative SVs, initialize with unique calls
        self.dropped_sv = set()                                                        # dropped SVs, initialize with empty set

        if not keep_hom_ref:
            non_ref_call =  set(df_info[df_info["AC"] > 0].index)
            sv_filter_ac = self.representative_sv.difference(non_ref_call)
            self.representative_sv = self.representative_sv.intersection(non_ref_call)
            self.dropped_sv = self.dropped_sv.union(sv_filter_ac)
            self.n_filter_ac = len(sv_filter_ac)
        else:
            self.n_filter_ac = 0

        if min_len_output > 0:
            large_call = set(df_info[df_info["SVLEN"].abs() >= min_len_output].index)
            sv_filter_len = self.representative_sv.difference(large_call)
            self.representative_sv = self.representative_sv.intersection(large_call)
            self.dropped_sv = self.dropped_sv.union(sv_filter_len)
            self.n_filter_len = len(sv_filter_len)
        else:
            self.n_filter_len = 0

        self.logger.info(f"No. of total SV: {len(self.all_sv)}")
        self.logger.info(f"No. of unique SV: {len(self.unique_call)}")
        self.logger.info(f"No. of unique SV selected as representative SV: {len(self.representative_sv)}")
        self.logger.info(f"No. of merging records: {len(self.merge_record)}")
    
    def add_representative(self, sv: str) -> None:
        """ Add 1 SV to representative SV set """
        self.representative_sv.add(sv)
    
    def union_representative(self, sv: set) -> None:
        """ Union SVs with representative SV set """
        self.representative_sv = self.representative_sv.union(sv)

    def add_dropped(self, sv: str) -> None:
        """ Add 1 SV to representative SV set """
        self.dropped_sv.add(sv)

    def union_dropped(self, sv: set) -> None:
        """ Union SVs with representative SV set """
        self.dropped_sv = self.dropped_sv.union(sv)

    def add_map(self, representative_sv: str, all_sv: list) -> None:
        """ Add mapping from representative SV to all SVs """
        self.sv_map[representative_sv] = all_sv

    def get_map(self, representative_sv: str) -> list:
        """ Get mapping result from representative SV to all SVs """
        return self.sv_map.get(representative_sv, [representative_sv])
    
    def initialize(self) -> None:
        """ Check if all representative and dropped SVs are found """
        if self.all_sv == self.representative_sv.union(self.dropped_sv):
            self.logger.info(f"No. of representative SVs filter by homozygous ref call: {self.n_filter_ac}")
            self.logger.info(f"No. of representative SVs filter by minimum length: {self.n_filter_len}")
            self.logger.info(f"{len(self.representative_sv)} representative SVs are kept")
        else:
            self.logger.error(f"All SV ({len(self.all_sv)}) is not equal to representative SV ({len(self.representative_sv)}) + dropped SV ({len(self.dropped_sv)})")
            raise SystemExit()
        
        overlap_sv = self.representative_sv.intersection(self.dropped_sv)
        if len(overlap_sv) == 0:
            self.ready = True
        else:
            for overlap_first in overlap_sv:
                break
            self.logger.error(f"{len(overlap_sv)} representative SVs are overlapped with dropped SVs: e.g., {overlap_first}")
            raise SystemExit()
        
        if len(self.representative_sv) >= len(self.dropped_sv):
            self.find_dropped = True
        else:
            self.find_dropped = False
        
        self.logger.info("Finish initialization")
    
    def is_representative(self, id: str) -> bool:
        """ check if SV is representative """
        if self.ready:
            if self.find_dropped:
                return id not in self.dropped_sv
            else:
                return id in self.representative_sv
        
        else:
            self.logger.error("Error: Please initialize before check if SV is representative")
            raise SystemExit()


def filter_merge_list(merge_list: list, keep_sv: set) -> list:
    """Keep merge records if any SV is in keep_sv"""
    keep_list = []
    for record in merge_list:
        for sv in record.split(','):
            if sv in keep_sv:
                keep_list.append(record)
                break
    
    return keep_list


def sv_distance(pos_1, svlen_1, pos_2, svlen_2):
    """ Fast calculate distance between 2 SVs for comparison usage """
    return (pos_2 - pos_1) ** 2 + (svlen_2 - svlen_1) ** 2


def find_by_max(df: pd.DataFrame, tag: str, idx: list, keep_hom_ref: bool, min_len_input: int=0, min_len_output: int=0) -> str:
    """ Select representative SV by max tag value """
    df_merged = df.loc[idx]

    if not keep_hom_ref:
        df_merged = df_merged[df_merged["AC"] > 0]
        # no representative SV when all merged SVs AC==0
        if len(df_merged.index) == 0:
            return "FILTER_AC0"

    if min_len_input > 0:
        df_merged = df_merged[abs(df_merged["SVLEN"]) >= min_len_input]
        # no representative SV when all merged SVs SVLEN < min_len_input
        if len(df_merged.index) == 0:
            return "FILTER_LEN"
    
    represent_sv = df_merged[tag].idxmax()

    if min_len_output > 0 and abs(df_merged.at[represent_sv, "SVLEN"]) < min_len_output:
        return "FILTER_LEN"
    else:
        return represent_sv


def find_by_freq(df: pd.DataFrame, idx: list, keep_hom_ref: bool, min_len_input: int=0, min_len_output: int=0) -> str:
    """ Select representative SV by most frequenct SV """
    df_merged = df.loc[idx]

    if not keep_hom_ref:
        df_merged = df_merged[df_merged["AC"] > 0]
        # no representative SV when all merged SVs AC==0
        if len(df_merged.index) == 0:
            return "FILTER_AC0"
    
    if min_len_input > 0:
        df_merged = df_merged[abs(df_merged["SVLEN"]) >= min_len_input]
        # no representative SV when all merged SVs SVLEN < min_len_input
        if len(df_merged.index) == 0:
            return "FILTER_LEN"
    
    # return called SV if only 1 call
    if len(df_merged.index) == 1:
        represent_sv = df_merged.index[0]        
    # select by frequency if more than 1 call
    else:
        df_merged_count = df_merged.groupby(["POS", "SVLEN"]).size().reset_index(name="count")
        df_merged_freq = df_merged_count[df_merged_count["count"] == df_merged_count["count"].max()]

        # find by most frequent POS and SVLEN
        if len(df_merged_freq.index) == 1:
            pos = df_merged_freq["POS"].values[0]
            svlen = df_merged_freq["SVLEN"].values[0]
            represent_sv = df_merged[(df_merged["POS"] == pos) & (df_merged["SVLEN"] == svlen)].index[0]
        # find by distance if fail to find most frequent POS and SVLEN
        else:
            pos_mean = df_merged["POS"].mean()
            svlen_mean = df_merged["SVLEN"].mean()
            df_merged["distance"] = df_merged.apply(lambda x: sv_distance(x.POS, x.SVLEN, pos_mean, svlen_mean), axis=1)
            represent_sv = df_merged["distance"].idxmin()
        
    # filter if less than min_len_output
    if min_len_output > 0 and abs(df_merged.at[represent_sv, "SVLEN"]) < min_len_output:
        return "FILTER_LEN"
    else:
        return represent_sv


def representSV_main(cmdargs) -> None:
    """Main function of representSV"""
    logger = logging.getLogger(__name__)

    args = parse_cmdargs(parser, cmdargs)

    # parse arguments
    file_invcf = args.invcf 
    file_outvcf = args.outvcf
    region = args.region
    file_merge = args.merge
    file_info = args.sv_info
    by_freq = args.by_freq
    by_max = args.by_max
    id_prefix = args.id_prefix
    save_id = args.save_id
    keep_hom_ref = args.keep_hom_ref
    min_len_input = args.min_len_input
    min_len_output = args.min_len_output
    verbosity = 1000

    # read input files
    invcf = read_vcf(file_invcf)
    sv_merge = [line.strip() for line in open(file_merge)]
    df_info = pd.read_csv(file_info, index_col="ID", sep="\t", na_values=[".", "NA"], 
                          dtype={'ID': str, 'CHR': str, 'SVTYPE': str})

    # extract region
    if region is not None:
        logger.info(f"Extract chromosome: {region}")
        df_info = df_info[df_info['CHR'] == region]
        keep_sv = set(df_info.index)
        sv_merge = filter_merge_list(sv_merge, keep_sv)

    # create cleaner
    invcf_cleaner = SV_cleaner(sv_merge, df_info, keep_hom_ref, min_len_output)

    # create output vcf ID generater
    if id_prefix is not None:
        outvcf_id = ID_generator(id_prefix, region)

    # find all representative SVs
    df_info_merged = df_info.loc[list(invcf_cleaner.merged_sv)] # only keep info for merged SVs

    for i, merged_sv in enumerate(invcf_cleaner.merge_record, start=1):
        # find representative SV
        if by_freq:
            representative_sv = find_by_freq(df_info_merged, merged_sv, keep_hom_ref, min_len_input, min_len_output)        
        elif by_max is not None:
            representative_sv = find_by_max(df_info_merged, by_max, merged_sv, keep_hom_ref, min_len_input, min_len_output)
        
        # update n_filter
        if representative_sv == "FILTER_AC0":
            invcf_cleaner.union_dropped(set(merged_sv))
            invcf_cleaner.n_filter_ac += 1
        elif representative_sv == "FILTER_LEN":
            invcf_cleaner.union_dropped(set(merged_sv))
            invcf_cleaner.n_filter_len += 1
        # update representative and drop SVs
        else:
            dropped_sv = set(merged_sv) - set([representative_sv])
            invcf_cleaner.add_representative(representative_sv)
            invcf_cleaner.union_dropped(dropped_sv)
            # update sv_map
            if save_id:
                invcf_cleaner.add_map(representative_sv, merged_sv)

        if i % verbosity == 0:
            logger.info(f"Finish {i} merging records")
    
    # check and initialization
    invcf_cleaner.initialize()

    # add vcf header
    new_header = invcf.header
    if save_id:
        new_header.info.add("ID_LIST", "1", "String", "Orignal ID list")
        new_header.info.add("REPRESENT_SV", "1", "String", "Orignal representative SV ID")

    # output representative SVs
    outvcf = pysam.VariantFile(file_outvcf, "w", header=new_header)
    if region is not None:
        invcf_iter = invcf.fetch(region=region)
    else:
        invcf_iter = invcf.fetch()

    for variant in invcf_iter:
        if invcf_cleaner.is_representative(variant.id):
            new_variant = variant.copy()

            if id_prefix is not None:
                new_variant.id = outvcf_id.get(variant.info["SVTYPE"])

            if save_id:
                new_variant.info["REPRESENT_SV"] = variant.id
                new_variant.info["ID_LIST"] = ",".join(invcf_cleaner.get_map(variant.id))

            outvcf.write(new_variant)
    
    invcf.close()
    outvcf.close()
    logger.info(f"Write output to {file_outvcf}")
    if file_outvcf.endswith("vcf.gz"):
        pysam.bcftools.index("-t", file_outvcf)

