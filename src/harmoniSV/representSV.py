#!/usr/bin/env python3

# Select representative SV after merging
# Created: 18/8/2022
# Author: Han Cao

import pysam
import pysam.bcftools
import argparse
import pandas as pd
import logging
import os

from utils import IdGenerator, read_vcf, OutVcf, vcf_to_df, read_manifest, manifest_to_df, ProgressLogger, parse_cmdargs

# parse arguments
parser = argparse.ArgumentParser(prog="harmonisv represent",
                                 description="Select the representative SV from merged SVs",
                                 add_help=False)
io_arg = parser.add_argument_group('Input/Output arguments')
io_arg.add_argument("-i", "--invcf", metavar="VCF", type=str, required=False,
                    help="input vcf")
io_arg.add_argument("-f", "--file-list", metavar="TSV", type=str, required=False,
                    help="list of VCF files, one VCF per line")
io_arg.add_argument("-o", "--outvcf", metavar="VCF", type=str, required=True,
                    help="output vcf")
io_arg.add_argument("--merge", metavar="FILE", type=str, required=True,
                    help="SV merging results, each line is comma-separated list of IDs of merged SVs, or the ID of unique SVs")
io_arg.add_argument("-r", "--region", metavar="chr", type=str, required=False,
                    help="genomic region to work on (require vcf with index)")

represent_arg = parser.add_argument_group('Representative SV selection arguments')
represent_arg.add_argument("--id-prefix", metavar="PREFIX", type=str, required=False,
                           help="Rename output SV ID as PREFIX.SVTYPE.NUMBER")
represent_arg.add_argument("--by-max", metavar="TAG", type=str, required=False,
                           help="Select representative SV by maximum INFO/TAG value")
represent_arg.add_argument("--by-freq", action="store_true", default=False,
                           help="Select representative SV by frequency of POS, SVLEN, if >1 SVs have same frequency, select the one closest to average POS, SVLEN")
represent_arg.add_argument("--save-id", action="store_true", default=False, required=False,
                           help="Save original IDs to INFO/ID_LIST")
represent_arg.add_argument("--keep-hom-ref", action="store_true", default=False, required=False,
                           help="Keep SVs with AC=0 before selecting representative SV (Default: False)")
represent_arg.add_argument("--min-len-input", metavar="N", type=int, required=False, default=35,
                           help="Remove SVLEN < min_len_input before selecting representative SV (Default: 35)")
represent_arg.add_argument("--min-len-output", metavar="N", type=int, required=False, default=50,
                           help="Remove SVLEN < min_len_output from output (Default: 50)")

optional_arg = parser.add_argument_group('optional arguments')
optional_arg.add_argument("-h", "--help", action='help', help='show this help message and exit')



class SVCleaner:
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

        self.represent_sv = self.unique_call.copy()                               # representative SVs, initialize with unique calls
        self.dropped_sv = set()                                                        # dropped SVs, initialize with empty set

        # remove SVs with AC=0 from initial representative SVs
        if not keep_hom_ref:
            non_ref_call =  set(df_info[df_info["AC"] > 0].index)
            sv_filter_ac = self.represent_sv.difference(non_ref_call)
            self.represent_sv = self.represent_sv.intersection(non_ref_call)
            self.dropped_sv = self.dropped_sv.union(sv_filter_ac)
            self.n_filter_ac = len(sv_filter_ac)
        else:
            self.n_filter_ac = 0

        # remove SVs with SVLEN < min_len_output from initial representative SVs
        if min_len_output > 0:
            large_call = set(df_info[df_info["SVLEN"].abs() >= min_len_output].index)
            sv_filter_len = self.represent_sv.difference(large_call)
            self.represent_sv = self.represent_sv.intersection(large_call)
            self.dropped_sv = self.dropped_sv.union(sv_filter_len)
            self.n_filter_len = len(sv_filter_len)
        else:
            self.n_filter_len = 0

        # report initial status
        self.logger.info(f"No. of total SV: {len(self.all_sv)}")
        self.logger.info(f"No. of unique SV: {len(self.unique_call)}")
        self.logger.info(f"No. of unique SV selected as representative SV: {len(self.represent_sv)}")
        self.logger.info(f"No. of merging records: {len(self.merge_record)}")
    
    def add_representative(self, sv: str) -> None:
        """ Add 1 SV to representative SV set """
        self.represent_sv.add(sv)
    
    def union_representative(self, sv: set) -> None:
        """ Union SVs with representative SV set """
        self.represent_sv = self.represent_sv.union(sv)

    def add_dropped(self, sv: str) -> None:
        """ Add 1 SV to representative SV set """
        self.dropped_sv.add(sv)

    def union_dropped(self, sv: set) -> None:
        """ Union SVs with representative SV set """
        self.dropped_sv = self.dropped_sv.union(sv)

    def add_map(self, representative_sv: str, merged_sv: list) -> None:
        """ Add mapping from representative SV to merged SVs """
        self.sv_map[representative_sv] = merged_sv

    def get_map(self, representative_sv: str) -> list:
        """ Get mapping result from representative SV to all SVs """
        return self.sv_map.get(representative_sv, [representative_sv])
    
    def initialize(self) -> None:
        """ Check if all representative and dropped SVs are found """
        if self.all_sv == self.represent_sv.union(self.dropped_sv):
            self.logger.info(f"No. of representative SVs filter by homozygous ref call: {self.n_filter_ac}")
            self.logger.info(f"No. of representative SVs filter by minimum length: {self.n_filter_len}")
            self.logger.info(f"{len(self.represent_sv)} representative SVs are kept")
        else:
            self.logger.error(f"All SV ({len(self.all_sv)}) is not equal to representative SV ({len(self.represent_sv)}) + dropped SV ({len(self.dropped_sv)})")
            raise SystemExit()
        
        overlap_sv = self.represent_sv.intersection(self.dropped_sv)
        if len(overlap_sv) == 0:
            self.ready = True
        else:
            for overlap_first in overlap_sv:
                break
            self.logger.error(f"{len(overlap_sv)} representative SVs are overlapped with dropped SVs: e.g., {overlap_first}")
            raise SystemExit()
        
        if len(self.represent_sv) >= len(self.dropped_sv):
            self.find_dropped = True
        else:
            self.find_dropped = False
        
        # set active set to work on
        self.represent_sv_active = self.represent_sv.copy()
        self.dropped_sv_active = self.dropped_sv.copy()
        
        self.logger.info("Finish initialization")
    
    def is_representative(self, id: str) -> bool:
        """ check if SV is representative """
        if self.ready:
            if self.find_dropped:
                return id not in self.dropped_sv_active
            else:
                return id in self.represent_sv_active
        
        else:
            self.logger.error("Error: Please initialize before check if SV is representative")
            raise SystemExit()
        
    def subset_active(self, subset: set) -> None:
        """ Subset active represent/dropped SVs """
        if self.ready:
            self.represent_sv_active = self.represent_sv.intersection(subset)
            self.dropped_sv_active = self.dropped_sv.intersection(subset)
        else:
            self.logger.error("Error: Please initialize before subset active SVs")
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


def find_by_max(df_merged: pd.DataFrame, tag: str, keep_hom_ref: bool, min_len_input: int=0, min_len_output: int=0) -> str:
    """ Select representative SV by max tag value """

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


def find_by_freq(df_merged: pd.DataFrame, keep_hom_ref: bool, min_len_input: int=0, min_len_output: int=0) -> str:
    """ Select representative SV by most frequenct SV """

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


def find_represent(df_info: pd.DataFrame, sv_merge: list, args: argparse.Namespace) -> SVCleaner:
    """ Find representative SVs from merged SVs """
    
    # create cleaner
    invcf_cleaner = SVCleaner(sv_merge, df_info, args.keep_hom_ref, args.min_len_output)

    # find all representative SVs
    df_info_merged = df_info.loc[list(invcf_cleaner.merged_sv)].copy() # only keep info for merged SVs
    progress_rep = ProgressLogger(item='SV merging records', verbosity=2000)

    for _, merged_sv in enumerate(invcf_cleaner.merge_record):
        # find representative SV
        df_merged = df_info_merged.loc[merged_sv].copy()
        if args.by_freq:
            representative_sv = find_by_freq(df_merged, args.keep_hom_ref, args.min_len_input, args.min_len_output)        
        elif args.by_max is not None:
            representative_sv = find_by_max(df_merged, args.by_max, args.keep_hom_ref, args.min_len_input, args.min_len_output)
        
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
            if args.save_id:
                invcf_cleaner.add_map(representative_sv, merged_sv)

        # update progress
        progress_rep.log()

    progress_rep.finish()
    invcf_cleaner.initialize()

    return invcf_cleaner


def outvcf_worker(invcf: pysam.VariantFile, outvcf: pysam.VariantFile, region: str, 
                  invcf_cleaner: SVCleaner, outvar_id: IdGenerator, save_id: bool, progress: ProgressLogger) -> None:
    """ Find representative SVs in 1 input VCF """
    
    # add headers
    invcf.header.info.add("REPRESENT_SV", "1", "String", "Orignal representative SV ID")
    invcf.header.info.add("ID_LIST", "1", "String", "Orignal ID list")

    # extract working region
    if region is not None:
        invcf_iter = invcf.fetch(region=region)
    else:
        invcf_iter = invcf.fetch()
    
    for variant in invcf_iter:
        # output representative SVs
        if invcf_cleaner.is_representative(variant.id):
            new_variant = variant.copy()
            new_variant.info["REPRESENT_SV"] = variant.id

            if outvar_id is not None:
                new_variant.id = outvar_id.get(variant.info["SVTYPE"])

            if save_id:
                new_variant.info["ID_LIST"] = ",".join(invcf_cleaner.get_map(variant.id))

            outvcf.write(new_variant)
            progress.log()
    


def representSV_main(cmdargs) -> None:
    """Main function of representSV"""
    logger = logging.getLogger(__name__)

    args = parse_cmdargs(parser, cmdargs)

    # parse info to be extracted
    info = ["SVTYPE", "SVLEN", "AC"]
    if args.by_max is not None:
        info.append(args.by_max)

    # read input files
    if args.file_list is not None:
        logger.info(f"Read VCF file list: {args.file_list}")
        df_manifest = read_manifest(args.file_list, header=None, default_info=info)
    elif args.invcf is not None:
        logger.info(f"Read input vcf: {args.invcf}")
        df_manifest = pd.DataFrame({'file': [args.invcf], 'info': [info]})
    else:
        logger.error("Error: Please specify either --invcf or --manifest")
        raise SystemExit()
    
    # parse VCFs to df
    logger.info("Parse input VCFs to dataframe")
    df_info = manifest_to_df(df_manifest, args.region)
    df_info.set_index("ID", inplace=True)

    # read SV merging results
    logger.info(f"Read SV merging results: {args.merge}")
    sv_merge = [line.strip() for line in open(args.merge, "r")]

    # extract region
    if args.region is not None:
        logger.info(f"Extract SV merging records within the region: {args.region}")
        keep_sv = set(df_info.index)
        sv_merge = filter_merge_list(sv_merge, keep_sv)

    # find representative SVs
    invcf_cleaner = find_represent(df_info, sv_merge, args)

    # add headers
    with pysam.VariantFile(df_manifest['file'].values[0]) as f:
        new_header = f.header
    new_header.info.add("REPRESENT_SV", "1", "String", "Orignal representative SV ID")
    new_header.info.add("ID_LIST", "1", "String", "Orignal ID list")

    # create output vcf ID generater
    if args.id_prefix is not None:
        outvar_id = IdGenerator(args.id_prefix, args.region)
    else:
        outvar_id = None
    
    # set up logger
    progress_out = ProgressLogger(item='representative SVs', verbosity=2000)

    # if one input VCF, write to output VCF and return
    if df_manifest.shape[0] == 1:
        invcf = read_vcf(df_manifest['file'].values[0], slient=True)
        outvcf = OutVcf(args.outvcf, new_header)
        outvcf_worker(invcf, outvcf, args.region, invcf_cleaner, outvar_id, args.save_id, progress_out)
        invcf.close()
        progress_out.finish()
        outvcf.close()
        
        return None
    
    # if multiple input VCFs, write to intermediate VCF
    else:
        outvcf_tmp = OutVcf(f'{args.outvcf}.tmp', new_header)

    # work on each input VCFs
    for idx, row in df_manifest.iterrows():
        # subset SVs from the current invcf
        invcf_set = set(df_info.loc[df_info["file_idx"] == idx].index)
        invcf_cleaner.subset_active(invcf_set)

        # read and work on invcf to write representative SVs
        file_invcf = row['file']
        invcf = read_vcf(file_invcf, slient=True)
        outvcf_worker(invcf, outvcf_tmp, args.region, invcf_cleaner, None, args.save_id, progress_out)
        invcf.close()
    
    progress_out.finish()
    outvcf_tmp.close()

    # if multiple input VCFs, sort output vcf
    logger.info(f"Sort VCF")
    pysam.bcftools.sort("-o", f'{args.outvcf}.tmp.sort', f'{args.outvcf}.tmp', catch_stdout=False)
    
    # write final vcf
    sorted_vcf = read_vcf(f'{args.outvcf}.tmp.sort', slient=True)
    outvcf = OutVcf(args.outvcf, sorted_vcf.header)

    for variant in sorted_vcf:
        if outvar_id is not None:
            variant.id = outvar_id.get(variant.info["SVTYPE"])
        outvcf.write(variant)
    
    sorted_vcf.close()
    outvcf.close()

    # remove tmp files
    os.remove(f'{args.outvcf}.tmp')
    os.remove(f'{args.outvcf}.tmp.sort')
    
    logger.info("Done")
