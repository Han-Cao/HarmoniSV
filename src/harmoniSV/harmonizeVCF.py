#!/usr/bin/env python3

# Harmonize SV VCFs from across samples and SV calling methods
# Created: 12/8/2022
# Author: Han Cao

import pysam
import pysam.bcftools
import argparse
import functools
import logging
import os

from utils import read_vcf, IdGenerator, SVTYPE_ALLOW, parse_cmdargs

# parse arguments
parser = argparse.ArgumentParser(prog="harmonisv harmonize",
                                 description="Harmonize SV VCFs from across samples and SV calling methods",
                                 add_help=False)

required = parser.add_argument_group('Required arguments')
required.add_argument("-i", "--invcf", metavar="VCF", type=str, required=True,
                      help="input VCF")
required.add_argument("-o", "--outvcf", metavar="VCF", type=str, required=True,
                      help="output VCF")

vcf_info_arg = parser.add_argument_group('VCF INFO manipulation')
vcf_info_arg.add_argument("--info", metavar="TAG", type=str, required=False,
                          help="Comma separated INFO tags to extract. INFO tags can be renamed by NEW=OLD, from high prioirty to low priority, e.g., 'NEW=OLD1,NEW=OLD2' means if OLD1 is present, use OLD1, otherwise use OLD2.")
vcf_info_arg.add_argument("--info-sum", metavar="TAG", type=str, required=False,
                          help="Comma separated INFO tags to extract and sum. Old tags with the same new tag name will be summed up, e.g., 'NEW=OLD1,NEW=OLD2' means 'INFO/NEW = INFO/OLD1 + INFO/OLD2'")
vcf_info_arg.add_argument("--format-to-info", metavar="TAG", type=str, required=False,
                          help="Comma separated FORMAT tags to sum across samples and add to INFO, from high prioirty to low priority, e.g., DP=DP means 'INFO/DP = sum(FORMAT/DP)'.")
vcf_info_arg.add_argument("--format-to-info-sum", metavar="TAG", type=str, required=False,
                          help="Comma separated FORMAT tags to sum across samples and tags and add to INFO, e.g., 'DP=DR,DP=DV' means 'INFO/DP = sum(FORMAT/DR) + sum(FORMAT/DV)'.")
vcf_info_arg.add_argument("--info-to-alt", metavar="TAG", type=str, required=False,
                          help="Comma separated INFO tags to fill in ALT, from high prioirty to low priority. This is useful if insertion sequence stored in INFO.")
# vcf_info_arg.add_argument("--sum", action="store_true", default=False,
#                           help="Change merge logic to sum, all tags must exist for all records, e.g., '--format-to-info DP=DR,DP=DV --sum' means 'INFO/DP = sum(FORMAT/DR + FORMAT/DV)'")
vcf_info_arg.add_argument("--keep-old", action="store_true", default=False,
                          help="Keep original INFO tags after renaming (default: False)")
vcf_info_arg.add_argument("--keep-all", action="store_true", default=False,
                          help="Keep all original INFO tags (default: False)")
vcf_info_arg.add_argument("--no-AC", action="store_true", default=False,
                          help="Disable automatically add AC to INFO (default: False)")

vcf_header_arg = parser.add_argument_group('VCF header manipulation')
vcf_header_arg.add_argument("--header", metavar="FILE", type=str, required=False,
                            help="New VCF header file to replace the header of input VCF.")
vcf_header_arg.add_argument("--header-str", metavar="string", type=str, required=False,
                            help="Semicolon separated INFO header string added to new header (metadata separated by comma), e.g., 'DP,1,Integar,Sequencing depth;AF,1,Float,Allele frequency'")

sv_arg = parser.add_argument_group('Structural variation format')
sv_arg.add_argument("--id-prefix", metavar="PREFIX", type=str, required=False,
                    help="Rename SV ID to PREFIX.raw_ID. Final ID should be Sample.Aligner.Caller.unique_id for downstream analysis")
sv_arg.add_argument("--rename-id", action="store_true", default=False,
                    help="Rename SV ID to PREFIX.SVTYPE.No., must use with --id-prefix")
sv_arg.add_argument("--svtype", metavar="SVTYPE", type=str, required=False, default="SVTYPE",
                    help="INFO tag stores the SVTYPE (default: SVTYPE). Will rename it to SVTYPE if not.")
sv_arg.add_argument("--INS", metavar="INS", type=str, required=False, default="INS",
                    help="Comma separated SVTYPE string for insertions, will be nomalized as INS (default: INS)")
sv_arg.add_argument("--DEL", metavar="DEL", type=str, required=False, default="DEL",
                    help="Comma separated SVTYPE string for deletions, will be nomalized as DEL (default: DEL)")
sv_arg.add_argument("--DUP", metavar="DUP", type=str, required=False, default="DUP",
                    help="Comma separated SVTYPE string for duplications, will be nomalized as DUP (default: DUP)")
sv_arg.add_argument("--INV", metavar="INV", type=str, required=False, default="INV",
                    help="Comma separated SVTYPE string for inversions, will be nomalized as INV (default: INV)")
sv_arg.add_argument("--CNV", metavar="CNV", type=str, required=False, default="CNV",
                    help="Comma separated SVTYPE string for copy number variations, will be nomalized as CNV (default: CNV)")
sv_arg.add_argument("--BND", metavar="BND", type=str, required=False, default="BND",
                    help="Comma separated SVTYPE string for breakends, will be nomalized as BND (default: BND)")                 
sv_arg.add_argument("--no-check", action="store_true", default=False,
                    help="Disable check and filter on ID, SVTYPE for downstream analysis (default: False)")

optional_arg = parser.add_argument_group('Optional arguments')
optional_arg.add_argument("-h", "--help", action='help', help='show this help message and exit')

# parse args.info pattern string
class TagParser:
    def __init__(self, pattern_str: str) -> None:
        self.pattern_str = pattern_str
        self.pattern_dict = self.parse_pattern_str(pattern_str)
        self.new_tags = self.pattern_dict.keys()
        self.old_tags = [tag for old_tag_list in self.pattern_dict.values() for tag in old_tag_list]

    def parse_pattern_str(self, pattern_str: str) -> dict:
        """ Parse pattern string """
        pattern_dict = {}
        for x in pattern_str.split(","):
            if "=" in x:
                new_tag, old_tag = x.split("=")
            else:
                new_tag = x
                old_tag = x

            # add or rename
            if new_tag not in pattern_dict:
                pattern_dict[new_tag] = [old_tag]
            # merge multiple old tags into one new tag
            else:
                pattern_dict[new_tag].append(old_tag)

        return pattern_dict
    
    def add_tag(self, new_tag: str, old_tag:str) -> None:
        """ Add new tag and map to old tags """
        if new_tag not in self.pattern_dict:
            self.pattern_dict[new_tag] = [old_tag]
        else:
            self.pattern_dict[new_tag].append(old_tag)

    def map_tag(self, new_tag:str) -> list:
        """ Map new tag to old tag list """
        return self.pattern_dict[new_tag]


class SVTypeMapper:
    """ Map SVTYPE strings"""
    def __init__(self, args: argparse.Namespace) -> None:
        self.svtype_dict = self.parse_svtype_str(args)
    
    def parse_svtype_str(self, args: argparse.Namespace) -> dict:
        """ Parse SVTYPE string """
        svtype_dict = {}

        if (    
                args.INS == "INS" and 
                args.DEL == "DEL" and
                args.DUP == "DUP" and
                args.INV == "INV" and
                args.CNV == "CNV" and
                args.BND == "BND"
            ):
            self.require_map = False
            return None
        else:
            self.require_map = True

            for x in args.INS.split(","):
                svtype_dict[x] = "INS"
            for x in args.DEL.split(","):
                svtype_dict[x] = "DEL"
            for x in args.DUP.split(","):
                svtype_dict[x] = "DUP"
            for x in args.INV.split(","):
                svtype_dict[x] = "INV"
            for x in args.CNV.split(","):
                svtype_dict[x] = "CNV"
            for x in args.BND.split(","):
                svtype_dict[x] = "BND"
            
            return svtype_dict
    
    def map_svtype(self, svtype: str) -> str:
        """ Map SVTYPE string """
        return self.svtype_dict[svtype]


def parse_tags(args) -> None:
    """ Parse tags to be extracted """
    info = args.info
    if info is None:
        info = f"SVTYPE={args.svtype}"
    info_parser = TagParser(info)

    if args.info_sum is not None:
        info_sum_parser = TagParser(args.info_sum)
    else:
        info_sum_parser = None

    if args.format_to_info is not None:
        format_parser = TagParser(args.format_to_info)
    else:
        format_parser = None
    
    if args.format_to_info_sum is not None:
        format_sum_parser = TagParser(args.format_to_info_sum)
    else:
        format_sum_parser = None

    if args.info_to_alt is not None:
        alt_parser = TagParser(args.info_to_alt)
    else:
        alt_parser = None

    return info_parser, info_sum_parser, format_parser, format_sum_parser, alt_parser


def sum_tag(a, b):
    """ sum two tuple / number """
    if isinstance(a, tuple) and isinstance(b, tuple):
        return tuple(map(sum, zip(a, b)))
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return a + b
    else:
        raise SystemExit(f"Error: {a} and {b} can not be added, only support tuple + tuple, number + number")

def is_tag_none(x) -> bool:
    """ Check if tuple or number is None"""
    if isinstance(x, tuple):
        return all([x is None for x in x])
    else:
        return x is None


def sum_info(info: pysam.VariantRecordInfo, tags: list):
    """ Sum of tag values in INFO """
    tag_values = [info[tag] for tag in tags]
    return functools.reduce(sum_tag, tag_values)


def sum_format(samples: pysam.VariantRecordSamples, tag: str):
    """ Sum of FORMAT tag values across samples """
    sample_list = list(samples)
    sample_values = [samples[x][tag] for x in sample_list if tag in samples[x]]
    # remove missing values
    # TODO: This don't support partially missing values like (1,None), not sure if this is a problem in real data
    sample_values = [x for x in sample_values if not is_tag_none(x)]
    if len(sample_values) > 0:
        return functools.reduce(sum_tag, sample_values)
    else:
        return None


# Only support bi-allelic
AC_DICT = {(None, None): 0,
           (0, 0): 0,
           (0, 1): 1,
           (1, 0): 1,
           (1, 1): 2}


def sum_AC(samples, AC_dict=AC_DICT):
    """ Sum GT to get AC """
    AC = 0
    for id in samples:
        AC += AC_dict[samples[id]["GT"]]
    
    return AC


def reheader(file_invcf: str, output_prefix: str, file_header: str) -> str:
    file_out = f"{output_prefix}.reheader.tmp"

    pysam.bcftools.reheader("-h", file_header,
                            "-o", file_out,
                            file_invcf,
                            catch_stdout=False)
    
    return file_out


def add_header(header: pysam.VariantHeader, info_parser: TagParser, header_str: str) -> pysam.VariantHeader:
    """ Create new vcf header """
    logger = logging.getLogger(__name__)

    if header_str is not None:
        for header_line in header_str.split(";"):
            header_line = header_line.split(",")
            if len(header_line) != 4:
                logger.error(f"Header line {header_line} should have 4 fields")
                raise SystemExit()

            if header_line[0] not in header.info:
                header.info.add(header_line[0], header_line[1], header_line[2], header_line[3])

    for new_tag in info_parser.new_tags:
        if new_tag in header.info:
            continue
        else:
            old_tag = info_parser.map_tag(new_tag)[0]
            old_header = header.info[old_tag]
            header.info.add(new_tag, old_header.number, old_header.type, old_header.description)
    # add mandatory tags if not exist
    if "SVTYPE" not in header.info:
        header.info.add("SVTYPE", "1", "String", "Type of structural variantion")
    # add INFO required for downstream analysis if missing
    if "AC" not in header.info:
        header.info.add("AC", "A", "Integer", "Allele count for each ALT allele")
    if "DP" not in header.info:
        header.info.add("DP", "1", "Integer", "Read depth")
    if "ID_OLD" not in header.info:
        header.info.add("ID_OLD", "1", "String", "Original variant ID")
    
    return header


def check_header(header: pysam.VariantHeader, parser: TagParser) -> None:
    """ Check if all new tags are defined in the header"""
    logger = logging.getLogger(__name__)
    for new_tag in parser.new_tags:
        if new_tag not in header.info:
            logger.error(f"New tag {new_tag} is not defined in the header. Please specify it by --header or --header-str.")
            raise SystemExit()



def harmonizeVCF_main(cmdargs) -> None:
    """ Harnomize vcf INFO fields """
    logger = logging.getLogger(__name__)

    args = parse_cmdargs(parser, cmdargs)

    # read input
    if args.header is not None:
        file_reheader_vcf = reheader(file_invcf=args.invcf, output_prefix=args.outvcf, file_header=args.header)
        invcf = read_vcf(file_reheader_vcf, check_n=0)
    else:
        invcf = read_vcf(args.invcf, check_n=0)
    
    # parse tags to be extracted
    info_parser, info_sum_parser, format_parser, format_sum_parser, alt_parser = parse_tags(args)

    # create new VCF header, add renamed INFO headers
    new_header = add_header(invcf.header, info_parser=info_parser, header_str=args.header_str)
    # check new VCF header
    if info_sum_parser is not None:
        check_header(new_header, info_sum_parser)
    if format_parser is not None:
        check_header(new_header, format_parser)
    if format_sum_parser is not None:
        check_header(new_header, format_sum_parser)

    # log for harmonize behavior
    # TODO: determine whether and where to log
    if args.keep_all:
        logger.info("Keep all original INFO tags")
    else:
        logger.info("Drop original INFO tags not included in --info string")
    if args.keep_old:
        logger.info("Keep old INFO names for renamed INFO")

    # write new vcf
    outvcf = pysam.VariantFile(args.outvcf, "w", header=new_header)
    
    # rename ID
    ID_set = set()
    if args.rename_id:
        outvcf_id = IdGenerator(args.id_prefix)

    svtype = args.svtype

    for variant in invcf.fetch():
        # check ID and SVTYPE
        if not args.no_check:
            if svtype not in variant.info:
                logger.warning(f"{svtype} is not found in {variant.id}, skip")
                continue
            if variant.info[svtype] not in SVTYPE_ALLOW:
                logger.warning(f"SVTYPE '{variant.info[svtype]}' of {variant.id} is not allowed, skip)")
                continue

            if not args.rename_id:
                if variant.id is None:
                    logger.warning(f"ID is missing at {variant.chrom}:{variant.pos}, skip. Please use --rename to update all ID")
                    continue
                if variant.id in ID_set:
                    logger.warning(f"Duplicate ID found: {variant.id}, only keep the first one. Please use --rename to update all ID")
                    continue
            
            ID_set.add(variant.id)

        new_variant = variant.copy()
        

        # clear all old INFO
        if not args.keep_all:
            new_variant.info.clear()

            # keep old tags if specified
            if args.keep_old:
                for old_tag in info_parser.old_tags:
                    if old_tag in variant.info:
                        new_variant.info[old_tag] = variant.info[old_tag]

        # rename ID
        if args.id_prefix is not None:
            new_variant.info['ID_OLD'] = variant.id
            if args.rename_id:
                new_variant.id = outvcf_id.get(variant.info["SVTYPE"])
            else:
                new_variant.id = f"{args.id_prefix}.{variant.id}"
        
        # calculate AC
        if not args.no_AC:
            new_variant.info["AC"] = sum_AC(variant.samples)

        # add new tags from INFO
        for new_tag in info_parser.new_tags:
            old_tag_list = info_parser.map_tag(new_tag)

            # set new tag as the first old tag
            if len(old_tag_list) == 1 or not args.sum:
                for old_tag in old_tag_list:
                    if old_tag in variant.info:
                        new_variant.info[new_tag] = variant.info[old_tag]
                        break
            # sum of old tags
            else:
                new_variant.info[new_tag] = sum_info(variant.info, old_tag_list)
        
        # add new tags from FORMAT
        if args.format_to_info is not None:
            for new_tag in format_parser.new_tags:
                format_tag_list = format_parser.map_tag(new_tag)

                # set new_tag as first sum(format tag)
                if len(format_tag_list) == 1 or not args.sum:
                    for format_tag in format_tag_list:
                        tag_sample_sum = sum_format(variant.samples, format_tag)
                        if tag_sample_sum is not None:
                            new_variant.info[new_tag] = tag_sample_sum
                            break
                # sum of format tags
                else:
                    format_tag_values = [sum_format(variant.samples, format_tag) for format_tag in format_tag_list]
                    # in case all FORMAT tag are missing
                    if None not in format_tag_values:
                        new_variant.info[new_tag] = functools.reduce(sum_tag, format_tag_values)


        # add INFO to ALT
        if args.info_to_alt is not None:
            for info_tag in alt_parser.old_tags:
                if info_tag in variant.info:
                    tag_value = variant.info[info_tag]
                    if isinstance(tag_value, tuple):
                        new_variant.alts = tag_value
                    else:
                        new_variant.alts = tuple([tag_value])
                    break

        outvcf.write(new_variant)

    invcf.close()
    outvcf.close()

    if args.header is not None:
        os.remove(file_reheader_vcf)

    logger.info(f"Write output to {args.outvcf}")
    if args.outvcf.endswith("vcf.gz"):
        pysam.bcftools.index("-t", args.outvcf)

