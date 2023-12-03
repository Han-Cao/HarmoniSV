#!/usr/bin/env python3

# Harmonize VCF headers
# Created: 12/8/2022
# Author: Han Cao

import pysam
import argparse
import logging

from utils import parse_cmdargs

# parse arguments
parser = argparse.ArgumentParser(prog="harmonisv harmonize-header",
                                 description="Harmonize VCF headers",
                                 add_help=False)

io_arg = parser.add_argument_group("Input/Output arguments")
io_arg.add_argument("-i", "--invcf", required=False, type=str, metavar="VCF",
                    help="Comma-separated list of input VCF files. Duplicate headers will use the first one, including SAMPLE header. For multi-sample VCF, please make sure all input VCFs have the same SAMPLE order.")
io_arg.add_argument("-f", "--file-list", required=False, type=str, metavar="FILE_LIST",
                    help="File containing a list of input VCF files, one VCF per line. VCFs from both -i and -f will be used.")
io_arg.add_argument("-o", "--output", required=True, type=str, metavar="OUTPUT",
                    help="Output VCF header file")

optional = parser.add_argument_group("Optional arguments")
optional.add_argument("-r", "--ref-vcf", required=False, type=str, metavar="VCF",
                      help="Reference VCF file, highest priority for duplicate headers")
optional.add_argument("-h", "--help", action="help", help="Show this help message and exit")


class UniqueHeader:
    """Class to store unique VCF headers"""
    def __init__(self) -> None:
        self.header = []
        self.id_dict = {'_GENERIC': []} # key: header key, value: list of IDs
        self.header_dict = {'_GENERIC': []} # key: header key, value: list of header lines


    def add_record(self, record: pysam.VariantHeaderRecord) -> None:
        """Add new header line if not exist"""
        record_key = record.key
        record_id = record.get('ID')

        # generic headers, no header key
        if record_id is None:
            if record_key not in self.id_dict['_GENERIC']:
                self.id_dict['_GENERIC'].append(record_key)
                self.header_dict['_GENERIC'].append(str(record))
        # INFO, FORMAT, FILTER, ALT, etc.
        else:
            # add new key (e.g., INFO, FORMAT)
            if record_key not in self.id_dict:
                self.id_dict[record_key] = []
                self.header_dict[record_key] = []
            # add new header line
            if record_id not in self.id_dict[record_key]:
                self.id_dict[record_key].append(record_id)
                self.header_dict[record_key].append(str(record))
    

    def add_header(self, header: pysam.VariantHeader) -> None:
        """Add new header"""
        for record in header.records:
            self.add_record(record)


    def get_header(self) -> list:
        """Get unique header"""
        # add generic headers first
        if len(self.id_dict['_GENERIC']) > 0:
            self.header.extend(self.header_dict['_GENERIC'])
        # add other headers
        for key in self.header_dict:
            if key != '_GENERIC':
                self.header.extend(self.header_dict[key])
        return self.header


def read_files(args: argparse.Namespace) -> list:
    """Read input files"""
    invcf_list = []
    if args.invcf is not None:
        invcf_list.extend(args.invcf.split(','))
    if args.file_list is not None:
        with open(args.file_list, 'r') as f:
            invcf_list.extend([line.strip() for line in f])
    if len(invcf_list) == 0:
        raise ValueError("Please provide input VCF files")
    
    return invcf_list


def harmonizeHeader_main(cmdargs) -> None:
    """Main function of harmonizeHeader"""

    logger = logging.getLogger(__name__)
    args = parse_cmdargs(parser, cmdargs)
    new_header = UniqueHeader()
    # read input files
    invcf_list = read_files(args)

    # add headers
    sample_list = [] # list of samples shown in the last line
    # first add reference if exist
    if args.ref_vcf is not None:
        logger.info(f"Reading reference VCF file: {args.ref_vcf}")
        ref_vcf = pysam.VariantFile(args.ref_vcf)
        new_header.add_header(ref_vcf.header)
        sample_list.extend(list(ref_vcf.header.samples))
        ref_vcf.close()
        
    # add all VCFs
    for vcf_file in invcf_list:
        if args.ref_vcf is not None:
            if vcf_file == args.ref_vcf:
                continue
        logger.info(f"Reading input VCF file: {vcf_file}")
        vcf = pysam.VariantFile(vcf_file)
        new_header.add_header(vcf.header)
        # add samples
        if len(sample_list) == 0:
            sample_list.extend(list(vcf.header.samples))
        vcf.close()
    
    # add final header line
    new_header_list = new_header.get_header()
    sample_str = '\t'.join(sample_list)
    new_header_list.append(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_str}\n')

    # write output
    logger.info(f"Writing output VCF header file: {args.output}")
    with open(args.output, 'w') as f:
        f.write(''.join(new_header_list))
    logger.info("Done")