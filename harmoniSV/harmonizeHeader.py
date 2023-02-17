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

required = parser.add_argument_group("required arguments")
parser.add_argument("-i", "--invcf", required=False, type=str,
                    help="Comma-separated list of input VCF files. Duplicate headers will use the first one.")
parser.add_argument("-o", "--output", required=True, type=str,
                    help="Output VCF header file")
parser.add_argument("-f", "--file-list", required=False, type=str,
                    help="File containing a list of input VCF files, one VCF per line")

optional = parser.add_argument_group("optional arguments")
optional.add_argument("-r", "--ref-vcf", required=False, type=str,
                      help="Reference VCF file, highest priority for duplicate headers")
optional.add_argument("-h", "--help", action="help", help="Show this help message and exit")


class UniqueHeader:
    """Class to store unique VCF headers"""
    def __init__(self) -> None:
        self.header = []
        self.id_dict = {'_GENERIC': []}
        self.header_dict = {'_GENERIC': []}


    def add_record(self, record: pysam.VariantHeaderRecord) -> None:
        """Add new header line if not exist"""
        record_key = record.key
        record_id = record.get('ID')

        # unstructured header
        if record_id is None:
            if record_key not in self.id_dict['_GENERIC']:
                self.id_dict['_GENERIC'].append(record_key)
                self.header_dict['_GENERIC'].append(str(record))
        else:
            if record_key not in self.id_dict:
                self.id_dict[record_key] = []
                self.header_dict[record_key] = []
            if record_id not in self.id_dict[record_key]:
                self.id_dict[record_key].append(record_id)
                self.header_dict[record_key].append(str(record))
    

    def add_header(self, header: pysam.VariantHeader) -> None:
        """Add new header"""
        for record in header.records:
            self.add_record(record)


    def get_header(self) -> list:
        """Get unique header"""
        if len(self.id_dict['_GENERIC']) > 0:
            self.header.extend(self.header_dict['_GENERIC'])
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
    sample_list = []
    if args.ref_vcf is not None:
        logger.info(f"Reading reference VCF file: {args.ref_vcf}")
        ref_vcf = pysam.VariantFile(args.ref_vcf)
        new_header.add_header(ref_vcf.header)
        sample_list.extend(list(ref_vcf.header.samples))
        ref_vcf.close()
        

    for vcf_file in invcf_list:
        if args.ref_vcf is not None:
            if vcf_file == args.ref_vcf:
                continue
        logger.info(f"Reading input VCF file: {vcf_file}")
        vcf = pysam.VariantFile(vcf_file)
        new_header.add_header(vcf.header)
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