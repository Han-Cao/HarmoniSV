#!/usr/bin/env python3

# Calculate genotype concordance between two VCFs
# Created: 8/2/2023
# Author: Han Cao

import argparse
import logging

import pysam

from utils import read_vcf, parse_cmdargs

# parse arguments
parser = argparse.ArgumentParser(prog="harmonisv concordance",
                                 description="Calculate genotype concordance between two VCFs", 
                                 add_help=False)
required = parser.add_argument_group('required arguments')
required.add_argument("-i", "--invcf", metavar="vcf", type=str, required=True,
                      help="input vcf file")
required.add_argument("-o", "--outvcf", metavar="vcf", type=str, required=True,
                      help="ouput vcf file")
required.add_argument("-r", "--ref", metavar="vcf", type=str, required=True,
                      help="reference vcf file to compare with")

optional_arg = parser.add_argument_group('optional arguments')
optional_arg.add_argument("--map", metavar="file", type=str, default=None,
                          help="tab-delimited file to map input vcf samples (1st col) to reference vcf samples (2nd col)")
optional_arg.add_argument("-h", "--help", action='help', help='show this help message and exit')


class ConcordanceCounter:
    """Accumulator to store concordance statistics"""
    def __init__(self):
        self.total = 0
        self.concordant = 0

    def add(self, n: int):
        """Add concordance statistics"""
        self.total += 1
        self.concordant += n

    def rate(self) -> float:
        """Calculate concordance rate"""
        return self.concordant / self.total



GT_DICT = {(None, None): -1,
           (None, ): -1,
           (0, 0): 0,
           (0, 1): 1,
           (1, 0): 1,
           (1, 1): 2}

def gt_concordance(gt_in: tuple, gt_ref: tuple) -> tuple:
    """Compare genotype concordance of two genotypes, return concordance by genotype and existence (i.e., 0/1 = 1/1)"""
    ac_in = GT_DICT[gt_in]
    ac_ref = GT_DICT[gt_ref]
    # if reference genotype is missing, return None and skip
    if ac_ref == -1:
        return (None, None)
    # concordance by genotype
    elif ac_in == ac_ref:
        return (1, 1)
    # if input gneotype is missing, return 0 for both
    elif ac_in == -1:
        return (0, 0)
    # concordance by existence only
    elif (ac_in > 0) and (ac_ref > 0):
        return (0, 1)
    # discordance for both
    else:
        return (0, 0)


def variant_concordance(invar: pysam.VariantRecord, 
                        refvar: pysam.VariantRecord, 
                        sample_map: dict,
                        counter_gt: ConcordanceCounter,
                        counter_exist: ConcordanceCounter) -> tuple:
    """Compare genotype concordance of two variants, return concordance by genotype and existence (i.e., 0/1 = 1/1)"""
    con_gt = 0
    con_exist = 0
    n_sample = len(sample_map)
    for sample1, sample2 in sample_map.items():
        gt_in = invar.samples[sample1]['GT']
        gt_ref = refvar.samples[sample2]['GT']
        con_gt_i, con_exist_i = gt_concordance(gt_in, gt_ref)

        # if reference genotype is missing, skip
        if con_gt_i is None:
            n_sample -= 1
            continue

        con_gt += con_gt_i
        con_exist += con_exist_i
        counter_gt.add(con_gt_i)
        counter_exist.add(con_exist_i)
    
    # if no sample can be compared, return None
    if n_sample == 0:
        return (None, None)
    else:
        return (con_gt / n_sample, con_exist / n_sample)


def compare_pos(var1: pysam.VariantRecord, var2: pysam.VariantRecord) -> str:
    """Compare if two variants are identical on the same chromosome"""
    # return when two variants are matched
    if (var1.chrom == var2.chrom) and (var1.pos == var2.pos) and (var1.alleles == var2.alleles):
        return "match"
    elif var1.pos < var2.pos:
        return "var1_before_var2"
    elif var1.pos > var2.pos:
        return "var1_after_var2"
    else:
        return "different_allele"


def multiallele_concordance(invar: pysam.VariantRecord,
                            refvcf: pysam.VariantFile,
                            outvcf: pysam.VariantFile,
                            sample_map: dict,
                            counter_gt: ConcordanceCounter,
                            counter_exist: ConcordanceCounter) -> None:
    """Handle multiallelic variants"""
    outvar = invar.copy()

    # loop all alleles in reference vcf
    position = f'{invar.chrom}:{invar.pos}-{invar.pos}'
    for refvar in refvcf.fetch(region=position):
        if invar.alleles == refvar.alleles:
            con_gt, con_exist = variant_concordance(invar, refvar, sample_map, counter_gt, counter_exist)
            outvar.info['CON_GT'] = con_gt
            outvar.info['CON_EXIST'] = con_exist
            outvcf.write(outvar)
            return None
    
    # if no allele can match, set concordance to missing
    outvar.info['CON_GT'] = None
    outvar.info['CON_EXIST'] = None
    outvcf.write(outvar)
    return None


def compare_chr(invcf: pysam.VariantFile, 
                refvcf: pysam.VariantFile, 
                outvcf: pysam.VariantFile,
                chr: str, 
                sample_map: dict,
                counter_gt: ConcordanceCounter,
                counter_exist: ConcordanceCounter,
                verbosity: int=2000):
    """Compare all variants on the same chromosome"""
    # initialize
    invcf_iter = invcf.fetch(region=chr)
    refvcf_iter = refvcf.fetch(region=chr)

    logger = logging.getLogger(__name__)

    # set flag to indicate if refvcf need next()
    # when invcf is before refvcf, set flag to False to make refvcf stop at current position
    ref_flag = True 

    # loop through all variants of invcf
    for i, invar in enumerate(invcf_iter):
        if i % verbosity == 0:
            logger.info(f"Processing variant at {invar.chrom}:{invar.pos} ({invar.id})...")

        outvar = invar.copy()
        
        # loop through all variants of refvcf
        loop_count = 0
        while True:
            loop_count += 1
            if loop_count > 10000:
                logger.error(f"Looping too many times at invcf: {invar.chrom}:{invar.pos} ({invar.id})")
                break
            # go to next variant of refvcf
            if ref_flag:
                try:
                    refvar = next(refvcf_iter)
                # if refvcf is at the end, set output concordance as None
                # iterate to next variant of invcf
                except StopIteration:
                    outvar.info['CON_GT'] = None
                    outvar.info['CON_EXIST'] = None
                    outvcf.write(outvar)
                    break

            # compare two variants
            comp = compare_pos(invar, refvar)
            if comp == "match":
                con_gt, con_exist = variant_concordance(invar, refvar, sample_map, counter_gt, counter_exist)
                outvar.info['CON_GT'] = con_gt
                outvar.info['CON_EXIST'] = con_exist
                outvcf.write(outvar)
                ref_flag = True
                break
            elif comp == "var1_before_var2":
                outvar.info['CON_GT'] = None
                outvar.info['CON_EXIST'] = None
                outvcf.write(outvar)
                ref_flag = False
                break
            elif comp == "var1_after_var2":
                ref_flag = True
                continue
            elif comp == "different_allele":
                ref_flag = False
                multiallele_concordance(invar, refvcf, outvcf, sample_map, counter_gt, counter_exist)
                break



def add_header(header: pysam.VariantHeader) -> pysam.VariantHeader:
    """Add header to output VCF"""
    header.info.add('CON_GT', '1', 'Float', 'Percentage of samples with concordant genotypes')
    header.info.add('CON_EXIST', '1', 'Float', 'Percentage of samples with concordant existence of variants')
    return header


def read_map(mapfile: str) -> dict:
    """Read sample map file"""
    sample_map = {}
    with open(mapfile, 'r') as f:
        for line in f:
            sample1, sample2 = line.strip().split('\t')
            sample_map[sample1] = sample2
    return sample_map


def concordanceVCF_main(cmdargs):
    """Main function to evaluate genotype concordance of two VCF files"""
    args = parse_cmdargs(parser, cmdargs)

    logger = logging.getLogger(__name__)
    logger.info('Calculate genotype concordance between two VCFs')

    # read input
    invcf = read_vcf(args.invcf)
    refvcf = read_vcf(args.ref)

    # generate output
    new_header = invcf.header
    new_header = add_header(new_header)
    outvcf = pysam.VariantFile(args.outvcf, 'w', header=new_header)

    # setup counters
    counter_gt = ConcordanceCounter()
    counter_exist = ConcordanceCounter()

    # get all contigs
    contig_invcf = invcf.header.contigs.keys()
    contig_refvcf = refvcf.header.contigs.keys()

    # get sample map
    if args.map is not None:
        sample_map = read_map(args.map)
    else:
        sample_invcf = list(invcf.header.samples)
        sample_refvcf = list(refvcf.header.samples)
        sample_overlap = [x for x in sample_invcf if x in sample_refvcf]
        sample_map = dict(zip(sample_overlap, sample_overlap))
    
    logger.info(f'Find {len(sample_map)} samples in both VCF files to compare.')

    # loop through all contigs
    for chr in contig_invcf:
        if chr in contig_refvcf:
            compare_chr(invcf, refvcf, outvcf, chr, sample_map, counter_gt, counter_exist)
        else:
            logger.warning(f'Contig {chr} is not found in reference VCF file. Set concordance to missing.')
            for variant in invcf.fetch(region=chr):
                variant.info['CON_GT'] = None
                variant.info['CON_EXIST'] = None
                outvcf.write(variant)
    
    # summarise concordance results
    logger.info(f'Overall concordance of genotypes: {counter_gt.rate()}')
    logger.info(f'Overall concordance of existence: {counter_exist.rate()}')
    

    # close files
    invcf.close()
    refvcf.close()
    outvcf.close()

    logger.info('Done')

    return 0