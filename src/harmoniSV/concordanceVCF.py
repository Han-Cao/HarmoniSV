#!/usr/bin/env python3

# Calculate genotype concordance between two VCFs
# Last update: 8-Oct-2024
# Author: Han Cao

import argparse
import logging

import pysam
import numpy as np

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
required.add_argument("-r", "--refvcf", metavar="vcf", type=str, required=True,
                      help="reference vcf file to compare with")

optional_arg = parser.add_argument_group('optional arguments')
optional_arg.add_argument("--region", metavar="chr[:start-end]", type=str, default=None,
                          help="region to compare (default: all)")
optional_arg.add_argument("--map", metavar="file", type=str, default=None,
                          help="tab-delimited file to map input vcf samples (1st col) to reference vcf samples (2nd col)")
optional_arg.add_argument('--compare', metavar="alleles", type=str, default="alleles", help="Use which information to determine if 2 variants are the same: alleles (default) or id")
optional_arg.add_argument("--include-missing", action="store_true", default=False,
                          help="Include missing genotypes in --invcf for comparison (default: False)")
optional_arg.add_argument("--pass-only", action="store_true", default=False,
                          help="Only consider FILTER == PASS variants in -r VCF (default: False)")
optional_arg.add_argument("-h", "--help", action='help', help='show this help message and exit')


class ConcordanceCounter:
    """Accumulator to store concordance statistics"""
    def __init__(self):
        self.total = 0
        self.concordant = 0

    def add(self, x: float):
        """Add concordance statistics"""
        self.total += 1
        self.concordant += x

    def rate(self) -> float:
        """Calculate concordance rate"""
        return float(self.concordant / self.total)


GT_DICT = {(0, 0): 0,
           (0, 1): 1,
           (1, 0): 1,
           (1, 1): 2}


def gt2code(gt: tuple) -> int:
    """ Convert gentoype to code for comparison """
    if None in gt:
        return -1
    else:
        return GT_DICT[gt]


def variant_concordance(invar: pysam.VariantRecord, 
                        refvar: pysam.VariantRecord, 
                        in_samples: list,
                        ref_samples: list,
                        include_missing: bool,
                        counter_gt: ConcordanceCounter,
                        counter_exist: ConcordanceCounter,
                        counter_weighted: ConcordanceCounter) -> tuple:
    """ Compare genotype concordance of two variants across samples """

    # vectorize genotypes
    invar_gt = [gt2code(invar.samples[sample]['GT']) for sample in in_samples]
    refvar_gt = [gt2code(refvar.samples[sample]['GT']) for sample in ref_samples]

    invar_gt = np.array(invar_gt, dtype=np.int8)
    refvar_gt = np.array(refvar_gt, dtype=np.int8)

    # find non-missing genotypes in both VCF
    if include_missing:
        idx_cmp = refvar_gt != -1
    else:
        idx_cmp = (invar_gt != -1) & (refvar_gt != -1)

    if np.sum(idx_cmp) == 0:
        return (None, None, None)

    invar_cmp = invar_gt[idx_cmp]
    refvar_cmp = refvar_gt[idx_cmp]

    # concordance by genotype
    con_gt = invar_cmp == refvar_cmp
    gc = np.mean(con_gt)
    # concordance by existence
    con_exist = (invar_cmp > 0) == (refvar_cmp > 0)
    gc_exist = np.mean(con_exist)
    # weighted genotype concordance
    idx_0 = refvar_cmp == 0
    idx_1 = refvar_cmp == 1
    idx_2 = refvar_cmp == 2
    gc_0 = np.mean(con_gt[idx_0]) if np.sum(idx_0) > 0 else 1
    gc_1 = np.mean(con_gt[idx_1]) if np.sum(idx_1) > 0 else 1
    gc_2 = np.mean(con_gt[idx_2]) if np.sum(idx_2) > 0 else 1
    gc_weighted = (gc_0 + gc_1 + gc_2) / 3

    counter_gt.add(gc)
    counter_exist.add(gc_exist)
    counter_weighted.add(gc_weighted)

    return float(gc), float(gc_exist), float(gc_weighted)


def find_variants(vcf:pysam.VariantFile, contig: str, pos: int, pass_only: bool) -> list:
    """ Find all variants at a given position """

    var_lst = []

    # start = pos - 1
    for variant in vcf.fetch(contig=contig, start=pos-1, stop=pos):
        if variant.pos == pos:
            if pass_only and ('PASS' not in variant.filter):
                continue

            assert len(variant.alleles) == 2
            var_lst.append(variant)
    
    return var_lst


def process_pos(prev_var: pysam.VariantRecord, 
                invar_lst: list, 
                refvcf: pysam.VariantFile, 
                outvcf: pysam.VariantFile, 
                compare: str,
                in_samples: list,
                ref_samples: list,
                include_missing: bool,
                pass_only: bool,
                counter_gt: ConcordanceCounter, 
                counter_exist: ConcordanceCounter,
                counter_weighted: ConcordanceCounter) -> None:
    """ Process all variants at the same position """

    # only 1 variant at the previous position:
    if len(invar_lst) == 0:
        invar_lst.append(prev_var)

    refvar_lst = find_variants(refvcf, prev_var.contig, prev_var.pos, pass_only)
    # key -> variant to index variants
    refvar_map = {}

    for refvar in refvar_lst:
        refvar_key = refvar.alleles if compare == "alleles" else refvar.id
        refvar_map[refvar_key] = refvar
    
    for invar in invar_lst:
        invar_key = invar.alleles if compare == "alleles" else invar.id
        if invar_key in refvar_map:
            refvar = refvar_map[invar_key]
            con_gt, con_exist, con_weighted = variant_concordance(invar, 
                                                                  refvar, 
                                                                  in_samples, 
                                                                  ref_samples,
                                                                  include_missing,
                                                                  counter_gt, 
                                                                  counter_exist, 
                                                                  counter_weighted)
        else:
            con_gt = None
            con_exist = None
            con_weighted = None
        
        outvar = invar.copy()
        outvar.info['CON_GT'] = con_gt
        outvar.info['CON_EXIST'] = con_exist
        outvar.info['CON_WEIGHTED'] = con_weighted
        outvcf.write(outvar)


def add_header(header: pysam.VariantHeader) -> pysam.VariantHeader:
    """Add header to output VCF"""
    header.info.add('CON_GT', '1', 'Float', 'Percentage of samples with concordant genotypes')
    header.info.add('CON_EXIST', '1', 'Float', 'Percentage of samples with concordant existence of variants')
    header.info.add('CON_WEIGHTED', '1', 'Float', 'Weighted genotype concordance')
    return header


def read_map(mapfile: str) -> tuple:
    """Read sample map file"""

    in_samples = []
    ref_samples = []
    with open(mapfile, 'r') as f:
        for line in f:
            sample1, sample2 = line.strip().split('\t')
            in_samples.append(sample1)
            ref_samples.append(sample2)

    return in_samples, ref_samples


def concordanceVCF_main(cmdargs):
    """Main function to evaluate genotype concordance of two VCF files"""
    args = parse_cmdargs(parser, cmdargs)

    logger = logging.getLogger(__name__)
    logger.info('Calculate genotype concordance between two VCFs')

    # read input
    invcf = read_vcf(args.invcf)
    refvcf = read_vcf(args.refvcf)

    # generate output
    new_header = invcf.header
    new_header = add_header(new_header)
    outvcf = pysam.VariantFile(args.outvcf, 'w', header=new_header)

    # setup counters
    counter_gt = ConcordanceCounter()
    counter_exist = ConcordanceCounter()
    counter_weighted = ConcordanceCounter()

    # get sample map
    if args.map is not None:
        in_samples, ref_samples = read_map(args.map)
    else:
        sample_invcf = list(invcf.header.samples)
        sample_refvcf = list(refvcf.header.samples)
        in_samples = [x for x in sample_invcf if x in sample_refvcf]
        ref_samples = in_samples
    
    # check key
    if args.compare not in ['alleles', 'id']:
        raise ValueError(f'--compare must be either "alleles" or "id"')
    
    logger.info(f'Find {len(in_samples)} samples in both VCF files to compare.')

    # loop through all variants of invcf
    verbosity = 10000
    prev_var = None
    invar_working = []
    if args.region is not None:
        invcf_iter = invcf.fetch(region=args.region)
    else:
        invcf_iter = invcf.fetch()
    for i, cur_var in enumerate(invcf_iter):
        if i % verbosity == 0:
            logger.info(f"Processing variant at {cur_var.chrom}:{cur_var.pos} ({cur_var.id})...")
        
        # make sure input is biallelic
        assert len(cur_var.alleles) == 2

        # get the first variant
        if prev_var is None:
            prev_var = cur_var.copy()
            continue

        # add variants at the same position
        if cur_var.pos == prev_var.pos:
            if len(invar_working) == 0:
                invar_working.append(prev_var)
            invar_working.append(cur_var)
        # seen all variants at the pos, process
        elif cur_var.pos > prev_var.pos:
            process_pos(prev_var, invar_working, refvcf, outvcf, args.compare, in_samples, ref_samples, args.include_missing, args.pass_only, counter_gt, counter_exist, counter_weighted)
            invar_working = []
        # switch to a new chromosome
        elif cur_var.chrom != prev_var.chrom:
            process_pos(prev_var, invar_working, refvcf, outvcf, args.compare, in_samples, ref_samples, args.include_missing, args.pass_only, counter_gt, counter_exist, counter_weighted)
            invar_working = []
        else:
            raise ValueError("Input VCF is not sorted by position")
        
        prev_var = cur_var.copy()
    
    # process the last position
    process_pos(prev_var, invar_working, refvcf, outvcf, args.compare, in_samples, ref_samples, args.include_missing, args.pass_only, counter_gt, counter_exist, counter_weighted)
   
    # summarise concordance results
    logger.info(f'Compared {counter_gt.total} variants')
    logger.info(f'Average genotype concordance: {counter_gt.rate()}')
    logger.info(f'Average existence concordance: {counter_exist.rate()}')
    logger.info(f'Average weighted genotype concordance: {counter_weighted.rate()}')
    
    # close files
    invcf.close()
    refvcf.close()
    outvcf.close()

    logger.info('Done')

    return 0