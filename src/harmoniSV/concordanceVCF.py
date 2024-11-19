#!/usr/bin/env python3

# Calculate genotype concordance between two VCFs
# Last update: 8-Oct-2024
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
optional_arg.add_argument('--compare', metavar="alleles", type=str, default="alleles", help="Use which information to determine if 2 variants are the same: alleles (default) or id")
optional_arg.add_argument("-h", "--help", action='help', help='show this help message and exit')


class WeightedGC:
    """ Weighted genotype concordance """
    def __init__(self) -> None:
        self.gc_lst = [[], [], []]
    

    def add(self, gc_code: int):
        """ 
        1 or -1 : 0/0 concordance or not
        2 or -2 : 0/1 concordance or not
        3 or -3 : 1/1 concordance or not
        """
        if gc_code is None:
            return
        res_con = gc_code > 0
        self.gc_lst[abs(gc_code) - 1].append(res_con)
    

    def concordance(self) -> float:
        """ Calculate weighted concordance """
        gc_1 = sum(self.gc_lst[0]) / len(self.gc_lst[0]) if len(self.gc_lst[0]) > 0 else 1
        gc_2 = sum(self.gc_lst[1]) / len(self.gc_lst[1]) if len(self.gc_lst[1]) > 0 else 1
        gc_3 = sum(self.gc_lst[2]) / len(self.gc_lst[2]) if len(self.gc_lst[2]) > 0 else 1
        return (gc_1 + gc_2 + gc_3) / 3


class ConcordanceCounter:
    """Accumulator to store concordance statistics"""
    def __init__(self):
        self.total = 0
        self.concordant = 0

    def add(self, n: float):
        """Add concordance statistics"""
        self.total += 1
        self.concordant += n

    def rate(self) -> float:
        """Calculate concordance rate"""
        return self.concordant / self.total



GT_DICT = {(0, 0): 0,
           (0, 1): 1,
           (1, 0): 1,
           (1, 1): 2}

def gt_concordance(gt_in: tuple, gt_ref: tuple) -> tuple:
    """Compare genotype concordance of two genotypes, return concordance by genotype, existence, and weighted concordance code"""
    
    # if any reference genotype is missing, return None and skip
    if None in gt_ref or None in gt_in:
        return (None, None, None)
    ac_ref = GT_DICT[gt_ref]

    # # if any input genotype is missing, return 0 for both
    # if None in gt_in:
    #     return (0, 0, -ac_ref - 1)
    
    ac_in = GT_DICT[gt_in]
    # concordance by genotype
    if ac_in == ac_ref:
        return (1, 1, ac_ref + 1)
    # concordance by existence only
    elif (ac_in > 0) and (ac_ref > 0):
        return (0, 1, -ac_ref - 1)
    # discordance for both
    else:
        return (0, 0, -ac_ref - 1)


def variant_concordance(invar: pysam.VariantRecord, 
                        refvar: pysam.VariantRecord, 
                        sample_map: dict,
                        counter_gt: ConcordanceCounter,
                        counter_exist: ConcordanceCounter,
                        counter_weighted: ConcordanceCounter) -> tuple:
    """Compare genotype concordance of two variants, return concordance by genotype and existence (i.e., 0/1 = 1/1)"""
    con_gt = 0
    con_exist = 0
    n_sample = len(sample_map)
    con_weighted = WeightedGC()

    for sample1, sample2 in sample_map.items():
        gt_in = invar.samples[sample1]['GT']
        gt_ref = refvar.samples[sample2]['GT']
        con_gt_i, con_exist_i, gc_code_i = gt_concordance(gt_in, gt_ref)

        # if input / reference genotype is missing, skip
        if con_gt_i is None:
            n_sample -= 1
            continue

        con_gt += con_gt_i
        con_exist += con_exist_i
        con_weighted.add(gc_code_i)

        counter_gt.add(con_gt_i)
        counter_exist.add(con_exist_i)
    
    # if no sample can be compared, return None
    if n_sample == 0:
        return (None, None)
    else:
        wgc = con_weighted.concordance()
        counter_weighted.add(wgc)
        return (con_gt / n_sample, con_exist / n_sample, wgc)


def find_variants(vcf:pysam.VariantFile, contig: str, pos: int) -> list:
    """ Find all variants at a given position """

    var_lst = []

    # start = pos - 1
    for variant in vcf.fetch(contig=contig, start=pos-1):
        if variant.pos == pos:
            var_lst.append(variant)
        elif variant.pos > pos:
            break
    
    return var_lst


def process_pos(prev_var: pysam.VariantRecord, invar_lst: list, refvcf: pysam.VariantFile, outvcf: pysam.VariantFile, compare: str,
                sample_map: dict, counter_gt: ConcordanceCounter, counter_exist: ConcordanceCounter, counter_weighted: ConcordanceCounter) -> None:
    """ Process all variants at the same position """

    # only 1 variant at the previous position:
    if len(invar_lst) == 0:
        invar_lst.append(prev_var)

    refvar_lst = find_variants(refvcf, prev_var.contig, prev_var.pos)
    # key -> variant to index variants
    refvar_map = {}

    for refvar in refvar_lst:
        refvar_key = refvar.alleles if compare == "alleles" else refvar.id
        refvar_map[refvar_key] = refvar
    
    for invar in invar_lst:
        invar_key = invar.alleles if compare == "alleles" else invar.id
        if invar_key in refvar_map:
            refvar = refvar_map[invar_key]
            con_gt, con_exist, con_weighted = variant_concordance(invar, refvar, sample_map, counter_gt, counter_exist, counter_weighted)
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
    counter_weighted = ConcordanceCounter()

    # get sample map
    if args.map is not None:
        sample_map = read_map(args.map)
    else:
        sample_invcf = list(invcf.header.samples)
        sample_refvcf = list(refvcf.header.samples)
        sample_overlap = [x for x in sample_invcf if x in sample_refvcf]
        sample_map = dict(zip(sample_overlap, sample_overlap))
    
    # check key
    if args.compare not in ['alleles', 'id']:
        raise ValueError(f'--compare must be either "alleles" or "id"')
    
    logger.info(f'Find {len(sample_map)} samples in both VCF files to compare.')

    # loop through all variants of invcf
    verbosity = 2000
    prev_var = None
    invar_working = []
    for i, cur_var in enumerate(invcf.fetch()):
        if i % verbosity == 0:
            logger.info(f"Processing variant at {cur_var.chrom}:{cur_var.pos} ({cur_var.id})...")
        
        # make input is biallelic
        if len(cur_var.alleles) > 2:
            raise ValueError("Input VCF is not biallelic")

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
            process_pos(prev_var, invar_working, refvcf, outvcf, args.compare, sample_map, counter_gt, counter_exist, counter_weighted)
            invar_working = []
        # switch to a new chromosome
        elif cur_var.chrom != prev_var.chrom:
            process_pos(prev_var, invar_working, refvcf, outvcf, args.compare, sample_map, counter_gt, counter_exist, counter_weighted)
            invar_working = []
        else:
            raise ValueError("Input VCF is not sorted by position")
        
        prev_var = cur_var.copy()
    
    # process the last position
    process_pos(prev_var, invar_working, refvcf, outvcf, args.compare, sample_map, counter_gt, counter_exist, counter_weighted)
   
    # summarise concordance results
    logger.info(f'Overall concordance of genotypes: {counter_gt.rate()}')
    logger.info(f'Overall concordance of existence: {counter_exist.rate()}')
    logger.info(f'Overall weighted genotype concordance: {counter_weighted.rate()}')
    
    # close files
    invcf.close()
    refvcf.close()
    outvcf.close()

    logger.info('Done')

    return 0