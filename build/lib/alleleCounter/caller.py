# modules
import os
import sys
import time
import json
import pysam
import logging
import itertools
import himut.cslib 
import himut.bamlib
import himut.haplib
import himut.mutlib 
import himut.filter
import himut.vcflib
from typing import Set, Dict, List, Tuple


class BAM:
    def __init__(self, line):
        # target
        self.tname = line.reference_name
        self.tstart = int(line.reference_start)
        self.tend = int(line.reference_end)
        # query
        self.qname = line.query_name
        self.qstart = line.query_alignment_start
        self.qend = line.query_alignment_end
        self.zmw = self.qname.split("/")[1]
        self.qseq = line.query_sequence
        self.qlen = line.query_length
        self.bq_int_lst = line.query_qualities
        hq_base_count = 0
        for i in self.bq_int_lst:
            if i == 93:
                hq_base_count +=1
        self.hq_base_fraction = hq_base_count/float(self.qlen)
        self.mapq = line.mapping_quality
        self.cs_tag = line.get_tag("cs") if line.has_tag("cs") else "."
        self.alignment_type = line.get_tag("tp") if line.has_tag("tp") else "."

MAJOR_HBIT = "0" 
MAJOR_HBIT_LST = ["00", "10"]


def init_state():
    return 0


def load_region(region, regions_file):
    if region is None and regions_file is None:
        chrom_lst = [str(_) for _ in range(1,23)]
        print("using default target regions: 1..22")
    elif region is None and regions_file is not None:
        chrom_lst = [line.strip() for line in open(regions_file).readlines()]
    elif region is not None and regions_file is None:
        chrom_lst = [region]
    elif region is not None and regions_file is not None:
        chrom_lst = [line.strip() for line in open(regions_file).readlines()]
    return chrom_lst 


def get_allele_counts(
    tsub: Tuple[str, int, str, str],
    tsub_read_lst: List[str],
    read_pos_alleles_map: Dict[str, Dict[int, List[str]]]
) -> Tuple[List[str], List[str], int, int, int]:

    bq_lst = []
    ref_count = 0
    alt_count = 0
    del_count = 0
    tsub_target_read_lst = []
    chrom, pos, ref, alt = tsub
    tref = (ref, ref) 
    talt = (ref, alt)
    for qread in tsub_read_lst:
        pos_allele_map = read_pos_alleles_map[qread]
        qallele = pos_allele_map[pos] if pos in pos_allele_map else ("-", "-", None) 
        qbq = qallele[2]
        qmut = qallele[0:2]
        if qmut == tref: 
            ref_count += 1
        elif qmut == talt:
            alt_count += 1
            bq_lst.append(qbq)
            tsub_target_read_lst.append(qread)
        else: 
            del_count += 1
    mean_bq = sum(bq_lst)/len(bq_lst)
    total_count = ref_count + alt_count + del_count
    return tsub_target_read_lst, ref_count, alt_count, total_count, mean_bq


def get_chrom_tsub_allelecounts(
    chrom: str,
    bam_file: str,
    tsub_lst: List[Tuple[str, int, str, str]],
) -> List[Tuple[str, int, str, str, int, int, int]]:

    tsub_count = len(tsub_lst) - 1
    index_lst = list(range(0, tsub_count, 20))
    if tsub_count not in index_lst: 
        index_lst.append(tsub_count)
    index_len = len(index_lst)

    tsub_allele_count_lst = []
    for i, index_start in enumerate(index_lst):
        if i + 1 == index_len: break
        index_end = index_lst[i+1] 
        read_tpos_alleles_map = himut.bamlib.get_region_read_alleles(bam_file, tsub_lst, chrom, index_start, index_end)
        for j in range(index_start, index_end):
            tsub = tsub_lst[j]
            _chrom, pos, ref, alt = tsub
            tsub_read_lst = himut.bamlib.get_tsub_intersecting_reads(bam_file, tsub)
            _tsub_target_read_lst, ref_count, alt_count, total_count, mean_bq = get_allele_counts(tsub, tsub_read_lst, read_tpos_alleles_map)
            vaf = alt_count/float(total_count)
            tsub_allele_count_lst.append((chrom, pos, ref, alt, ref_count, alt_count, total_count, vaf, mean_bq))
    return tsub_allele_count_lst

def get_major_minor_hbit(state, hbit_count_map):

    major_minor_hbit_candidate_lst = []
    if state == 1:
        hbit = "0"
        hbit_count = hbit_count_map["0"]
        hbit_complement = himut.haplib.HBIT_COMPLEMENT_MAP[hbit]
        hbit_complement_count = hbit_count_map[hbit_complement]
        if hbit_count != 0 and hbit_complement_count != 0:
            if hbit_count > hbit_complement_count:
                major_minor_hbit_candidate_lst.append((hbit, hbit_complement))
            else:
                major_minor_hbit_candidate_lst.append((hbit_complement, hbit))
    elif state == 2:
        for hbit in MAJOR_HBIT_LST:
            hbit_count = hbit_count_map[hbit]
            hbit_complement = himut.haplib.HBIT_COMPLEMENT_MAP[hbit]
            hbit_complement_count = hbit_count_map[hbit_complement]
            if hbit_count != 0 and hbit_complement_count != 0:
                if hbit_count > hbit_complement_count:
                    major_minor_hbit_candidate_lst.append((hbit, hbit_complement))
                else:
                    major_minor_hbit_candidate_lst.append((hbit_complement, hbit))
    hap_count = len(major_minor_hbit_candidate_lst)
    return hap_count, major_minor_hbit_candidate_lst


def get_ref_alt_hbit(state, tsub_target_read_lst, read_hbit_map, hbit_count_map):

    ref_alt_hbit_candidate_lst = []
    alt_hbit_lst = [read_hbit_map[target_read] for target_read in tsub_target_read_lst if target_read in read_hbit_map]
    if state == 1:
        alt_hbit_lst = [_ for _ in alt_hbit_lst if _ == MAJOR_HBIT]
        for alt_hbit in alt_hbit_lst:
            ref_hbit = himut.haplib.HBIT_COMPLEMENT_MAP[alt_hbit] 
            ref_hbit_count = hbit_count_map[ref_hbit]
            alt_hbit_count = hbit_count_map[alt_hbit]
            if ref_hbit_count != 0 and alt_hbit_count != 0:
                ref_alt_hbit_candidate_lst.append((ref_hbit, alt_hbit, ref_hbit_count, alt_hbit_count))
    elif state == 2:
        alt_hbit_lst = [_ for _ in alt_hbit_lst if _ in MAJOR_HBIT_LST]
        for alt_hbit in alt_hbit_lst:
            ref_hbit = himut.haplib.HBIT_COMPLEMENT_MAP[alt_hbit] 
            ref_hbit_count = hbit_count_map[ref_hbit]
            alt_hbit_count = hbit_count_map[alt_hbit]
            if ref_hbit_count != 0 and alt_hbit_count != 0:
                ref_alt_hbit_candidate_lst.append((ref_hbit, alt_hbit, ref_hbit_count, alt_hbit_count))
    ref_alt_hbit_candidate_lst = list(set(ref_alt_hbit_candidate_lst))
    hap_count = len(ref_alt_hbit_candidate_lst)
    return hap_count, ref_alt_hbit_candidate_lst


def get_alt_hbit_counts(alt_hbit, alt_hbit_count, tsub_target_read_lst, read_hbit_map):
    target_read_alt_hbit_count = 0
    for target_read in tsub_target_read_lst:
        if target_read in read_hbit_map and read_hbit_map[target_read] == alt_hbit:
            target_read_alt_hbit_count += 1
    query_read_alt_hbit_count = alt_hbit_count - target_read_alt_hbit_count
    return target_read_alt_hbit_count, query_read_alt_hbit_count


def phase_subs(
    chrom: str,
    tsub_lst: List[Tuple[str, int, str, str]],
    sample_hetsnps: Set[Tuple[str, int, str, str]],
    bam_file: str,
    min_ref_count: int,
    min_alt_count: int,
    md_threshold: float,
    min_hap_count: int,
    min_hap_proportion: float
) -> List[Tuple[str, int, str, str, int, int, int, float, float]]:

    tsub_count = len(tsub_lst) - 1
    index_lst = list(range(0, tsub_count, 20))
    if tsub_count not in index_lst: 
        index_lst.append(tsub_count)
    index_len = len(index_lst)

    phased_tsub_lst = []
    unphaseable_tsub_count = 0
    hetsnp_count = len(sample_hetsnps)
    hetsnp_positions = [hetsnp[1] for hetsnp in sample_hetsnps]
    for i, index_start in enumerate(index_lst):
        if i + 1 == index_len: 
            break
        index_end = index_lst[i+1] 
        read_tpos_alleles_map = himut.bamlib.get_region_read_alleles(bam_file, tsub_lst, chrom, index_start, index_end)

        for j in range(index_start, index_end):  
            tsub = tsub_lst[j]
            tsub_read_lst = himut.bamlib.get_tsub_intersecting_reads(bam_file, tsub)
            tsub_target_read_lst, ref_count, alt_count, total_count, mean_bq = get_allele_counts(tsub, tsub_read_lst, read_tpos_alleles_map)
            if ref_count < min_ref_count or total_count > md_threshold or alt_count < min_alt_count: # ab and md filter
                continue

            hetsnp_state, upstream_hetsnp, downstream_hetsnp = himut.haplib.get_tsub_adjacet_hetsnps(tsub[1], tsub_target_read_lst, read_tpos_alleles_map, sample_hetsnps, hetsnp_count, hetsnp_positions)
            if hetsnp_state == 0: # no adjacent hetsnps
                unphaseable_tsub_count += 1
                continue
            
            read_hbit_map, hbit_count_map = himut.haplib.get_read_hbits(tsub_read_lst, hetsnp_state, upstream_hetsnp, downstream_hetsnp, read_tpos_alleles_map)
            hap_count, major_minor_hbit_candidate_lst = get_major_minor_hbit(hetsnp_state, hbit_count_map)
            target_read_hap_count, ref_alt_hbit_candidate_lst = get_ref_alt_hbit(hetsnp_state, tsub_target_read_lst, read_hbit_map, hbit_count_map)
            if hap_count == 1 and target_read_hap_count == 1: # one possible haplotype # target reads support a single haplotype 
                major_hbit, minor_hbit = major_minor_hbit_candidate_lst[0]
                ref_hbit, alt_hbit, ref_hbit_count, alt_hbit_count = ref_alt_hbit_candidate_lst[0]

                if ((ref_hbit == major_hbit and alt_hbit == minor_hbit) or (ref_hbit == minor_hbit and alt_hbit == major_hbit)): # haplotypes of the read with the substitution is supported
                    hbit_total_count = sum(hbit_count_map.values())
                    ref_alt_hbit_proportion = (ref_hbit_count + alt_hbit_count)/float(hbit_total_count)
                    target_read_alt_hbit_count, query_read_alt_hbit_count = get_alt_hbit_counts(alt_hbit, alt_hbit_count, tsub_target_read_lst, read_hbit_map)
                    if ref_hbit_count >= min_hap_count and query_read_alt_hbit_count >= min_hap_count and ref_alt_hbit_proportion >= min_hap_proportion: # PASS
                        vaf = alt_count/float(total_count)
                        _chrom, pos, ref, alt = tsub_lst[j]
                        phased_tsub_lst.append((chrom, pos, ref, alt, ref_count, alt_count, total_count, vaf, mean_bq))
                    else: # does not meet minimum haplotype filtering conditions
                        continue
                else: # read with the substitution has a deletion in the upstream or the downstream hetsnp position
                    continue
            else: # no haplotype # multiple haplotype # multiple haplotype candidates for target reads # bad regions
                # print("bad")
                # print("tsub:{}, upstream_hetsnp:{}, downstream_hetsnp: ({})".format(tsub, upstream_hetsnp, downstream_hetsnp))
                # print(ref_alt_hbit_candidate_lst, major_minor_hbit_candidate_lst)
                # print(json.dumps(hbit_count_map,indent=4))
                continue
    return phased_tsub_lst
