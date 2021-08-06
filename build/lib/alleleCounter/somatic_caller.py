import os
import sys
import time
import pysam
import natsort
import logging
import himut.cslib 
import himut.mutlib 
import himut.bamlib
import himut.filter
import himut.vcflib
import himut.caller
import multiprocessing as mp
from typing import Set, List, Tuple


genome_tsub_allelecount_lst = []
def collect_substitutions(result):
    genome_tsub_allelecount_lst.extend(result)


def get_chrom_phased_somatic_substitutions(
    chrom: str,
    bam_file: str, 
    min_mapq: int,
    min_hq_base_fraction: float,
    min_sequence_identity: float,
    pon_subs: Set[List[Tuple[str, int, str, str]]],
    common_snps: Set[List[Tuple[str, int, str, str]]],
    sample_snps: Set[List[Tuple[str, int, str, str]]],
    sample_hetsnps: Set[List[Tuple[str, int, str, str]]],
    min_bq: int,
    mismatch_window: int,
    min_mismatch_count: int,
    min_ref_count: int, 
    min_alt_count: int,
    md_threshold: float,
    min_hap_count: int,
    min_hap_proportion: float,
    verbose: bool
) -> List[Tuple[str, int, str, str, int, int, int, float, float]]:

    counter = 0
    chrom_tsub_lst = []
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for line in alignments.fetch(chrom):
        read = himut.caller.BAM(line)
        if read.mapq < min_mapq or read.hq_base_fraction < min_hq_base_fraction or read.alignment_type != "P": # select primary alignments 
            continue

        cstuple_lst = himut.cslib.cs2tuple(read.cs_tag, read.qstart, read.qseq) 
        sequence_identity = himut.cslib.blast_sequence_identity(cstuple_lst, read.qstart, read.qend, read.qlen)
        if sequence_identity < min_sequence_identity:  # read filter
            continue

        mut_state, bq_state, tsub_lst, qsub_lst, qbq_lst, tindel_lst, tmbs_lst, tcomplex_lst = himut.cslib.cs2mut(
            cstuple_lst,
            read.tname,
            read.tstart,
            read.zmw,
            read.qstart,
            min_bq,
            read.bq_int_lst,
        )
        if not bq_state or not mut_state[0]: # reads with substitutions
            continue

        contamination_state = himut.filter.get_sample_contamination_state(tsub_lst, common_snps, sample_snps)
        if contamination_state:
            continue

        tsub_state, tsub_candidate_lst, qsub_candidate_lst, qbq_candidate_lst = himut.filter.get_somatic_substitution_candidates(tsub_lst, qsub_lst, min_bq, qbq_lst, pon_subs, sample_snps)
        if not tsub_state: # no substitution
            continue

        mismatch_state, tsub_candidate_lst = himut.filter.get_mismatch_state(tsub_candidate_lst, qsub_candidate_lst, qbq_candidate_lst, tindel_lst, tmbs_lst, tcomplex_lst, read.qlen, min_bq, mismatch_window, min_mismatch_count)
        if mismatch_state: # substitutions with adjacent mismatches 
            chrom_tsub_lst.extend(tsub_candidate_lst)
            counter += 1
            # if counter > 30:
            #     break
    alignments.close()

    chrom_tsub_lst = natsort.natsorted(list(set(chrom_tsub_lst)))
    if verbose:
        chrom_tsub_allelecount_lst = himut.caller.get_chrom_tsub_allelecounts(chrom, bam_file, chrom_tsub_lst) 
        chrom_phased_tsub_allelecount_lst = himut.caller.phase_subs(chrom, chrom_tsub_lst, sample_hetsnps, bam_file, min_ref_count, min_alt_count, md_threshold, min_hap_count, min_hap_proportion)
        return chrom_tsub_allelecount_lst, chrom_phased_tsub_allelecount_lst
    else:
        chrom_phased_tsub_allelecount_lst = himut.caller.phase_subs(chrom, chrom_tsub_lst, sample_hetsnps, bam_file, min_ref_count, min_alt_count, md_threshold, min_hap_count, min_hap_proportion)
        return chrom_phased_tsub_allelecount_lst


def call_somatic_substitutions(
    bam_file: str, 
    region: str, 
    regions_file: str, 
    ploidy: str, 
    sample_is_reference: bool, 
    min_mapq: int, 
    min_sequence_identity: float,
    min_hq_base_fraction: float,
    common_snps: str,
    sample_snps: str,
    panel_of_normals: str, 
    min_bq: int, 
    mismatch_window: int,
    min_mismatch_count: int,
    min_ref_count: int,
    min_alt_count: int,
    min_hap_count: int,
    min_hap_proportion: float,
    # sigprofile: str,
    verbose: bool,
    threads: int, 
    version: str, 
    vcf_file: str
) -> None: # delete hbool later

    state = 1
    if not vcf_file.endswith(".vcf"):
        state = 0
        print("vcf file must have the .vcf suffix")
    if bam_file is None or not os.path.exists(os.path.abspath(bam_file)):
        state = 0
        print("BAM file is missing")
    if common_snps is None or not os.path.exists(os.path.abspath(common_snps)):
        state = 0
        print("common SNPs mut file is missing")
    if sample_snps is None or not os.path.exists(os.path.abspath(sample_snps)):
        state = 0
        print("germline SNPs VCF file is missing")
    if panel_of_normals is None or not os.path.exists(os.path.abspath(panel_of_normals)):
        state = 0
        print("germline SNPs VCF file is missing")
    if state == 0:
        sys.exit()


    start = time.time()/60
    print("calling and phasing substitutions with {} threads".format(threads))
    p = mp.Pool(threads)
    chrom_lst = himut.caller.load_region(region, regions_file, )
    tname_lst, tname2tsize = himut.bamlib.get_tname2tsize(bam_file)
    chrom2pon_subs = himut.vcflib.load_pon_subs(panel_of_normals, tname_lst, chrom_lst)
    chrom2common_snps = himut.mutlib.load_mut_file(common_snps, tname_lst, chrom_lst)
    chrom2md_threshold = himut.bamlib.get_target_md_threshold(bam_file, chrom_lst, tname2tsize, threads)
    chrom2sample_snps = himut.vcflib.load_sample_snps(sample_snps, tname_lst, chrom_lst, ploidy, sample_is_reference)
    chrom2sample_hetsnps = himut.vcflib.load_sample_hetsnps(sample_snps, tname_lst, chrom_lst, ploidy, sample_is_reference)
    get_chrom_phased_substitution_arg_lst = [(
        chrom, 
        bam_file, 
        min_mapq, 
        min_hq_base_fraction, 
        min_sequence_identity,
        chrom2pon_subs[chrom], 
        chrom2common_snps[chrom], 
        chrom2sample_snps[chrom],
        chrom2sample_hetsnps[chrom],
        min_bq,
        mismatch_window,
        min_mismatch_count,
        min_ref_count, 
        min_alt_count,
        chrom2md_threshold[chrom],
        min_hap_count, 
        min_hap_proportion,
        verbose
    ) for chrom in chrom_lst]
    p.starmap_async(get_chrom_phased_somatic_substitutions, get_chrom_phased_substitution_arg_lst, callback = collect_substitutions)
    chrom2sample_hetsnps.clear()
    chrom2md_threshold.clear()
    chrom2sample_snps.clear()
    chrom2common_snps.clear()
    chrom2pon_subs.clear()
    p.close()
    p.join()
    print("finished phasing and calling substitutions with {} threads".format(threads))

    print("returning substitutions")
    if verbose:
        vcf_header = himut.bamlib.get_vcf_header(bam_file, version)
        phased_vcf_file = open(vcf_file.replace(".vcf", ".phased.vcf"), "w")
        unphased_vcf_file = open(vcf_file.replace(".vcf", ".unphased.vcf"), "w") 
        phased_vcf_file.write("{}\n".format(vcf_header))
        unphased_vcf_file.write("{}\n".format(vcf_header))
        for chrom_tsub_tuple_lst in genome_tsub_allelecount_lst:
            for tsub in chrom_tsub_tuple_lst[0]:
                chrom, pos, ref, alt, ref_count, alt_count, total_count, vaf, mean_bq = tsub
                unphased_vcf_file.write("{}\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT:BQ:DP:AD:VAF\t./.:{:.0f}:{}:{},{}:{:.2f}\n".format(chrom, pos, ref, alt, mean_bq, total_count, ref_count, alt_count, vaf))
            for tsub in chrom_tsub_tuple_lst[1]:
                chrom, pos, ref, alt, ref_count, alt_count, total_count, vaf, mean_bq = tsub
                phased_vcf_file.write("{}\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT:BQ:DP:AD:VAF\t./.:{:.0f}:{}:{},{}:{:.2f}\n".format(chrom, pos, ref, alt, mean_bq, total_count, ref_count, alt_count, vaf))
        phased_vcf_file.close()
        unphased_vcf_file.close()

    else:
        vcf_header = himut.bamlib.get_vcf_header(bam_file, version)
        phased_vcf_file = open(vcf_file, "w")
        phased_vcf_file.write("{}\n".format(vcf_header))
        for chrom_tsub_tuple_lst in genome_tsub_allelecount_lst:
            for tsub in chrom_tsub_tuple_lst:
                chrom, pos, ref, alt, ref_count, alt_count, total_count, vaf, mean_bq = tsub
                phased_vcf_file.write("{}\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT:BQ:DP:AD:VAF\t./.:{:.0f}:{}:{},{}:{:.2f}\n".format(chrom, pos, ref, alt, mean_bq, total_count, ref_count, alt_count, vaf))
        phased_vcf_file.close()
    print("finished returning substitutions")

    end = time.time()/60
    duration = end - start
    print("himut somatic mutation detection took {} minutes".format(duration))