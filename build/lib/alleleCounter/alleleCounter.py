import os
import sys
import time
import pysam
import alleleCounter.cslib 
import alleleCounter.bamlib
## import alleleCounter.caller
import multiprocessing as mp
from typing import Dict, List, Tuple


class BAM:
    def __init__(self, line):
        # target
        self.tname = line.reference_name
        self.tstart = line.reference_start
        self.tend = line.reference_end
        # query
        self.qname = line.query_name
        self.qstart = line.query_alignment_start
        self.qend = line.query_alignment_end
        self.qseq = line.query_sequence
        self.qlen = len(line.query_sequence)
        self.bq_int_lst = line.query_qualities
        self.mapq = line.mapping_quality
        self.cs_state = line.has_tag("cs")
        self.cs_tag = line.get_tag("cs") if line.has_tag("cs") else "."
        self.alignment_type = line.get_tag("tp") if line.has_tag("tp") else "."


genome_allelecount_lst = []
def collect_allelecounts(result):
    genome_allelecount_lst.extend(result)


def load_loci(loci: str) -> Dict[str, Tuple[str, int, int]]:
    loci_hsh = {}
    for line in open(loci):
        chr, start, end = line.strip().split()
        if chr not in loci_hsh:
            loci_hsh[chr] = []
        loci_hsh[chr].append((chr, int(start), int(end)))
    return loci_hsh


def load_region(
    region: str, 
    region_file: str
) -> List[str]:
    if region is not None and region_file is None:
        region_lst = [region]
    elif region is None and region_file is not None:
        region_lst = [line.strip() for line in open(region_file)]
    elif region is not None and region_file is not None:
        print("--region and --region_file parameters are both provided") 
        sys.exit(0)
    elif region is None and region_file is None:
        print("--region or --region_file parameters are required if --loci parameter is not provided") 
        sys.exit(0)
    return region_lst


def bam2alleleCounts(
    alignments, 
    loci_lst: List[Tuple[str, int, str]],
    min_bq: int,
    min_mapq: int,
) -> List[Tuple[str, int, int, int, int, int, int, int, float]]:

    counter = 0
    pos2counts = {}
    for chrom, loci_start, loci_end in loci_lst:
        for pos in range(loci_start, loci_end):
            pos2counts[(chrom, pos)] = {}
            pos2counts[(chrom, pos)]["A"] = 0
            pos2counts[(chrom, pos)]["T"] = 0
            pos2counts[(chrom, pos)]["G"] = 0
            pos2counts[(chrom, pos)]["C"] = 0
            pos2counts[(chrom, pos)]["-"] = 0
            pos2counts[(chrom, pos)]["bq"] = []
        counter += 1 
        if counter == 2:
            break
    
    counter = 0
    if min_bq == 1:
        for chrom, loci_start, loci_end in loci_lst:
            for line in alignments.fetch(chrom, loci_start, loci_end):
                read = BAM(line)
                if read.cs_state and read.mapq >= min_mapq and read.alignment_type == "P":
                    cstuple_lst = alleleCounter.cslib.cs2tuple(read.cs_tag, read.qstart, read.qseq)
                    pos2alleles = alleleCounter.cslib.cs2alleles(cstuple_lst, read.tstart, read.qstart, read.bq_int_lst)
                    for pos in range(loci_start, loci_end):
                        if pos not in pos2alleles: continue
                        allele_state, _ref, alt, bq = pos2alleles[pos]
                        if allele_state == 1 or allele_state == 2:
                            if (chrom, pos) in pos2counts:
                                pos2counts[(chrom, pos)][alt] += 1
                                pos2counts[(chrom, pos)]["bq"].append(bq)
                        elif allele_state == 4:
                            if (chrom, pos) in pos2counts:
                                pos2counts[(chrom, pos)]["-"] += 1
                                pos2counts[(chrom, pos)]["bq"].append(bq)
            counter += 1 
            if counter == 2:
                break
    return pos2counts

def count(
    bam_file: str, 
    loci: str,
    region: str, 
    region_file: str, 
    min_bq: int, 
    min_mapq: int, 
    threads: int, 
    version: str, 
    out_file: str
) -> None: 

    state = 1
    if bam_file is None or not os.path.exists(os.path.abspath(bam_file)):
        state = 0
        print("BAM file is missing")
    if state == 0:
        sys.exit()

    start = time.time()/60
    print("loading loci list")
    if loci is not None and (region is None or region_file is None):
        chrom_loci_hsh = load_loci(loci)
    elif loci is None and (region is not None or region_file is not None):
        target_lst = load_region(region, region_file)
        chrom_loci_hsh = {target: [] for target in target_lst}
        tname2tsize_map = alleleCounter.bamlib.get_tname2tsize(bam_file)
        for target in target_lst:
            tsize = tname2tsize_map[target]
            chrom_start_lst = range(0, tsize, 100000)
            for chrom_start in chrom_start_lst:
                chrom_end = chrom_start + 100000
                chrom_loci_hsh[target].append((target, chrom_start, chrom_end))

    print("alleleCounter started counting alleles with {} threads".format(threads))
    if threads == 1:
        alignments = pysam.AlignmentFile(bam_file, "rb", header=True)
        for chrom, loci in chrom_loci_hsh.items():
            result = bam2alleleCounts(alignments, chrom_loci_hsh[chrom], min_bq, min_mapq)
            genome_allelecount_lst.append(result)
        alignments.close()
    elif threads > 1:
        # p = mp.Pool(threads)
        # p.join()
        # p.close()
        pass
    print("alleleCounter finished counting alleles {} threads".format(threads))

    print("returning alleleCounts")
    o = open(out_file, "w")
    o.write("{}\n".format("\t".join(["CHROM", "POS", "A", "T", "G", "C", "DEL", "TOTAL", "BQ"])))
    for chrom_allelecount in genome_allelecount_lst:
        for coord in chrom_allelecount:
            chrom, pos = coord
            A = chrom_allelecount[coord]["A"] 
            T = chrom_allelecount[coord]["T"] 
            G = chrom_allelecount[coord]["G"] 
            C = chrom_allelecount[coord]["C"] 
            del_count = chrom_allelecount[coord]["-"] 
            total_count = sum([A, T, G, C, del_count])
            if total_count == 0:
                if del_count != 0: 
                    o.write("{}\t{}\t0\t0\t0\t0\t{}\t0\t0\n".format(chrom, pos, del_count))
            else:
                mean_bq = sum(chrom_allelecount[coord]["bq"])/ len(chrom_allelecount[coord]["bq"])
                o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}\n".format(chrom, pos, A, T, G, C, del_count, total_count, mean_bq))
    o.close()
    print("finished returning alleleCounts")

    end = time.time()/60
    duration = end - start
    print("alleleCounter took {} minutes".format(duration))
