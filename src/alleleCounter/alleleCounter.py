import os
import sys
import time
import pysam
import natsort
import alleleCounter.cslib
import alleleCounter.bamlib
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


def load_region(region: str, region_file: str) -> List[str]:
    if region is not None and region_file is None:
        region_lst = [region]
    elif region is None and region_file is not None:
        region_lst = [line.strip() for line in open(region_file)]
    elif region is not None and region_file is not None:
        print("--region and --region_file parameters are both provided")
        sys.exit(0)
    elif region is None and region_file is None:
        print(
            "--region or --region_file parameters are required if --loci parameter is not provided"
        )
        sys.exit(0)
    return region_lst


def bam2alleleCounts(
    bam_file,
    loci_lst: List[Tuple[str, int, str]],
    min_bq: int,
    min_mapq: int,
) -> List[Tuple[str, int, int, int, int, int, int, int, float]]:

    pos2counts = {}
    chrom = loci_lst[0][0]
    o = open("{}.alleleCounts".format(chrom), "w") 
    o.write(
        "{}\n".format(
            "\t".join(["CHROM", "POS", "A", "T", "G", "C", "DEL", "TOTAL", "BQ", "MAPQ"])
        )
    )

    alignments = pysam.AlignmentFile(bam_file, "rb")
    if min_bq == 1:
        for chrom, loci_start, loci_end in loci_lst:
            pos2counts = {}
            for pos in range(loci_start, loci_end):
                pos2counts[(chrom, pos)] = {}
                pos2counts[(chrom, pos)]["A"] = 0
                pos2counts[(chrom, pos)]["T"] = 0
                pos2counts[(chrom, pos)]["G"] = 0
                pos2counts[(chrom, pos)]["C"] = 0
                pos2counts[(chrom, pos)]["-"] = 0
                pos2counts[(chrom, pos)]["bq"] = []
                pos2counts[(chrom, pos)]["mapq"] = []

            for line in alignments.fetch(chrom, loci_start, loci_end):
                read = BAM(line)
                if (
                    read.cs_state
                    and read.mapq >= min_mapq
                    and read.alignment_type == "P"
                ):
                    cstuple_lst = alleleCounter.cslib.cs2tuple(
                        read.cs_tag, read.qstart, read.qseq
                    )
                    pos2alleles = alleleCounter.cslib.cs2alleles(
                        cstuple_lst, read.tstart, read.qstart, read.bq_int_lst
                    )
                    for pos in range(loci_start, loci_end):
                        if pos not in pos2alleles:
                            continue
                        allele_state, _ref, alt, bq = pos2alleles[pos]
                        if allele_state == 1 or allele_state == 2:
                            if (chrom, pos) in pos2counts:
                                pos2counts[(chrom, pos)][alt] += 1
                                pos2counts[(chrom, pos)]["bq"].append(bq)
                                pos2counts[(chrom, pos)]["mapq"].append(read.mapq)
                        elif allele_state == 4:
                            if (chrom, pos) in pos2counts:
                                pos2counts[(chrom, pos)]["-"] += 1
                                pos2counts[(chrom, pos)]["bq"].append(bq)
                                pos2counts[(chrom, pos)]["mapq"].append(read.mapq)
            for coord in pos2counts:
                chrom, pos = coord
                A = pos2counts[coord]["A"]
                T = pos2counts[coord]["T"]
                G = pos2counts[coord]["G"]
                C = pos2counts[coord]["C"]
                del_count = pos2counts[coord]["-"]
                total_count = sum([A, T, G, C, del_count])
                if total_count == 0:
                    if del_count != 0:
                        mean_mapq = sum(pos2counts[coord]["mapq"]) / len(
                            pos2counts[coord]["mapq"]
                        )
                        o.write(
                            "{}\t{}\t0\t0\t0\t0\t{}\t0\t0\t{:.2f}\n".format(chrom, pos, del_count, mean_mapq)
                        )
                else:
                    mean_bq = sum(pos2counts[coord]["bq"]) / len(
                        pos2counts[coord]["bq"]
                    )
                    mean_mapq = sum(pos2counts[coord]["mapq"]) / len(
                        pos2counts[coord]["mapq"]
                    )
                    o.write(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}\t{:.2f}\n".format(
                            chrom, pos, A, T, G, C, del_count, total_count, mean_bq, mean_mapq
                        )
                    )
    alignments.close()
    o.close()


def count(
    bam_file: str,
    loci: str,
    region: str,
    region_file: str,
    min_bq: int,
    min_mapq: int,
    threads: int,
    version: str,
    out_file: str,
) -> None:

    state = 1
    if bam_file is None or not os.path.exists(os.path.abspath(bam_file)):
        state = 0
        print("BAM file is missing")
    if state == 0:
        sys.exit()

    start = time.time() / 60
    print("loading loci list")
    if loci is not None and (region is None or region_file is None):
        chrom_loci_hsh = load_loci(loci)
        target_lst = list(chrom_loci_hsh.keys())
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
        for _chrom, chrom_loci in chrom_loci_hsh.items():
            bam2alleleCounts(bam_file, chrom_loci, min_bq, min_mapq)
    elif threads > 1:
        p = mp.Pool(threads)
        chrom_bam2allelecount_arg_lst = [
            (bam_file, chrom_loci, min_bq, min_mapq)
            for _chrom, chrom_loci, in chrom_loci_hsh.items()
        ]
        p.starmap_async(
            bam2alleleCounts,
            chrom_bam2allelecount_arg_lst,
        )
        p.close()
        p.join()
    print("alleleCounter finished counting alleles {} threads".format(threads))

    print("merging alleleCounts")
    o = open(out_file, "w")
    o.write(
        "{}\n".format(
            "\t".join(["CHROM", "POS", "A", "T", "G", "C", "DEL", "TOTAL", "BQ"])
        )
    )
    target_lst = natsort.natsorted(target_lst)
    for target in target_lst:
        for line in open("{}.alleleCounts".format(target)).readlines():
            if line.startswith("CHROM"): continue
            o.write("{}".format(line))
        os.remove("{}.alleleCounts".format(target))
    o.close()  
    print("finished merging alleleCounts")

    end = time.time() / 60
    duration = end - start
    print("alleleCounter took {} minutes".format(duration))
