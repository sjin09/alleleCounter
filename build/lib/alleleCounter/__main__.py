__version__ = "0.0.1"
__author__ = "Sangjin Lee"

# modules
import sys
import multiprocessing as mp
import alleleCounter.alleleCounter 
from alleleCounter.parse_args import parse_args

def main():
    thread_count = mp.cpu_count()
    options = parse_args(program_version=__version__)
    if options.threads > thread_count:
        options.threads = thread_count
    alleleCounter.alleleCounter.count(
        options.reads, # input # bamfile
        options.loci, # loci: chrom\tpos
        options.region, # target contigs/scaffolds/chromosomes
        options.region_list, # target contigs/scaffolds/chromosomes fofn (file of file names)
        options.min_bq, # int: 1 - 93
        options.min_mapq, # int: 0 - 60
        options.threads, # number of threads
        __version__, # true/false print other metrics
        options.txt, # output # himut vcf file
    )
if __name__ == "__main__":
    main()
    sys.exit(0)