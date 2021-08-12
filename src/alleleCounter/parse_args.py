# modules
import sys
import argparse
from distutils.util import strtobool

# argparse
def parse_args(program_version, arguments=sys.argv[1:]):
    # main_arguments
    parser = argparse.ArgumentParser(
        add_help=True,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="alleleCounter return A, T, G, C and deletion counts for each reference position",
    )
    parser.add_argument(
        "-i",
        "--reads",
        type=str,
        required=True,
        help="minimap2 (parameters: -ax asm20 -cs) aligned SAM/BAM files",
    )
    parser.add_argument(
        "--min_bq",
        type=int,
        default=1,
        required=False,
        help="base quality score threshold (default = 1)",
    )
    parser.add_argument(
        "--min_mapq",
        type=int,
        default=60,
        required=False,
        help="minimum mapping quality score (default = 60)",
    )
    parser.add_argument(
        "--loci",
        type=str,
        required=False,
        help="target loci (chr\tpos\tend)",
    )
    parser.add_argument(
        "--region",
        type=str,
        required=False,
        help="target chromosome",
    )
    parser.add_argument(
        "--region_list",
        type=str,
        required=False,
        help="target chromsome list (one per line, default = 1..22)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        required=False,
        help="maximum number of threads to be used (default = 1)",
    )
    parser.add_argument(
        "-o",
        "--txt",
        type=str,
        required=True,
        help="file to return alleleCounts",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="\n%(prog)s {version}\n".format(version=program_version),
    )
    # no arguments
    if len(arguments) == 0:
        parser.print_help()
        parser.exit()
    else:
        return parser.parse_args(arguments)
