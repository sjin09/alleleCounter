import sys
import pysam

def get_tname2tsize(bamfile: str):
    tname2tsize_hsh = {}
    bamfile = pysam.AlignmentFile(bamfile, "rb", header=True)
    bam_header_lst = str(bamfile.header).strip().split("\n")
    for h in bam_header_lst:
        if h.startswith("@SQ"):
            _tag, tname, tsize = h.split("\t")
            tname = tname.replace("SN:", "")
            tsize = tsize.replace("LN:", "")
            tsize = int(tsize)
            tname2tsize_hsh[tname] = tsize
    bamfile.close()
    return tname2tsize_hsh

