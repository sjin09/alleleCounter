import re
from typing import Dict, List, Tuple


def cs2lst(cs_tag: str):
    cslst = [cs for cs in re.split("(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+)", cs_tag)]
    cslst = [cs.upper() for cs in cslst if cs != ""]
    return cslst


def cs2tuple(
    cs_str: str, 
    qpos: int, 
    qseq: str
) -> List[Tuple[int, str, str, int, int]]:

    cstuple_lst = []
    cs_lst = cs2lst(cs_str)
    for cs in cs_lst:
        m = cs[1:]
        mlen = len(m)
        qstart = qpos
        if cs.startswith("="):  # match # --cs=long
            cs = ":{}".format(mlen)
            t = (1, m, m, mlen, mlen)
        elif cs.startswith(":"):  # match # --cs=short
            mlen = int(m)
            qend = qpos + mlen
            m = qseq[qstart:qend]
            t = (1, m, m, mlen, mlen)
        elif cs.startswith("*"):  # snp # target and query
            mlen = 1
            ref, alt = list(m)
            t = (2, ref, alt, 1, 1)
        elif cs.startswith("+"):  # insertion # query
            ref = qseq[qpos - 1] 
            alt = ref + m
            t = (3, ref, alt, 0, mlen)
        elif cs.startswith("-"):  # deletion # target
            alt = qseq[qpos - 1]
            ref = alt + m
            t = (4, ref, alt, mlen, 0)
            mlen = 0
        cstuple_lst.append(t)
        qpos += mlen
    return cstuple_lst


def cs2alleles(
    cstuple_lst: List[Tuple[int, str, str, int, int]], 
    tpos: int, 
    qpos: int, 
    qbq_lst: List[int]
) -> Dict[int, Tuple[str, str, int]]:

    pos2alleles = {}
    for cstuple in cstuple_lst:
        state, ref, alt, ref_len, alt_len, = cstuple
        if state == 1: # match
            for i, (ref_base, alt_base) in enumerate(zip(ref, alt)):
                pos2alleles[tpos + i + 1] = (1, ref_base, alt_base, qbq_lst[qpos + i])
        elif state == 2: # snp
            pos2alleles[tpos + 1] = (2, ref, alt, qbq_lst[qpos])
        elif state == 3: # insertion
            pos2alleles[tpos] = ("1", ref[0], alt[0], qbq_lst[qpos])
        elif state == 4: # deletion
            pos2alleles[tpos] = ("1", ref[0], alt[0], qbq_lst[qpos])
            for j, k in enumerate(ref):
                if j == 0:
                    continue
                pos2alleles[tpos+1] = (4, k, "-", 0)
        tpos += ref_len
        qpos += alt_len
    return pos2alleles