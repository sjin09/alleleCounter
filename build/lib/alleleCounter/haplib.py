import bisect
import natsort
import itertools
from typing import Set, Dict, List, Tuple

# global
BIT_LST = ["0", "1", "-"]
BIT_COMB = itertools.product(BIT_LST, repeat=2)
BIT_COMB_LST = ["{}{}".format(i, j) for i,j in BIT_COMB]
HBIT_COMPLEMENT_MAP = {
	"0": "1",
	"1": "0",
	"-": "0",
	"00": "11",
	"01": "10",
	"0-": "10",
	"10": "01",
	"11": "00",
	"1-": "00",
	"-0": "01",
	"-1": "00",
	"--": "00"
}
ALT_HBIT_COMPLEMENT_MAP = {
	"0": ["1", "2", "-"],
	"1": ["0", "2", "-"],
	"-": ["0", "1", "2"],
	"00": ['1-', '11', '-1', '--'],
	"01": ['1-', '10', '-0', '--'], 
	"0-": ['10', '11', '-0', '-1'],
	"10": ['0-', '01', '-1', '--'],
	"11": ['00', '0-', '-0', '--'],
	"1-": ['00', '01', '-0', '-1'],
	"-0": ['0-', '01', '1-', '11'],
	"-1": ['00', '0-', '1-', '10'],
	"--": ['00', '01', '10', '11']
}


def get_hapbits():
	bit_lst = ["0", "1", "-"]
	bit_comb_lst = []
	bit_comb = itertools.product(bit_lst, repeat=2)
	for i,j in bit_comb:
		bit_comb_lst.append("{}{}".format(i,j))
	return bit_lst, bit_comb_lst


def get_hbit_complement(hbit):
	hbit_complement = HBIT_COMPLEMENT_MAP[hbit]
	return hbit_complement


def get_alt_hbit_complement(alt_hbit, hbit_count_hsh):
	ref_hbit_candidate = HBIT_COMPLEMENT_MAP[alt_hbit]
	ref_hbit_candidate_lst = ALT_HBIT_COMPLEMENT_MAP[alt_hbit]	
	ref_hbit_candidate_count_lst = [hbit_count_hsh[ref_hbit] for ref_hbit in ref_hbit_candidate_lst]
	ref_hbit_max_count = max(ref_hbit_candidate_count_lst)
	if ref_hbit_max_count == 0:
		return 0, ref_hbit_candidate, 0
	else:
		counter = sum([1 for j in ref_hbit_candidate_count_lst if j == ref_hbit_max_count])
		if counter == 1:
			index = ref_hbit_candidate_count_lst.index(ref_hbit_max_count)
			ref_hbit_candidate = ref_hbit_candidate_lst[index]
			return 1, ref_hbit_candidate, ref_hbit_max_count
		else:
			index_lst = [i for i, j in enumerate(ref_hbit_candidate_count_lst) if j == ref_hbit_max_count]
			ref_hbit_candidate = ref_hbit_candidate_lst[index_lst[0]]
			return 0, ref_hbit_candidate, ref_hbit_max_count


def get_tsub_adjacet_hetsnps(pos, tsub_target_read_lst, read_tpos_alleles_map, hetsnps, hetsnp_count, hetsnp_positions):

	target_hetsnp2count_map = {} 
	tsub_index = bisect.bisect_left(hetsnp_positions, pos)
	for target_read in tsub_target_read_lst:
		state = 1
		iteration = 0
		target_tsub_index = tsub_index
		target_read_tpos_lst = list(read_tpos_alleles_map[target_read].keys())
		target_read_start = target_read_tpos_lst[0]
		target_read_end = target_read_tpos_lst[-1]
		while state:
			if target_tsub_index == 0:
				i = 0
				j = 1
			elif target_tsub_index == hetsnp_count:
				i = target_tsub_index - 2
				j = target_tsub_index - 1
			else:
				i = target_tsub_index - 1
				j = target_tsub_index
			upstream_hetsnp = hetsnps[i]
			downstream_hetsnp = hetsnps[j]
			upos = upstream_hetsnp[1]
			dpos = downstream_hetsnp[1]
			if upos > target_read_end: # read doesn't overlap or span the hetsnp
				break # finish search
			elif target_read_start < upos and dpos < target_read_end: # read spans both upstream and downstream hetsnps
				for snp in [upstream_hetsnp, downstream_hetsnp]:
					if snp not in target_hetsnp2count_map:
						target_hetsnp2count_map[snp] = 1
					else:
						target_hetsnp2count_map[snp] += 1
				break # finish search
			else: # continue search
				iteration += 1
				target_tsub_index += 1 # sliding window
				if target_read_start < upos and upos < target_read_end: # read overlaps with upstream hetsnp
					if upstream_hetsnp not in target_hetsnp2count_map:
						target_hetsnp2count_map[upstream_hetsnp] = 1
					else:
						target_hetsnp2count_map[upstream_hetsnp] += 1
			if iteration > 10:
				state = 0
				break

	# # return
	target_hetsnp_lst = list(target_hetsnp2count_map.keys())
	target_hetsnp_candidate_count = len(target_hetsnp_lst)
	if target_hetsnp_candidate_count == 0: # reads do not overlap with hetsnps
		return 0, 0, 0
	elif target_hetsnp_candidate_count == 1:
		upstream_hetsnp = target_hetsnp_lst[0]
		return 1, upstream_hetsnp, 0
	elif target_hetsnp_candidate_count > 1:
		consensus_count = max(target_hetsnp2count_map.values())
		target_hetsnp_consensus_lst = [target_hetsnp for target_hetsnp, counts in target_hetsnp2count_map.items() if counts == consensus_count]
		target_hetsnp_consensus_lst = natsort.natsorted(target_hetsnp_consensus_lst)
		target_hetsnp_consensus_count = len(target_hetsnp_consensus_lst)
		if target_hetsnp_consensus_count == 1: # upstream hetsnp is present within the read
			upstream_hetsnp = target_hetsnp_lst[0]
			return 1, upstream_hetsnp, 0
		elif target_hetsnp_consensus_count > 1: # upstream and downstream hetsnps are present within the reads
			upstream_hetsnp, downstream_hetsnp = target_hetsnp_lst[0:2]
		return 2, upstream_hetsnp, downstream_hetsnp


def get_read_hbits(
	read_lst: List[str], 
	hetsnp_state: int, 
	upstream_hetsnp: Tuple[str, int, str, str], 
	downstream_hetsnp: Tuple[str, int, str, str], 
	read_tpos_alleles_map: Dict[str, Dict[int, Tuple[str, str, int]]]
) -> Tuple[Dict[str, Tuple[str, str]], Dict[Tuple[str, str], int]]:

	read_hbits_map = {}
	_chrom, upos, uref, ualt = upstream_hetsnp
	uhet = (uref, ualt)
	uhom = (uref, uref)
	hbit_count_map = {bit: 0 for bit in BIT_LST} if hetsnp_state == 1 else {bits: 0 for bits in BIT_COMB_LST}
	for qread in read_lst:
		if qread not in read_tpos_alleles_map:
			continue
		else:
			tpos_alleles_map = read_tpos_alleles_map[qread]
			if upos not in tpos_alleles_map: 
				ubit = '-'
			else:
				qread_usnp = tpos_alleles_map[upos][0:2]
				if qread_usnp == uhet: # upstream hetsnp
				    ubit = "1"
				elif qread_usnp == uhom: # upstream honsnp
				    ubit = "0"
				else:
					ubit = "-"
		if hetsnp_state == 1:
			hbit_count_map[ubit] += 1
			read_hbits_map[qread] = ubit
		elif hetsnp_state == 2:
			_chrom, dpos, dref, dalt = downstream_hetsnp
			dhet = (dref, dalt)
			dhom = (dref, dref)
			if dpos not in tpos_alleles_map:
				dbit = "-"
			else:
				qread_dsnp = tpos_alleles_map[dpos][0:2]
				if qread_dsnp == dhet: # downstream snp
				    dbit = "1"
				elif qread_dsnp == dhom:
				    dbit = "0"
				else:
					dbit = "-"
			read_hbits = "{}{}".format(ubit, dbit)
			hbit_count_map[read_hbits] += 1
			read_hbits_map[qread] = read_hbits
	return read_hbits_map, hbit_count_map

# def get_target_hbit_count(target_read_lst, read_hbit_hsh, hbit_count_hsh):
# 	target_hbit_lst = [read_hbit_hsh[read] for read in target_read_lst]
# 	target_uniq_hbit_count = len(list(set(target_hbit_lst)))
# 	if target_uniq_hbit_count == 1: # single alt haplotype
# 		target_hbit = target_hbit_lst[0]
# 		target_hbit_count = len(target_read_lst)
# 		alt_hbit_count = hbit_count_hsh[target_hbit] - target_hbit_count
# 		return 1, target_hbit, alt_hbit_count
# 	else: # multiple alt haplotypes
# 		return 0, 0, 0
		