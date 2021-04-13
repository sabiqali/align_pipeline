import sys
import pysam
import parasail

def percentage_identity(cigar_exp):
    m = 0
    nm = 0
    for chr in cigar_exp:
        if chr == '|':
            m = m + 1
        else:
            nm = nm + 1
    pi = float(m/(nm + m))
    return pi

read_file = sys.argv[1]
left_flank_file = sys.argv[2]
right_flank_file = sys.argv[3]

for read in pysam.FastxFile(left_flank_file):
    left_flank = read.sequence

for read in pysam.FastxFile(right_flank_file):
    right_flank = read.sequence

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
scoring_matrix = parasail.matrix_create("ACGT", 5, -1)
left_flank_missing = 0
right_flank_missing = 0
oligo_missing = 0
all_good = 0
for read in pysam.FastxFile(read_file):
    seq = read.sequence
    reverse_complement = "".join(complement.get(base, base) for base in reversed(read.sequence))
    if len(read.sequence) > 1000 and len(read.sequence) < 3000:
        result_left = parasail.ssw(left_flank, read.sequence, 5, 4, scoring_matrix)
        result_left_traceback = parasail.sw_trace_striped_16(left_flank, read.sequence, 5, 4, scoring_matrix)
        result_right = parasail.ssw(right_flank, read.sequence, 5, 4, scoring_matrix)
        result_right_traceback = parasail.sw_trace_striped_16(right_flank, read.sequence, 5, 4, scoring_matrix)

        result_left_reverse = parasail.ssw(left_flank, reverse_complement, 5, 4, scoring_matrix)
        result_left_traceback_reverse = parasail.sw_trace_striped_16(left_flank, reverse_complement, 5, 4, scoring_matrix)
        result_right_reverse = parasail.ssw(right_flank, reverse_complement, 5, 4, scoring_matrix)
        result_right_traceback_reverse = parasail.sw_trace_striped_16(right_flank, reverse_complement, 5, 4, scoring_matrix)
        #print(result_left_traceback.cigar.beg_ref)
        #print(seq[result_left_traceback.end_ref:result_right_traceback.cigar.beg_ref])
        if percentage_identity(result_left_traceback.traceback.comp) < 0.6:
            left_flank_missing = left_flank_missing + 1
        if percentage_identity(result_right_traceback.traceback.comp) < 0.6:
            right_flank_missing = right_flank_missing + 1
        #if len(seq[result_left.ref_end1:result_right.ref_begin1]) < 50 or len(seq[result_left.ref_end1:result_right.ref_begin1]) > 180:
        #if len(seq[result_left_traceback.end_ref:result_right_traceback.cigar.beg_ref]) < 50 or len(seq[result_left_traceback.end_ref:result_right_traceback.cigar.beg_ref]) > 180:
        if len(seq[result_left.ref_end1:result_right.ref_begin1]) < 50 or len(seq[result_left.ref_end1:result_right.ref_begin1]) > 180:
            oligo_missing = oligo_missing + 1
        if percentage_identity(result_left_traceback.traceback.comp) > 0.6 and percentage_identity(result_right_traceback.traceback.comp) > 0.6 and len(seq[result_left.ref_end1:result_right.ref_begin1]) > 50 and len(seq[result_left.ref_end1:result_right.ref_begin1]) < 180:
            all_good = all_good + 1
        elif percentage_identity(result_left_traceback_reverse.traceback.comp) > 0.6 and percentage_identity(result_right_traceback_reverse.traceback.comp) > 0.6 and len(seq[result_left_reverse.ref_end1:result_right_reverse.ref_begin1]) > 50 and len(seq[result_left_reverse.ref_end1:result_right_reverse.ref_begin1]) < 180:
            all_good = all_good + 1

print("left_flank_missing\tright_flank_missing\toligo_missing\tall_good")
print("%d\t%d\t%d\t%d"%(left_flank_missing,right_flank_missing,oligo_missing,all_good))