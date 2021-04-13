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

input_file = sys.argv[1]
left_file = sys.argv[2]
right_file = sys.argv[3]

for read in pysam.FastxFile(left_file):
    left = read.sequence

for read in pysam.FastxFile(right_file):
    right = read.sequence

print(right)
input_fd = open(input_file)

scoring_matrix = parasail.matrix_create("ACGT", 5, -1)

header = input_fd.readline()
c = 0
tot = 0
for line in input_fd:
    pi = line.rstrip().split()[3]
    pi_reverse = line.rstrip().split()[6]
    left_flank = line.rstrip().split()[7]
    right_flank = line.rstrip().split()[8]
    if float(pi) > 0.6 or float(pi_reverse) > 0.6:
        result_left = parasail.ssw(left_flank, left, 5, 4, scoring_matrix)
        result_left_traceback = parasail.sw_trace_striped_16(left_flank, left, 5, 4, scoring_matrix)
        result_right = parasail.ssw(right_flank, right, 5, 4, scoring_matrix)
        result_right_traceback = parasail.sw_trace_striped_16(right_flank, right, 5, 4, scoring_matrix)
        c = c + 1
        #print("%s\t%f\t%s\t%f"%(left_flank, percentage_identity(result_left_traceback.traceback.comp), right_flank, percentage_identity(result_right_traceback.traceback.comp)))
        if percentage_identity(result_left_traceback.traceback.comp) > 0.6 and percentage_identity(result_right_traceback.traceback.comp) > 0.6:
            tot = tot + 1

    
print(c)
print(tot)
#print("%s\t%s" % (pi, pi_reverse))