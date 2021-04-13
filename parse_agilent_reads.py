import sys
import pysam
import re
import parasail
from Bio import pairwise2
from Bio.Seq import Seq

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
out_file = sys.argv[2]
left_flank_file = sys.argv[3]
right_flank_file = sys.argv[4]

out_file_fd = open(out_file,'w')

#for read in pysam.FastxFile(read_file):
#    print(read.sequence)


#left_flank = "CGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGTCTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAACTGTCAGACCAAGTTTACTCATATATACTTTAGATTGATTTAAAACTTCATTTTTAATTTAAAAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGACCAAAATC"
#for i in range(0,len(left_flank)):
#    c1 = 0
#    c2 = 0
#    c3 = 0
#    flank = left_flank[i:len(left_flank)]
#    for read in pysam.FastxFile(read_file):
#        seq = read.sequence
        #if seq.find("CGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGTCTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAACTGTCAGACCAAGTTTACTCATATATACTTTAGATTGATTTAAAACTTCATTTTTAATTTAAAAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGACCAAAATC") != -1:
#        if seq.find(flank) != -1:
#            if seq.find("CCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCT") != -1:
#                oligo = seq[seq.find(flank):seq.find("CCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCT")]
#                c1 = c1 + 1
                #print(seq.find("CCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCT"))
                #print( seq.find("CGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGC"))
#                print(oligo)
#            c3 = c3 + 1
#        else:
#            c2 = c2 + 1
#    print(flank)
#    print("Not Found:%d\nFound both:%d\nFound:%d\n"%(c2,c1,c3))
#    if c3 != 0:
#        break

#for read in pysam.FastxFile(read_file):
#    seq = read.sequence
#    print(">%s"%(read.name))
#    c1 = 0
#    c2 = 0
#    c3 = 0
#    left_flank = "CGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGTCTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAACTGTCAGACCAAGTTTACTCATATATACTTTAGATTGATTTAAAACTTCATTTTTAATTTAAAAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGACCAAAATC"
#    for i in range(0,len(left_flank)):
#        lflank = left_flank[i:len(left_flank)]
#        if seq.find(lflank) != -1:
#            c3 = c3 + 1
#            break
#        else:
#            c2 = c2 + 1
    #print(seq.count(lflank))
    #print("Not Found:%d\nFound both:%d\nFound:%d\n"%(c2,c1,c3))
 #   c1 = 0
#    c2 = 0
#    c3 = 0
#    right_flank = "CCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCCGAAGGTAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTTCTTCTAGTGTAGCCGTAGTTAGGCCACCACTTCAAGAACTCTGTAGCACCGCCTACATACCTCGCTCTGCTAATCCTGTTACCAGTGGCTGCTGCCAGTGGCGATAAGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAGCAACATTGAAAAAGGAAGAGTATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTATTGA"
#    for i in range(0,len(right_flank)):
#        rflank = right_flank[0:len(right_flank)-i]
#        if seq.find(rflank) != -1:
#            c3 = c3 + 1
#            break
#        else:
#            c2 = c2 + 1
    #print(seq.count(rflank))
#    #print("Not Found:%d\nFound both:%d\nFound:%d\n"%(c2,c1,c3))
#    idxlflank = [i.start() for i in re.finditer(lflank, left_flank)] 
#    print(idxlflank)
#    idxrflank = [i.start() for i in re.finditer(rflank, right_flank)]
#    print(idxrflank)
#    print(seq[seq.find(lflank):seq.find(rflank)])

for read in pysam.FastxFile(left_flank_file):
    left_flank = read.sequence

for read in pysam.FastxFile(right_flank_file):
    right_flank = read.sequence

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
#right_flank = "CGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGTCTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAACTGTCAGACCAAGTTTACTCATATATACTTTAGATTGATTTAAAACTTCATTTTTAATTTAAAAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGACCAAAATC"
#left_flank = "CCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCCGAAGGTAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTTCTTCTAGTGTAGCCGTAGTTAGGCCACCACTTCAAGAACTCTGTAGCACCGCCTACATACCTCGCTCTGCTAATCCTGTTACCAGTGGCTGCTGCCAGTGGCGATAAGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAGCAACATTGAAAAAGGAAGAGTATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTATTGA"   
scoring_matrix = parasail.matrix_create("ACGT", 5, -1)
print("read_name\tstrand\t5p_flank_start\t5p_flank_end\t5p_flank_percent_identity\t3p_flank_start\t3p_flank_end\t3p_flank_percent_identity\toligo_sequence\tleft_flank_missing\tright_flank_missing\toligo_missing\tall_good")
for read in pysam.FastxFile(read_file):
    seq = read.sequence
    left_flank_missing = 0
    right_flank_missing = 0
    oligo_missing = 0
    all_good = 0
    pi_left = 0
    pi_right = 0
    oligo = ""
    left_begin = 0 
    right_begin = 0
    left_end = 0
    right_end = 0
    strand = "+"
    if len(seq) < 1000:
        continue

    reverse_complement = "".join(complement.get(base, base) for base in reversed(read.sequence))
    #seq1 = Seq(left_flank)
    #seq2 = Seq(read.sequence)
    #test_alignments = pairwise2.align.localds(seq1, seq2, parasail.blosum62, -10, -1)
    #print(test_alignments)
    result_left = parasail.ssw(left_flank, read.sequence, 5, 4, scoring_matrix)
    result_left_traceback = parasail.sw_trace_striped_16(left_flank, read.sequence, 5, 4, scoring_matrix)
    result_right = parasail.ssw(right_flank, read.sequence, 5, 4, scoring_matrix)
    result_right_traceback = parasail.sw_trace_striped_16(right_flank, read.sequence, 5, 4, scoring_matrix)

    result_left_reverse = parasail.ssw(left_flank, reverse_complement, 5, 4, scoring_matrix)
    result_left_traceback_reverse = parasail.sw_trace_striped_16(left_flank, reverse_complement, 5, 4, scoring_matrix)
    result_right_reverse = parasail.ssw(right_flank, reverse_complement, 5, 4, scoring_matrix)
    result_right_traceback_reverse = parasail.sw_trace_striped_16(right_flank, reverse_complement, 5, 4, scoring_matrix)
    #traceback = result.traceback
    #traceback_right = result_right.traceback
    if percentage_identity(result_left_traceback.traceback.comp) > 0.6 and percentage_identity(result_right_traceback.traceback.comp) > 0.6 and len(seq[result_left.ref_end1:result_right.ref_begin1]) > 50 and len(seq[result_left.ref_end1:result_right.ref_begin1]) < 180:
        #print(read.name)
        #print(seq[result_left.ref_end1])
        #print(left_flank[result_left.read_end1])
        #print(result_left_traceback.begin)
        #print(result_left.ref_begin1)
        #print(result_left.ref_end1)
        #print(result_left.read_begin1)
        #print(result_left.read_end1)
        #cigar = result_left.cigar
        #print(result_left.cigar.decode)    
        #print(result_left_traceback.traceback.comp)
        #print(percentage_identity(result_left_traceback.traceback.comp))
        #print(result_right.ref_begin1)
        #print(result_right.ref_end1)
        #print(result_right.read_begin1)
        #print(result_right.read_end1)
        #print(result_right_traceback.traceback.comp)
        #print(percentage_identity(result_right_traceback.traceback.comp))
        all_good = 1
        strand = "+"
        out_file_fd.write(">%s\n"%(read.name))
        out_file_fd.write("%s\n"%(seq[result_left.ref_end1:result_right.ref_begin1]))
        oligo = seq[result_left.ref_end1:result_right.ref_begin1]
        pi_left = percentage_identity(result_left_traceback.traceback.comp)
        pi_right = percentage_identity(result_right_traceback.traceback.comp)
        left_begin = result_left.ref_begin1
        left_end = result_left.ref_end1
        right_begin = result_right.ref_begin1
        right_end = result_right.ref_end1
        #print("%s\t%d\t%d\t%f\t%d\t%d\t%f\t%s"%(read.name,result_left.ref_begin1,result_left.ref_end1,percentage_identity(result_left_traceback.traceback.comp),result_right.ref_begin1,result_right.ref_end1,percentage_identity(result_right_traceback.traceback.comp),seq[result_left.ref_end1:result_right.ref_begin1]))
        #print("\n")
    elif percentage_identity(result_left_traceback_reverse.traceback.comp) > 0.6 and percentage_identity(result_right_traceback_reverse.traceback.comp) > 0.6 and len(reverse_complement[result_left_reverse.ref_end1:result_right_reverse.ref_begin1]) > 50 and len(reverse_complement[result_left_reverse.ref_end1:result_right_reverse.ref_begin1]) < 180:
        all_good = 1
        strand = "-"
        out_file_fd.write(">%s\n"%(read.name))
        out_file_fd.write("%s\n"%(reverse_complement[result_left_reverse.ref_end1:result_right_reverse.ref_begin1]))
        oligo = reverse_complement[result_left_reverse.ref_end1:result_right_reverse.ref_begin1]
        pi_left = percentage_identity(result_left_traceback_reverse.traceback.comp)
        pi_right = percentage_identity(result_right_traceback_reverse.traceback.comp)
        left_begin = result_left_reverse.ref_begin1
        left_end = result_left_reverse.ref_end1
        right_begin = result_right_reverse.ref_begin1
        right_end = result_right_reverse.ref_end1

    if all_good == 0:
        pi_left = percentage_identity(result_left_traceback.traceback.comp)
        pi_right = percentage_identity(result_right_traceback.traceback.comp)
        left_begin = result_left.ref_begin1
        left_end = result_left.ref_end1
        right_begin = result_right.ref_begin1
        right_end = result_right.ref_end1
        if percentage_identity(result_left_traceback.traceback.comp) < 0.6 or percentage_identity(result_left_traceback_reverse.traceback.comp) < 0.6 :
            left_flank_missing = 1
        if percentage_identity(result_right_traceback.traceback.comp) < 0.6 or percentage_identity(result_right_traceback_reverse.traceback.comp) < 0.6:
            right_flank_missing = 1
        if len(seq[result_left.ref_end1:result_right.ref_begin1]) < 50 or len(seq[result_left.ref_end1:result_right.ref_begin1]) > 180 or len(reverse_complement[result_left_reverse.ref_end1:result_right_reverse.ref_begin1]) < 50 or len(reverse_complement[result_left_reverse.ref_end1:result_right_reverse.ref_begin1]) > 180:
            oligo_missing = 1
    print("%s\t%s\t%d\t%d\t%f\t%d\t%d\t%f\t%s\t%d\t%d\t%d\t%d"%(read.name,strand,left_begin,left_end,pi_left,right_begin,right_end,pi_right,oligo,left_flank_missing,right_flank_missing,oligo_missing,all_good))
    #print(read.name)
    #print(len(traceback.query))
    #print(result.score)
    #print(traceback.ref)
    #print(traceback.comp)
    #print(traceback.query)
    #print(len(result1.traceback.query))
    #print(result1.score)
    #print(result1.traceback.ref)
    #print(result1.traceback.comp)
    #print(result1.traceback.query)


