import math
import sys

import edlib

# File of the matchings/overlaps
overlap_file = sys.argv[1]

# FASTA file of the input reads
input_reads = sys.argv[2]

# Error threshold between good ann bad matchings/overlaps
limit_error = float(sys.argv[3])    # 3.0

def rev_comp(sequence):
    conv_dict = {'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C'}
    return ''.join(conv_dict[c] for c in sequence.rstrip()[::-1])

# Coordinates 0-based
# Read indexes 1-based

with open(input_reads,'r') as read_input:
    read_input_list = read_input.readlines()

with open(overlap_file,'r') as overlap_input:
    overlap_list = overlap_input.readlines()

read_input_list = read_input_list[1::2]

parameter_records = overlap_list[:overlap_list.index("//\n")]

for par in parameter_records:
    print(par.rstrip())

print('# CHECKING PARAMETERS:')
print('# Script: ' + sys.argv[0])
print('# Overlap file: ' + str(overlap_file))
print('# Read file: ' + str(input_reads))
print("//")

count_0err = 0
count_ok = 0
count_no = 0

count = 0

min_overlap_length = -1
max_overlap_length = 0
total_length = 0

equal_reads = 0

overlap_list = overlap_list[overlap_list.index("//\n")+1:]

for record in overlap_list:
    
    record = record.rstrip()
    
    sys.stderr.write(record + "\n")
    
    record = record.split()
    
    if record[1] != record[7]:
        read1 = read_input_list[int(record[1])-1]
        read2 = read_input_list[int(record[7])-1]

        start1 = int(record[3])
        end1 = int(record[4])
    
        start2 = int(record[9])
        end2 = int(record[10])
 
        ov1 = read1[start1:end1+1]
        ov2 = read2[start2:end2+1]
   
        strand1 = record[5]
        strand2 = record[11]

        if strand1 == '1':
            ov1 = rev_comp(ov1)
            
        if strand2 == '1':
            ov2 = rev_comp(ov2)

        if max(len(ov1), len(ov2)) > max_overlap_length:
            max_overlap_length = max(len(ov1), len(ov2))

        if min_overlap_length == -1 or max(len(ov1), len(ov2)) < min_overlap_length:
            min_overlap_length = max(len(ov1), len(ov2))
    
        total_length = total_length + max(len(ov1), len(ov2))

        err = 0
        if ov1 != ov2:
            ed_res = edlib.align(
                ov1, ov2,
                k = int(math.ceil( max(len(ov1), len(ov2)) * 2 * limit_error / 100))
            )
            err = ed_res['editDistance'] / max(len(ov1), len(ov2))*100

        if err == 0:
            count_0err = count_0err + 1
        elif ed_res['editDistance'] != -1 and err <= limit_error:
            count_ok = count_ok + 1
        else:
            count_no = count_no + 1

        count = count + 1
    else:
        equal_reads = equal_reads + 1

mean_length = total_length / count

print("Zero error: " + str(count_0err))
print("<= " + str(limit_error) + " error: " + str(count_ok))
print("> " + str(limit_error) + " error: " + str(count_no))
print("Minimum matching length: " + str(min_overlap_length))
print("Maximum matching length: " + str(max_overlap_length))
print("Mean matching length: " + str(mean_length))
print("Pairs of the same reads: " + str(equal_reads))


