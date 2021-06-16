import Levenshtein
import sys

# File of the matchings/overlaps
mini_file = sys.argv[1]

# FASTA file of the input reads
input_reads = sys.argv[2]

# Error threshold between good ann bad matchings/overlaps
limit_error = float(sys.argv[3])    # 3.0

def rev_comp(sequence):
    conv_dict = {'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C'}
    return ''.join(conv_dict[c] for c in sequence.rstrip()[::-1])

# Coordinates 0-based
# Read indexes 1-based

print('# CHECKING PARAMETERS:')
print('# Script: ' + sys.argv[0])
print('# Minimap file: ' + str(mini_file))
print('# Read file: ' + str(input_reads))
print("//")

with open(input_reads,'r') as read_input:
    read_input_list = read_input.readlines()

with open(mini_file,'r') as overlap_input:
    overlap_list = overlap_input.readlines()

header_input_list = read_input_list[0::2]
read_input_list = read_input_list[1::2]

header_input_list = [h.rstrip()[1:] for h in header_input_list]
read_input_list = [r.rstrip() for r in read_input_list]

count_0err = 0
count_ok = 0
count_no = 0

count = 0

min_overlap_length = -1
max_overlap_length = 0
total_length = 0
equal_reads = 0

dupl_set = set()
present = 0
present_0 = 0
present_ok = 0
present_no = 0

for record in overlap_list:
    
    record = record.rstrip()
    
    sys.stderr.write(record + "\n")

    record = record.split()
    
    if record[0] != record[5]:
        
        pair = (record[0], record[5])
        pair = tuple(sorted(pair))
        
        is_present = False
        if pair in dupl_set:
            present = present + 1
            is_present = True
        
        dupl_set.add(pair)
        
        read1 = read_input_list[header_input_list.index(record[0])]
        read2 = read_input_list[header_input_list.index(record[5])]

        start1 = int(record[2])
        end1 = int(record[3])
    
        start2 = int(record[7])
        end2 = int(record[8])
    
        ov1 = read1[start1:end1]
        ov2 = read2[start2:end2]
    
        if max(len(ov1), len(ov2)) > max_overlap_length:
            max_overlap_length = max(len(ov1), len(ov2))

        if min_overlap_length == -1 or max(len(ov1), len(ov2)) < min_overlap_length:
            min_overlap_length = max(len(ov1), len(ov2))
    
        total_length = total_length + max(len(ov1), len(ov2))
   
        strand = record[4]

        if strand == '-':
            ov2 = rev_comp(ov2)
            
        err = 0
        if ov1 != ov2:
            err = Levenshtein.distance(ov1, ov2)/max(len(ov1), len(ov2))*100

        if err == 0:
            count_0err = count_0err + 1
            if is_present:
                present_0 = present_0 + 1
        elif err <= limit_error:
            count_ok = count_ok + 1
            if is_present:
                present_ok = present_ok + 1
        else:
            count_no = count_no + 1
            if is_present:
                present_no = present_no + 1

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
print("Duplicated: " + str(present))
print("Duplicated_0: " + str(present_0))
print("Duplicated_ok: " + str(present_ok))
print("Duplicated_no: " + str(present_no))

