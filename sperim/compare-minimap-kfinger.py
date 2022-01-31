# Confronto minimap-kfinger (PAF files)

import Levenshtein
import sys

kfinger_file = sys.argv[1]
mini_file = sys.argv[2]
input_reads = sys.argv[3]

def rev_comp(sequence):
    conv_dict = {'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C'}
    return ''.join(conv_dict[c] for c in sequence.rstrip()[::-1])

# Le coordinate sono 0-based
# Gli indici dei reads sono 1-based

print('# COMPARING PARAMETERS:')
print('# Script: ' + sys.argv[0])
print('# Kfinger file: ' + str(kfinger_file))
print('# Minimap file: ' + str(mini_file))
print('# Read file: ' + str(input_reads))
print("//")

with open(input_reads,'r') as read_input:
    read_input_list = read_input.readlines()

with open(kfinger_file,'r') as k_overlap_input:
    k_overlap_list = k_overlap_input.readlines()

with open(mini_file,'r') as m_overlap_input:
    m_overlap_list = m_overlap_input.readlines()

header_input_list = read_input_list[0::2]
read_input_list = read_input_list[1::2]

header_input_list = [h.rstrip()[1:] for h in header_input_list]
read_input_list = [r.rstrip() for r in read_input_list]

k_set = set()

ov_list = k_overlap_list

for record in ov_list:
    
    record = record.rstrip()
    record = record.split()
    
    if record[0] != record[5]:
        pair = (record[0], record[5])
        pair = tuple(sorted(pair))
        k_set.add(pair)

ov_list = m_overlap_list

count = 0

for record in ov_list:
    
    record = record.rstrip()
    record1 = record.split()
    
    if record1[0] != record1[5]:
        pair = (record1[0], record1[5])
        pair = tuple(sorted(pair))
        if not(pair in k_set):
            read1 = read_input_list[header_input_list.index(record1[0])]
            read2 = read_input_list[header_input_list.index(record1[5])]
        
            start1 = int(record1[2])
            end1 = int(record1[3])
        
            start2 = int(record1[7])
            end2 = int(record1[8])
        
            ov1 = read1[start1:end1]
            ov2 = read2[start2:end2]
        
            strand = record1[4]
        
            if strand == '-':
                ov2 = rev_comp(ov2)

            err = 0
            if ov1 != ov2:
                err = Levenshtein.distance(ov1, ov2)/max(len(ov1), len(ov2))*100
            
            if err == 0:
                count = count + 1
                print(record)
                print(str(end2-start2))

sys.stderr.write("Count: " + str(count) + "\n")
      

