import sys

# Post-processing of a kfinger output obtained from non-combined fingerprints (with read duplication)

# kfinger coordinates are 0-based (read indices are 1-based)

# A flag "_1" at the end of a read index means "original read", otherwise a flag "_0" means the r-c version of the read.

finger_alignments = sys.argv[1]

check_consistent_overlaps = bool(int(sys.argv[2]))  # 0

strict_consistency = bool(int(sys.argv[3]))     # 0

# For a pair of reads (in the original set), the first overlap encountered (in their duplication) is retained if check_consistent_overlaps is False,
# otherwise the consistency is checked. Only pairs of reads aligned 0,0 and 1,1 but not (1,0) or (0,1), OR aligned 1,0 and 0,1 but not (1,1) or (0,0) are retained.
# A weak consinstency accepts even one between 0,0 and 1,1 (with non 1,0 and 0,1), or 1,0 and 0,1 (with non 0,0 or 1,1).

# Output path
output_file_name = sys.argv[4]

with open(finger_alignments,'r') as finger_input:
    finger_input_list = finger_input.readlines()

if finger_input_list[-1][:5] == 'Time:':
    finger_input_list.pop(-1)

parameter_records = finger_input_list[:finger_input_list.index("//\n")]

finger_input_list = finger_input_list[finger_input_list.index("//\n")+1:]

output_file = open(output_file_name, "w")

for par in parameter_records:
    output_file.write(par)

output_file.write('# POST-PROCESSING PARAMETERS:\n')
output_file.write('# Script: ' + sys.argv[0] + "\n")
output_file.write('# Overlap file: ' + str(finger_alignments) + "\n")
output_file.write('# Check consistency: ' + str(check_consistent_overlaps) + "\n")
if check_consistent_overlaps:
    output_file.write('\t# Strict_consistency: ' + str(check_consistent_overlaps) + "\n")

output_file.write('# Output file: ' + str(output_file_name) + "\n")
output_file.write("//" + "\n")

# The first and the seventh fields of a KFINGER record are the ids of the two overlapping reads in the input file (0_based index + underscore + duplication flag). The 6-th field if the strand flag of the first read (always 0 if the fingerprints are not combined), the 12-th field is the strand flag of the second read (0 in any case).

finger_dict = dict()

min_retained_overlap_length = -1
max_retained_overlap_length = 0
mean_retained_overlap_length = 0

count_retained = 0
total_retained_length= 0

for alignment in finger_input_list:
    record = alignment.split()
    
    #Indices are in the form [index]_[flag]
    #[index] is the 0-based index of the read
    # NB: [flag] 1 refers to the original read, [flag] 0 refers to the reverse version
    
    first_index = "_".join(record[0].split('_')[:-1])
    
    first_rank = int(record[1])
    if first_rank % 2 != 0:
        first_rank = first_rank // 2 + 1
    else:
        first_rank = first_rank // 2

    second_index = "_".join(record[6].split('_')[:-1])

    second_rank = int(record[7])
    if second_rank % 2 != 0:
        second_rank = second_rank // 2 + 1
    else:
        second_rank = second_rank // 2

    alignments = finger_dict.get((first_index, second_index), [])
    
    #alignments.append(alignment)
    
    #alignments.append([record[0].split('_')[1], record[6].split('_')[1], alignment])
    alignments.append([record[0].split('_')[-1], record[6].split('_')[-1], alignment, first_rank, second_rank])
    
    finger_dict[(first_index, second_index)] = alignments


for key in finger_dict:
    #first_index = int(key[0])
    #second_index = int(key[1])
    first_index = key[0]
    second_index = key[1]

    records = finger_dict[key][0][2].split()
        
    ok = False

    first_rank = finger_dict[key][0][3]
    second_rank = finger_dict[key][0][4]

    # When consistency is checked, strands (0,0) are compatible with (1,1) and strands (1,0) are compatible with (0,1) (and viceversa)
    # NB: 1 refers to the original read, 0 refers to the reverse version
    if check_consistent_overlaps == True:
        strand_set = set([(al_list[0], al_list[1]) for al_list in finger_dict[key]])
        if strict_consistency == True:
            if ('0','0') in strand_set and ('1','1') in strand_set:
                if not(('0','1') in strand_set or ('1','0') in strand_set):
                    ok = True
            if ('0','1') in strand_set and ('1','0') in strand_set:
                if not(('0','0') in strand_set or ('1','1') in strand_set):
                    ok = True
        else:
            if ('0','0') in strand_set or ('1','1') in strand_set:
                if not(('0','1') in strand_set or ('1','0') in strand_set):
                    ok = True
            if ('0','1') in strand_set or ('1','0') in strand_set:
                if not(('0','0') in strand_set or ('1','1') in strand_set):
                    ok = True
    else:
        ok = True

    if ok == True:
        #MODIFICATO QUI
        # Original read
        if str(finger_dict[key][0][0]) == '1':
            records[5] = str(0)
        else:
            records[5] = str(1)
                
        if str(finger_dict[key][0][1]) == '1':
            records[11] = str(0)
        else:
            records[11] = str(1)

        records[0] = str(first_index)

        records[1] = str(first_rank)

        records[6] = str(second_index)

        records[7] = str(second_rank)

        start1 = int(records[3])
        end1 = int(records[4])

        start2 = int(records[9])
        end2 = int(records[10])

        length1 = int(records[2])
        length2 = int(records[8])

        if records[5] == '1':
            records[3] = str(length1-end1-1)
            records[4] = str(length1-start1-1)
           
        if records[11] == '1':
            records[9] = str(length2-end2-1)
            records[10] = str(length2-start2-1)

        output_file.write('\t'.join(records) + "\n")

        overlap_length = end1-start1+1

        count_retained = count_retained + 1
            
        if min_retained_overlap_length == -1 or min_retained_overlap_length > overlap_length:
            min_retained_overlap_length = overlap_length

        if max_retained_overlap_length < overlap_length:
            max_retained_overlap_length = overlap_length
                    
        total_retained_length = total_retained_length + overlap_length


mean_retained_overlap_length = total_retained_length/count_retained

sys.stderr.write("Minimum retained length: " + str(min_retained_overlap_length) + "\n")
sys.stderr.write("Maximum retained length: " + str(max_retained_overlap_length) + "\n")
sys.stderr.write("Mean retained length: " + str(mean_retained_overlap_length) + "\n")

