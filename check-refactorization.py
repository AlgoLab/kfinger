#!/usr/bin/env python
# coding: utf-8

import time

import sys
import Levenshtein

from utils import compute_fingerprint
from utils.factorizations import CFL, ICFL_recursive, CFL_icfl, CFL_icfl_cfl
from utils.fingerprint_utils import compute_long_fingerprint_by_list

# PRIMARY PARAMETERS

# Input file of the overlap records
raw_overlap_file = sys.argv[1]

# Input of file of the original fingerprints
input_orig_finger_file = sys.argv[2]

# Input of file of the original factorizations
input_orig_fact_file = sys.argv[3]

# Factorization type
fact_type = sys.argv[4] # CFL_ICFL_COMB-30

if not(fact_type in ('CFL_ICFL-30', 'CFL_ICFL_COMB-30')):
    sys.stderr.write(" Factorization not yet supported for refactorizing!\n")

# k-fingers size
k = int(sys.argv[5]) # 5

# Minimum total length for retaining a k-finger
min_total_length = int(sys.argv[6]) # 10

min_coverage = float(sys.argv[7]) # 0.0

# 1: the suffix-prefix overlaps are produced in output
# 0: the matchings are produced in output
extend_suffix_prefix_overlap = bool(int(sys.argv[8])) # 0

# The following four parameters are active if search_gap is True
downstream_k = int(sys.argv[9])       # 3 for error-free datasets, 2 otherwise
min_downstream_length = int(sys.argv[10])  # 10 for error-free datasets, 20 otherwise
gap_length_tolerance = int(sys.argv[11])   # 0 for error-free datasets, 20 otherwise

# Output path
output_path = sys.argv[12]

# SECONDARY PARAMETERS

noerr_reads = False # False

# Consider k-fingers having at most max_number_of_one 1-long factors (-1 if a k-finger can contain any number of 1s)
max_number_of_one = -1   # -1

# 0, if the matching substrings are limited to the leftmost k-finger extended on the right; not used if noerr_reads is True
search_gap = True   # True

#Separator for the input fingerprints
split_separator = '|'       # '|'

print_error = False         # False

output_file_name = output_path + "refact-matchings.txt"

output_file = open(output_file_name, "w")

def get_matching_prefixes(f1, f2, k, downstream_k, min_downstream_length, search_gap, gap_length_tolerance, multiple_gaps):
    
    #Find the longest common prefix of f1 and f2 (supposing that the first k integers are already identical)
    init_common_pos = k
    while init_common_pos < min(len(f1), len(f2)) and f1[init_common_pos] == f2[init_common_pos]:
        init_common_pos = init_common_pos + 1
    
    #Variable init_common_pos contains the 0-based position after the end of the common prefix
    
    return_list = [init_common_pos, init_common_pos]

    if search_gap == False:
        return return_list

    # SEARCH THE LONGEST GAP

    # Initialize an empty dictionary
    common_right_dict = dict()
    
    #Fill the dictionary with the k-fingers of f1 for k = downstream_k
    for i in range(len(f1))[init_common_pos:len(f1)-downstream_k]:
        if sum(tuple(f1[i:i+downstream_k])) >= min_downstream_length:
            common_right_dict[tuple(f1[i:i+downstream_k])] = i

    # The dictionary contains all the k-fingers of f1 for k = downstream_k (as keys) together their starting positions as values

    # Scan f2 from the end to position init_common_pos in order to search the common k-fingers for k = downstream_k
    i = len(f2)-downstream_k
    found = False
    while i >= init_common_pos and found == False:
        key = tuple(f2[i:i+downstream_k])
        if  key in common_right_dict and abs(sum(f1[:common_right_dict[key]]) - sum(f2[:i])) <= gap_length_tolerance:
            found = True
        else:
            i = i - 1

    if found:
        #Find the longest common factor of f1 and f2 after the gap
        i = i + downstream_k
        p = common_right_dict[key] + downstream_k
        while (p < len(f1) and i < len(f2)) and (f1[p] == f2[i]):
            p = p + 1
            i = i + 1
        
        return_list = [p, i]

        if multiple_gaps:
            k = p-1
            q = i-1
            
            stop = False
            while stop == False:
                while k >= init_common_pos and q >= init_common_pos and (f1[k] == f2[q]):
                    k = k - 1
                    q = q - 1
            
                found2 = False
                while q >= init_common_pos and found2 == False:
                    key = tuple(f2[j:q+downstream_k])
                    if  key in common_right_dict and abs(sum(f1[:common_right_dict[key]]) - sum(f2[:q])) <= gap_length_tolerance:
                        found2 = True
                        k = common_right_dict[key]
                        return_list.append([k,q])
                    
                    q = q - 1
                       
                if k <= init_common_pos or q <= init_common_pos:
                    stop = True
                
    return return_list



with open(raw_overlap_file,'r') as input_overlap:
    record_list = input_overlap.readlines()

parameter_records = record_list[:record_list.index("//\n")]

for par in parameter_records:
    output_file.write(par)

output_file.write('# REFACTORIZATION MATCHING PARAMETERS:\n')
output_file.write('# Script: ' + sys.argv[0] + "\n")
output_file.write('# Raw overlap file: ' + str(raw_overlap_file) + "\n")
output_file.write('# Fingerprint file: ' + str(input_orig_finger_file) + "\n")
output_file.write('# Factorization file: ' + str(input_orig_fact_file) + "\n")
output_file.write('# k: ' + str(k) + "\n")
output_file.write('# Minimum total kfinger length: ' + str(min_total_length) + "\n")
output_file.write('# Maximum number of 1s for considering a k-finger (-1: no checking): ' + str(max_number_of_one) + "\n")
output_file.write('# Minimum matching rate: ' + str(min_coverage) + "\n")
output_file.write('# Extend to suffix-prefix overlap: ' + str(extend_suffix_prefix_overlap) + "\n")

if noerr_reads:
    output_file.write('# Search gap: ' + str(False) + "\n")
else:
    output_file.write('# Search gap: ' + str(search_gap) + "\n")

if search_gap and noerr_reads == False:
    output_file.write('\t# Minimum number of factors required after the gap: ' + str(downstream_k) + "\n")
    output_file.write('\t# Minimum total length of the factors required after the gap: ' + str(min_downstream_length) + "\n")
    output_file.write('\t# Maximum tolerance in the gap length: ' + str(gap_length_tolerance) + "\n")
                  
output_file.write('# Factorization algorithm: ' + str(fact_type) + "\n")
                  
if print_error:
    output_file.write('# Error in the last field\n' + "\n")

output_file.write('# No error reads: ' + str(noerr_reads) + "\n")
output_file.write('# Output file: ' + str(output_file_name) + "\n")

output_file.write("//\n")

overlap_records = record_list[record_list.index("//\n")+1:]

with open(input_orig_finger_file,'r') as input_orig_finger:
    orig_fing_list = input_orig_finger.readlines()

split_length = sum(int(l) for l in " ".join(orig_fing_list[0].split()[1:]).split(split_separator)[0].split())

dict_tuples = [(row.split()[0], row.split()[1:]) for row in orig_fing_list]
orig_finger_dict = dict(dict_tuples)

with open(input_orig_fact_file,'r') as input_orig_fact:
    orig_fact_list = input_orig_fact.readlines()

dict_tuples = [(row.split()[0], row.split()[1:]) for row in orig_fact_list]
orig_fact_dict = dict(dict_tuples)


for i in range(len(overlap_records)):
    
    fields = overlap_records[i].split()
    
    first_read_index = fields[1]
    second_read_index = fields[7]
    
    #sys.stderr.write("Read " + fields[0] + " read " + fields[6] + "\n")
    sys.stderr.write("Read " + first_read_index + " read " + second_read_index + "\n")

    first_read_length = fields[2]
    second_read_length = fields[8]

    sharing = fields[12]
    
    read_ids = (fields[0], fields[6])
    
    cut_read_index = int(fields[13])-1
    cut_prefix = int(fields[14])

    start_index = int(fields[15])
    
    cut_read_seq = "".join(f for f in orig_fact_dict[read_ids[cut_read_index]] if f != '|')
    
    cut_read_seq = cut_read_seq[cut_prefix:]
    
    [fingerprint, l_factorization] = compute_fingerprint(sequence=cut_read_seq, split=split_length, type_factorization=fact_type, fact_file ='create')
    
    cut_refact_finger = [int(l) for l in fingerprint.rstrip().split() if l != split_separator]
    
    if cut_read_index == 1:
        whole_read_id = read_ids[0]
    else:
        whole_read_id = read_ids[1]
    
    whole_finger = [int(l) for l in orig_finger_dict[whole_read_id] if l != split_separator]

    if noerr_reads:
        found = False
        
        j = 0
        while j < min(len(cut_refact_finger),len(whole_finger[start_index:])) and cut_refact_finger[j] == whole_finger[start_index+j]:
            j = j + 1
        
        pref_length = j

        if pref_length >= k and sum(cut_refact_finger[:pref_length]) >= min_total_length:
            found = True

            if cut_read_index == 0:
                first_start_matching = cut_prefix
                first_end_matching = first_start_matching + sum(cut_refact_finger[:pref_length])
                
                second_start_matching = sum(whole_finger[:start_index])
                second_end_matching = second_start_matching + sum(whole_finger[start_index:start_index+pref_length])
            else:
                second_start_matching = cut_prefix
                second_end_matching = second_start_matching + sum(cut_refact_finger[:pref_length])
            
                first_start_matching = sum(whole_finger[:start_index])
                first_end_matching = first_start_matching + sum(whole_finger[start_index:start_index+pref_length])

            current_length = first_end_matching - first_start_matching + 1

        if found == False or current_length < (int(fields[4]) - int(fields[3]) + 1):
            first_start_matching = int(fields[3])
            first_end_matching = int(fields[4]) + 1
            second_start_matching = int(fields[9])
            second_end_matching = int(fields[10]) + 1
            found = True
    else:
        k_finger_dict = dict()
        #Fill the dictionary with the k-fingers of cut_refact_finger
        for i in range(len(cut_refact_finger))[:len(cut_refact_finger)-k][::-1]:
            if sum(tuple(cut_refact_finger[i:i+k])) >= min_total_length:
                k_finger_dict[tuple(cut_refact_finger[i:i+k])] = i

        j = -1
        found = False
        i = start_index
        while i < len(whole_finger) and found == False:
            k_finger = tuple(whole_finger[i:i+k])
            if k_finger in k_finger_dict:
                j = k_finger_dict[k_finger]
                found = True
            else:
                i = i + 1

        if found:
            return_list = get_matching_prefixes(whole_finger[i:], cut_refact_finger[j:], k, downstream_k, min_downstream_length, search_gap, gap_length_tolerance, False)

            if cut_read_index == 0:
                first_start_matching = cut_prefix + sum(cut_refact_finger[:j])
                first_end_matching = first_start_matching + sum(cut_refact_finger[j:j+return_list[1]])

                second_start_matching = sum(whole_finger[:i])
                second_end_matching = second_start_matching + sum(whole_finger[i:i+return_list[0]])
            else:
                second_start_matching = cut_prefix + sum(cut_refact_finger[:j])
                second_end_matching = second_start_matching + sum(cut_refact_finger[j:j+return_list[1]])
        
                first_start_matching = sum(whole_finger[:i])
                first_end_matching = first_start_matching + sum(whole_finger[i:i+return_list[0]])

            current_length = first_end_matching - first_start_matching + 1

        if found == False or current_length < (int(fields[4]) - int(fields[3]) + 1):
            first_start_matching = int(fields[3])
            first_end_matching = int(fields[4]) + 1
            second_start_matching = int(fields[9])
            second_end_matching = int(fields[10]) + 1
            found = True
    
    overlap_ok = False
    
    if found:
        upstream_length = min(first_start_matching, second_start_matching)
        downstream_length = min(int(first_read_length)-first_end_matching, int(second_read_length)-second_end_matching)
        matching_length = first_end_matching - first_start_matching + 1
        full_overlap_length = upstream_length + matching_length + downstream_length

        if matching_length / full_overlap_length >= min_coverage:
                overlap_ok = True

        if overlap_ok and extend_suffix_prefix_overlap:
            first_start_matching = first_start_matching - upstream_length
            second_start_matching = second_start_matching - upstream_length
            first_end_matching = first_end_matching + downstream_length
            second_end_matching = second_end_matching + downstream_length

    if overlap_ok:
        if first_end_matching - first_start_matching + 1 < min_total_length or second_end_matching - second_start_matching + 1 < min_total_length:
            sys.stderr.write("Error!\n")
            exit()
        
        out_string = ""

        out_string = out_string + fields[0] + "\t" + str(first_read_index) + "\t" + str(first_read_length)
    
        out_string = out_string + "\t" + str(first_start_matching) + "\t" + str(first_end_matching-1)
    
        out_string = out_string + "\t0"
    
        out_string = out_string + "\t" + fields[6] + "\t" + str(second_read_index) + "\t" + str(second_read_length)
    
        out_string = out_string + "\t" + str(second_start_matching) + "\t" + str(second_end_matching-1)
    
        out_string = out_string + "\t0"
        out_string = out_string + "\t" + sharing

        if print_error:
            first_read_sequence = "".join(f for f in orig_fact_dict[read_ids[0]] if f != "|")
            second_read_sequence = "".join(f for f in orig_fact_dict[read_ids[1]] if f != "|")
        
            first_overlap = first_read_sequence[first_start_matching:first_end_matching]
            second_overlap = second_read_sequence[second_start_matching:second_end_matching]
       
            error = 0.0
            if first_overlap != second_overlap:
                error = Levenshtein.distance(first_overlap, second_overlap)/max(len(first_overlap), len(second_overlap))*100

            out_string = out_string + "\t" + str(error)

        output_file.write(out_string + "\n")

output_file.close()





