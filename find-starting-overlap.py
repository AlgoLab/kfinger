#!/usr/bin/env python
# coding: utf-8

# ### Start overlap detection (only for single-stranded factorizations)

# Double-stranded factorizations are not handled

import time

import sys
import Levenshtein
import itertools

from collections import Counter,defaultdict

# PRIMARY PARAMETERS

# Input file of fingerprints
input_file_name = sys.argv[1]

# Input of file factorizations
input_file_name_f = sys.argv[2]

# k-fingers size
k = int(sys.argv[3]) # 7

# Minimum total length for retaining a k-finger
min_total_length = int(sys.argv[4])  # 40

# Minimum number of common unique kfingers between two reads to be considered candidate overlapping reads
min_number_unique_kfingers = int(sys.argv[5])   # 6

# Minimum length for retaining an input read
min_read_length = int(sys.argv[6])   # 0

# Maximum length for retaining an input read (0, if any length must be retained)
max_read_length = int(sys.argv[7])   # 0

min_coverage = float(sys.argv[8])   # 0.0

# 1: the suffix-prefix overlaps are produced in output
# 0: the matchings are produced in output
extend_suffix_prefix_overlap = bool(int(sys.argv[9]))   #  0

# The following four parameters are active if search_gap is True
downstream_k = int(sys.argv[10])       # 3 for error-free datasets, 2 otherwise
min_downstream_length = int(sys.argv[11])  # 10 for error-free datasets, 20 otherwise
gap_length_tolerance = int(sys.argv[12])   # 0 for error-free datasets, 15 otherwise

# Path of the output directory
# output_path = sys.argv[13]
output_file_name = sys.argv[13]

# SECONDARY PARAMETERS

# Consider k-fingers having at most max_number_of_one 1-long factors (-1 if a k-finger can contain any number of 1s)
max_number_of_one = -1   # -1

# False: the matching substring is limited to the leftmost k-finger extended on the right as long as possible
# True: the matching substring is obtained by serching the rightmost gap after the leftmost k-finger conserving the length on the two reads
search_gap = True   # True

multiple_gaps = False       # False

#Separator for the input fingerprints|factorizations
split_separator = '|'   # '|'

# True: the error of the matching substring|suffix-prefix overlap (depending on print_suffix_prefix_overlap) is produced in output
print_error = False     # False

# Maximum number of overlap produced in output
only_first_overlaps = -1     # -1, all the overlaps are produced as output

# Also the FASTA file of the cut reads is produced
# print_read = False   # False

print_log = False   # False

# output_file_name = output_path + "raw-matchings.txt"
output_file = open(output_file_name, "w")

# if print_read:
#     read_file_name = output_path + "cut-reads.fa"
#     read_file = open(read_file_name, "w")

output_file.write('# RAW MATCHING PARAMETERS:\n')
output_file.write('# Script: ' + sys.argv[0] + "\n")
output_file.write('# Fingerprint file: ' + str(input_file_name) + "\n")
output_file.write('# Factorization file: ' + str(input_file_name_f) + "\n")
output_file.write('# k: ' + str(k) + "\n")
output_file.write('# Minimum total kfinger length: ' + str(min_total_length) + "\n")
output_file.write('# Minimum number of unique shared kfingers: ' + str(min_number_unique_kfingers) + "\n")
output_file.write('# Minimum length for considering a read: ' + str(min_read_length) + "\n")
output_file.write('# Maximum length for considering a read: ' + str(max_read_length) + "\n")
output_file.write('# Maximum number of 1s for considering a k-finger (-1: no checking): ' + str(max_number_of_one) + "\n")
output_file.write('# Minimum matching rate: ' + str(min_coverage) + "\n")
output_file.write('# Extend to suffix-prefix overlap: ' + str(extend_suffix_prefix_overlap) + "\n")
output_file.write('# Search gap: ' + str(search_gap) + "\n")
if search_gap:
    output_file.write('\t# Minimum number of factors required after the gap: ' + str(downstream_k) + "\n")
    output_file.write('\t# Minimum total length of the factors required after the gap: ' + str(min_downstream_length) + "\n")
    output_file.write('\t# Maximum tolerance in the gap length: ' + str(gap_length_tolerance) + "\n")
    output_file.write('\t# Multiple gaps: ' + str(multiple_gaps) + "\n")
if print_error:
    output_file.write('# Error in the last field\n')

if only_first_overlaps != -1:
    output_file.write('# Only the first ' + str(only_first_overlaps) + ' overlaps are produced\n')

output_file.write('# Output file: ' + str(output_file_name) + "\n")


def get_matching_prefixes(f1, f2, k, downstream_k, min_downstream_length, search_gap, gap_length_tolerance, multiple_gaps):

    lenf1 = len(f1)
    lenf2 = len(f2)
    #Find the longest common prefix of f1 and f2 (supposing that the first k integers are already identical)
    init_common_pos = k
    while init_common_pos < lenf1 and init_common_pos < lenf2 and f1[init_common_pos] == f2[init_common_pos]:
        init_common_pos += 1
    
    #Variable init_common_pos contains the 0-based position after the end of the common prefix
    
    return_list = [init_common_pos, init_common_pos]

    if search_gap == False:
        return return_list

    # SEARCH THE LONGEST GAP

    #Fill the dictionary with the k-fingers of f1 for k = downstream_k
    sumf1 = tuple(itertools.accumulate(f1))
    common_right_dict = {
        tuple(f1[i:i+downstream_k]):i
        for i in range(init_common_pos, lenf1-downstream_k)
        if sumf1[i+downstream_k-1] - sumf1[i-1] >= min_downstream_length
    }

    # The dictionary contains all the k-fingers of f1 for k = downstream_k (as keys) together their starting positions as values

    # Scan f2 from the end to position init_common_pos in order to search the common k-fingers for k = downstream_k
    i = lenf2-downstream_k
    found = False
    sumf2 = tuple(itertools.accumulate(f2))
    while i >= init_common_pos:
        key = tuple(f2[i:i+downstream_k])
        if  key in common_right_dict and abs(sumf1[common_right_dict[key]-1] - sumf2[i-1]) <= gap_length_tolerance:
            found = True
            break
        else:
            i -= 1

    if found:
        #Find the longest common factor of f1 and f2 after the gap
        i = i + downstream_k
        p = common_right_dict[key] + downstream_k
        while (p < lenf1 and i < lenf2) and (f1[p] == f2[i]):
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
            
                while q >= init_common_pos:
                    key = tuple(f2[j:q+downstream_k])
                    if  key in common_right_dict and abs(sum(f1[:common_right_dict[key]]) - sum(f2[:q])) <= gap_length_tolerance:
                        k = common_right_dict[key]
                        return_list.append([k,q])
                        break
                    
                    q = q - 1
                       
                if k <= init_common_pos or q <= init_common_pos:
                    stop = True
                
    return return_list

with open(input_file_name_f,'r') as input_file_f:
    file_input_f = input_file_f.readlines()

with open(input_file_name,'r') as input_file:
    file_input = tuple(read.split() for read in input_file)

# Dizionario id_dict --> KEY: indice 0-based del read nel file di input; VALUE: lista di dimensione 1 contenente l'ID del read nel file di input
id_dict = [ [n[0], -1] for n in file_input ]

# Dizionario finger_dict --> KEY: k-finger (tupla di interi); VALUE: lista degli indici 0-based dei reads che contengono la k-finger chiave (la lista NON è un multi-insieme in quanto le k-fingers che occorrerono più volte in un read vengono rimosse)

# Ai valori (liste di dimensione 1) del dizionario id_dict viene a questo punto aggiunto il numero di k-fingers che il read (chiave) contiene (NOTA: il valore non viene al momento usato); le liste dopo il for seguente avranno dimensione 2

finger_dict = defaultdict(list)

count_retained_by_length = 0

max_length = 0

# Dizionario leftmost_dict: KEY: indice 0-based di un read i; VALUE: dizionario interno
# dizionario interno --> KEY: indice 0-based di un read j; VALUE: lista [off_i, off_j] dei due offset (in termini di nucleotidi) delle leftmost kfingers
# in comune tra i e j
# NOTA: si ha sempre i < j
leftmost_dict = {}

# Dizionario leftmost_offset_dict: KEY: k-finger (tupla di interi); VALUE: lista degli offset delle k-fingers chiave (corrispondente ai reads nel dizionario finger_dict)
leftmost_offset_dict = defaultdict(list)
# Dizionario leftmost_index_dict: KEY: k-finger (tupla di interi); VALUE: lista degli indici di inizio delle k-fingers chiave (corrispondente ai reads nel dizionario finger_dict)
leftmost_index_dict = defaultdict(list)

split_length = sum(int(l) for l in " ".join(file_input[0][1:]).split(split_separator)[0].split())

for i in range(len(file_input)):
    
    sys.stderr.write(str(i)+'\n')
    
    read = file_input[i]

    id = read[0]
    finger = read[1:]
    
    finger = [int(f) for f in finger if f != split_separator]
    read_length = sum(finger)
    
    if read_length > max_length:
        max_length = read_length
    
    if read_length >= min_read_length and (max_read_length == 0 or read_length <= max_read_length):

        id_dict[i][1] = len(finger)-k+1
        
        temp_list_kfinger = []

        for j in range(len(finger)-k+1):
            kfinger_tuple = tuple(finger[j:j+k])
            ok_one = True
            if max_number_of_one != -1 and kfinger_tuple.count(1) > max_number_of_one:
                ok_one = False
            
            if ok_one and sum(kfinger_tuple) >= min_total_length:
                finger_dict[kfinger_tuple].append(i)
                leftmost_offset_dict[kfinger_tuple].append(sum(finger[:j]))
                leftmost_index_dict[kfinger_tuple].append(j)

                temp_list_kfinger.append(kfinger_tuple)

        c = Counter(temp_list_kfinger)
        for key in c:
            if c[key] > 1:
                value = finger_dict[key]
                value[len(value)-c[key]:len(value)] = []
                finger_dict[key] = value
                
                value_off = leftmost_offset_dict[key]
                value_off[len(value_off)-c[key]:len(value_off)] = []
                leftmost_offset_dict[key] = value_off

                value_index = leftmost_index_dict[key]
                value_index[len(value_index)-c[key]:len(value_index)] = []
                leftmost_index_dict[key] = value_index

        for p in range(len(temp_list_kfinger)):
            key = temp_list_kfinger[p]
            
            if c[key] == 1:
                value = finger_dict[key]
                value_off = leftmost_offset_dict[key]
                value_index = leftmost_index_dict[key]

                for j in range(len(value)-1):
                    read = value[j]
                    offset = value_off[j]
                    index = value_index[j]
                    d = leftmost_dict.get(read, dict())
                    if i not in d:
                        d[i] = [offset, value_off[len(value)-1], index, value_index[len(value)-1]]
                        leftmost_dict[read] = d

        count_retained_by_length += 1

del leftmost_offset_dict
del leftmost_index_dict

output_file.write('# Number of retained reads: ' + str(int(count_retained_by_length/2)) + ' out of ' + str(int(len(file_input)/2)) + "\n")
output_file.write('//' + "\n")

# Dizionario sharing_dict --> KEY: tupla di due indici 0-based di reads; VALUE: numero di k-fingers che i due reads chiave hanno in comune (vengono contate solo le k-fingers che sono uniche nei due reads; il primo indice di read è sempre minore del secondo)

sharing_dict = defaultdict(int)

# Il blocco seguente rimuove i duplicati anche se di fatto non ci sono perché filtrati sopra (anzi è proprio inutile usare
for kfinger in finger_dict:
    
    #sys.stderr.write(str(kfinger)+'\n')

    shared_list = finger_dict[kfinger]
    
    for i in range(len(shared_list)):
        for j in range(i+1, len(shared_list)):
            if shared_list[i] != shared_list[j]:
                sharing_dict[(shared_list[i], shared_list[j])] += 1

# Dizionario filtered_sharing_dict --> KEY: tupla di due indici 0-based di reads; VALUE: numero di k-fingers che i due reads chiave hanno in comune (non compaiono le coppie di reads che condividono un numero di k-fingers uniche al di sotto della soglia minima min_number_unique_kfingers)

#sys.stderr.write(str('Filtering...\n'))

filtered_sharing_dict = {
    pair:numb
    for pair,numb in sharing_dict.items()
    if numb >= min_number_unique_kfingers
}
del sharing_dict

ordered_keys = sorted(filtered_sharing_dict.keys())

read_dict = {}

count = 0

# Le coppie di reads (candidate a essere in overlap, cioè che condividono una soglia minima di k-fingers uniche) vengono considerate in ordine crescente
for pair in ordered_keys:
    # Indice del primo read della coppia (quello minore)
    first_read_index = pair[0]
    
    # Recupero tutte le informazioni del primo read della coppia (quello di indice minore).

    first_read_f = file_input_f[first_read_index]
    first_read_list_f = first_read_f.split()[1:]
    first_read_finger_f = [factor for factor in first_read_list_f[:] if factor != split_separator]
    
    first_read_finger = [int(l) for l in file_input[first_read_index][1:] if l != '|']
   
    # Sequenza del read di indice minore
    first_read_sequence = ''.join(first_read_finger_f)

    # Indice del secondo read della coppia (quello maggiore)
    second_read_index = pair[1]
    
    sys.stderr.write('Reads ' + str(first_read_index) + ' ' + str(id_dict[first_read_index][0]) + ' ' + str(second_read_index) + ' ' + str(id_dict[second_read_index][0]) + '\n')
    
    # Recupero tutte le informazioni del secondo read della coppia (quello di indice maggiore).s

    second_read_list_f = file_input_f[second_read_index].split()[1:]
    second_read_finger_f = [factor for factor in second_read_list_f[:] if factor != '|']

    second_read_finger = [int(l) for l in file_input[second_read_index][1:] if l != '|']

    # Sequenza del read di indice maggiore
    second_read_sequence = ''.join(second_read_finger_f)
    
    first_downstream_fingerprint = first_read_finger[leftmost_dict[first_read_index][second_read_index][2]:]
    second_downstream_fingerprint = second_read_finger[leftmost_dict[first_read_index][second_read_index][3]:]
    
    return_list = get_matching_prefixes(first_downstream_fingerprint, second_downstream_fingerprint, k, downstream_k, min_downstream_length, search_gap, gap_length_tolerance, multiple_gaps)
    
    overlap_length1 = sum(first_downstream_fingerprint[:return_list[0]])
    overlap_length2 = sum(second_downstream_fingerprint[:return_list[1]])
    
    left_most_first_read_offset = leftmost_dict[first_read_index][second_read_index][0]
    left_most_second_read_offset = leftmost_dict[first_read_index][second_read_index][1]
    
    first_left_substr = left_most_first_read_offset
    second_left_substr = left_most_second_read_offset
    first_right_substr = left_most_first_read_offset + overlap_length1
    second_right_substr = left_most_second_read_offset + overlap_length2
    
    upstream_length = min(first_left_substr, second_left_substr)
    downstream_length = min(len(first_read_sequence)-first_right_substr, len(second_read_sequence)-second_right_substr) - 1
    
    overlap_ok = True
    
    overlap_rate = max(overlap_length1, overlap_length2) / (max(overlap_length1, overlap_length2) + upstream_length + downstream_length)
    if overlap_rate < min_coverage:
        overlap_ok = False
    
    if overlap_ok:
        #Extension to the suffix-prefix overlap
        if extend_suffix_prefix_overlap:
            first_left_substr = first_left_substr - upstream_length
            second_left_substr = second_left_substr - upstream_length
            first_right_substr = first_right_substr + downstream_length
            second_right_substr = second_right_substr + downstream_length

        first_overlap = first_read_sequence[first_left_substr:first_right_substr]
        second_overlap = second_read_sequence[second_left_substr:second_right_substr]
    
        first_start_prefix_split = left_most_first_read_offset % split_length
        second_start_prefix_split = left_most_second_read_offset % split_length
    
        start_index = -1
    
        if left_most_first_read_offset < left_most_second_read_offset:
            which_read = 1
        
            if left_most_first_read_offset < second_start_prefix_split:
                start_cut = left_most_first_read_offset
                upstream_split = left_most_second_read_offset
                start_index = leftmost_dict[first_read_index][second_read_index][3]
            else:
                start_cut = (left_most_first_read_offset - second_start_prefix_split) % split_length
                split_count = (left_most_first_read_offset - second_start_prefix_split - start_cut) // split_length
                upstream_split = left_most_second_read_offset - second_start_prefix_split - split_count * split_length
                window_list = " ".join(file_input[second_read_index]).split(split_separator)

            #upstream_split = left_most_second_read_offset - second_start_prefix_split - split_count * split_length
            #initial_start_cut_seq = first_read_sequence[start_cut:start_cut+split_length]
            initial_start_cut_seq = first_read_sequence[start_cut:]
            #window_cut_seq = second_read_sequence[upstream_split:upstream_split+split_length]
            window_cut_seq = second_read_sequence[upstream_split:]
        else:
            which_read = 2
            if left_most_second_read_offset < first_start_prefix_split:
                start_cut = left_most_second_read_offset
                upstream_split = left_most_first_read_offset
                start_index = leftmost_dict[first_read_index][second_read_index][2]
            else:
                start_cut = (left_most_second_read_offset - first_start_prefix_split) % split_length
                split_count = (left_most_second_read_offset - first_start_prefix_split - start_cut) // split_length
                upstream_split = left_most_first_read_offset - first_start_prefix_split - split_count * split_length
                window_list = " ".join(file_input[first_read_index]).split(split_separator)

            #initial_start_cut_seq = second_read_sequence[start_cut:start_cut+split_length]
            initial_start_cut_seq = second_read_sequence[start_cut:]
            #window_cut_seq = first_read_sequence[upstream_split:upstream_split+split_length]
            window_cut_seq = first_read_sequence[upstream_split:]

        if start_index == -1:
            window_num = upstream_split // split_length
            for i in range(window_num):
                start_index = start_index + len(window_list[i].split())

        # Caso in cui la sottostringa di match ha lo stesso start ed end sui due reads (cioé di fatto un read è prefisso dell'altro). E start è diverso da 0 solo per i valori scelti per i parametri
        if start_index == -1 and start_cut == 0:
            start_index = 0

        out_string = id_dict[first_read_index][0] + "\t" + str(first_read_index+1) + "\t" + str(len(first_read_sequence))

        out_string = out_string + "\t" + str(first_left_substr) + "\t" + str(first_right_substr-1)
            
        out_string = out_string + "\t0"
            
        out_string = out_string + "\t" + id_dict[second_read_index][0] + "\t" + str(second_read_index+1) + "\t" + str(len(second_read_sequence))

        out_string = out_string + "\t" + str(second_left_substr) + "\t" + str(second_right_substr-1)

        out_string = out_string + "\t0"
        out_string = out_string + "\t" + str(filtered_sharing_dict[pair])

        out_string = out_string + "\t" + str(which_read) + "\t" + str(start_cut) + "\t" + str(start_index)

        if print_error:
            error = 0.0
            if first_overlap != second_overlap:
                error = Levenshtein.distance(first_overlap, second_overlap)/max(len(first_overlap), len(second_overlap))*100
       
            out_string = out_string + "\t" + str(error)

        output_file.write(out_string + "\n")

#         if print_read:
#             if which_read == 1:
#                 read_file.write(">" + id_dict[first_read_index][0] + "_" + str(start_cut) + "\n")
#                 read_file.write(first_read_sequence[start_cut:] + "\n")
#             else:
#                 read_file.write(">" + id_dict[second_read_index][0] + "_" + str(start_cut) + "\n")
#                 read_file.write(second_read_sequence[start_cut:] + "\n")

        count = count + 1

    if overlap_ok and print_log:
        
        print('Reads ' + str(first_read_index) + ' ' + str(id_dict[first_read_index][0]) + ' ' + str(second_read_index) + ' ' + str(id_dict[second_read_index][0]))
        
        read_dict[id_dict[first_read_index][0]] = first_read_sequence
        read_dict[id_dict[second_read_index][0]] = second_read_sequence

        print("Matching substrings:")
        print(first_overlap)
        print(second_overlap)
        
        print("Read " + str(which_read) + " must be cut from position " + str(start_cut))
        print(initial_start_cut_seq)
        print("Starting window of the other read:")
        print(window_cut_seq)
        if which_read == 1:
            print(" ".join(str(l) for l in second_read_finger_f[start_index:]))
            print(" ".join(str(l) for l in second_read_finger[start_index:]))
        else:
            print(" ".join(str(l) for l in first_read_finger_f[start_index:]))
            print(" ".join(str(l) for l in first_read_finger[start_index:]))

        check_length = split_length - first_start_prefix_split
        print_first_read_finger_with_sep = []
        print_first_read_factor_with_sep = []
        len_sum = 0
        f_list = first_read_finger[leftmost_dict[first_read_index][second_read_index][2]:]
        factor_list = first_read_finger_f[leftmost_dict[first_read_index][second_read_index][2]:]
        for i in range(len(f_list)):
            print_first_read_finger_with_sep.append(f_list[i])
            print_first_read_factor_with_sep.append(factor_list[i])
            len_sum = len_sum + f_list[i]
            if len_sum == check_length:
                print_first_read_finger_with_sep.append("|")
                print_first_read_factor_with_sep.append("|")
                check_length = split_length
                len_sum = 0
            else:
                if len_sum > check_length:
                    sys.stderr.write('Error!\n')
                    exit()

        check_length = split_length - second_start_prefix_split
        print_second_read_finger_with_sep = []
        print_second_read_factor_with_sep = []
        len_sum = 0
        f_list = second_read_finger[leftmost_dict[first_read_index][second_read_index][3]:]
        factor_list = second_read_finger_f[leftmost_dict[first_read_index][second_read_index][3]:]
        for i in range(len(f_list)):
            print_second_read_finger_with_sep.append(f_list[i])
            print_second_read_factor_with_sep.append(factor_list[i])
            len_sum = len_sum + f_list[i]
            if len_sum == check_length:
                print_second_read_finger_with_sep.append("|")
                print_second_read_factor_with_sep.append("|")
                check_length = split_length
                len_sum = 0
            else:
                if len_sum > check_length:
                    sys.stderr.write('Error!\n')
                    exit()

        #Portions of the fingerprints after and before the leftmost common kfingers with separator bar

        print("Fingerprints in the matching substrings:")
        print("\t".join(str(i) for i in first_downstream_fingerprint[:return_list[0]]))
        print("\t".join(str(i) for i in second_downstream_fingerprint[:return_list[1]]))
        
        print("Fingerprints from the leftmost k-finger:")
        print("\t".join(str(i) for i in print_first_read_finger_with_sep))
        print("\t".join(str(i) for i in print_second_read_finger_with_sep))
        print("\t".join(str(f) for f in print_first_read_factor_with_sep))
        print("\t".join(str(f) for f in print_second_read_factor_with_sep))
        
        if multiple_gaps:
            print("Other gaps:")
            for gap in return_list[2:]:
                print("Gap:")
                print("\t".join(str(i) for i in first_downstream_fingerprint[:gap[0]]))
                print("\t".join(str(i) for i in second_downstream_fingerprint[:gap[1]]))
                
        print("//")

    if count == only_first_overlaps:
        exit()

output_file.close()

# if print_read:
#     read_file.close()
