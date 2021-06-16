# Controlla gli errori dei matching di kfinger

import sys

overlap_file = sys.argv[1]

# Output path
paf_file = sys.argv[2]

# Le coordinate sono 0-based
# Gli indici dei reads sono 1-based

with open(overlap_file,'r') as overlap_input:
    overlap_list = overlap_input.readlines()

paf_input = open(paf_file,'w')

parameter_records = overlap_list[:overlap_list.index("//\n")]

#for par in parameter_records:
    #paf_input.write(par)

#paf_input.write('# CONVERSION PARAMETERS:'  + "\n")
##paf_input.write('# Script: ' + sys.argv[0]  + "\n")
#paf_input.write('# Overlap file: ' + str(overlap_file)  + "\n")
#paf_input.write('# PAF file: ' + str(paf_file)  + "\n")
#paf_input.write("//\n")

overlap_list = overlap_list[overlap_list.index("//\n")+1:]

for record in overlap_list:
    
    record = record.rstrip()
    
    sys.stderr.write(record + "\n")
    
    record = record.split()
    
    if record[1] != record[7]:
        id1 = record[0]
        length1 = record[2]
        
        start1 = record[3]
        end1 = int(record[4])+1
    
        id2 = record[6]
        length2 = record[8]
    
        start2 = record[9]
        end2 = int(record[10])+1
        
        strand = '+'
 
        if record[5] != record[11]:
            strand = '-'

        new_record = [id1, length1, start1, str(end1), strand, id2, length2, start2, str(end2), length1, length1, '255']
        paf_input.write("\t".join(new_record) + "\n")
            



