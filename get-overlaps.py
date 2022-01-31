# Ottiene gli overlap prefisso-suffisso dai matchings (considerati solo single-stranded)

import sys

matching_file = sys.argv[1]

min_coverage = float(sys.argv[2])   # 0.80

output_file_name = sys.argv[3]

# Le coordinate sono 0-based
# Gli indici dei reads sono 1-based

with open(matching_file,'r') as matching_input:
    matching_list = matching_input.readlines()

parameter_records = matching_list[:matching_list.index("//\n")]

output_file = open(output_file_name, "w")

for par in parameter_records:
    output_file.write(par)

output_file.write("# GET FULL OVERLAPS:" + "\n")
output_file.write("# Matching file: " + matching_file + "\n")
output_file.write("# Min coverage: " + str(min_coverage) + "\n")
output_file.write("# Output file: " + output_file_name + "\n")
output_file.write("//" + "\n")

matching_list = matching_list[matching_list.index("//\n")+1:]

for record in matching_list:
    
    record = record.rstrip()
    
    sys.stderr.write(record + "\n")
    
    record = record.split()
    
    length1 = int(record[2])
    
    start1 = int(record[3])
    end1 = int(record[4])
    
    length2 = int(record[8])

    start2 = int(record[9])
    end2 = int(record[10])

    upstream_length = min(start1, start2)
    downstream_length = min(length1-end1, length2-end2) - 1
    
    overlap_rate = max(end1-start1+1, end2-start2+1) / (max(end1-start1+1, end2-start2+1) + upstream_length + downstream_length)
    
    if overlap_rate >= min_coverage:
        record[3] = str(start1 - upstream_length)
        record[4] = str(end1 + downstream_length)
        
        record[9] = str(start2 - upstream_length)
        record[10] = str(end2 + downstream_length)

        output_file.write("\t".join(record) + "\n")


