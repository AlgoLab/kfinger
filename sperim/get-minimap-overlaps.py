import Levenshtein
import sys

minimap_file = sys.argv[1]
min_coverage = float(sys.argv[2])   # 0.80
output_path = sys.argv[3]

with open(minimap_file,'r') as minimap_input:
    minimap_list = minimap_input.readlines()

output_file_name = output_path + "full-minimap-overlaps.paf"
output_file = open(output_file_name, "w")

for record in minimap_list:
    
    record = record.rstrip()
    
    sys.stderr.write(record + "\n")
    
    record = record.split()
    
    start1 = int(record[2])
    end1 = int(record[3])
    
    start2 = int(record[7])
    end2 = int(record[8])
    
    length1 = int(record[1])
    length2 = int(record[6])
    
    upstream_length = min(start1, start2)
    downstream_length = min(length1-end1, length2-end2)
    
    overlap_rate = max(end1-start1, end2-start2) / (max(end1-start1, end2-start2) + upstream_length + downstream_length)
    
    if overlap_rate >= min_coverage:
        record[2] = str(start1 - upstream_length)
        record[3] = str(end1 + downstream_length)
        
        record[7] = str(start2 - upstream_length)
        record[8] = str(end2 + downstream_length)

        output_file.write("\t".join(record) + "\n")





