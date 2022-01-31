import bz2
import logging
import gzip
from math import ceil
import pickle
import sys


import edlib

logging.basicConfig(
    format='%(levelname)s:%(message)s',
    #level=logging.DEBUG,
    level=logging.INFO,
)

MIN_FRAC_MAPPED = 0.95
MIN_OVLP_RATE = 0.80
FRAC_BAD_OVERLAPS = 0.15

rc = str.maketrans("ACGT", "TGCA")
def revcompl(seq):
    return seq.translate(rc)[::-1]

def opener(filename):
    f = open(filename, 'rb')
    if (f.read(2) == b'\x1f\x8b'):
        f.close()
        return gzip.open(filename, 'rt')
    else:
        f.close()
        return open(filename, 'rt')

def cigar_ops(cg):
    num = 0
    for c in cg:
        if c.isdigit():
            num = num*10 + int(c)
        else:
            yield (num, c)
            num = 0

def extract_region(seq, cg, start, stop):
    logging.debug("extract_region len(seq)=%d, start=%d, stop=%d", len(seq), start, stop)
    sstart = None
    sstop = None
    cum_query = 0
    cum_ref = 0
    num = 0
    c = None
    cgops = cigar_ops(cg)
    # Find start
    try:
        while cum_ref < start:
            (num, c) = next(cgops)
            assert c in ["M", "D", "I"]
            if c == "M" or c == "D":
                cum_ref += num
            if c == "M" or c == "I":
                cum_query += num
        if cum_ref == start:
            sstart = cum_query
        if cum_ref > start:
            if c == "M":
                sstart = cum_query - (cum_ref-start)
            else: # "D"
                sstart = cum_query
    except StopIteration:
        pass
    # Find stop
    try:
        while cum_ref < stop:
            (num, c) = next(cgops)
            assert c in ["M", "D", "I"]
            if c == "M" or c == "D":
                cum_ref += num
            if c == "M" or c == "I":
                cum_query += num
        if cum_ref == stop:
            sstop = cum_query
        if cum_ref > stop:
            if c == "M":
                sstop = cum_query - (cum_ref-stop)
            else: # "D"
                sstop = cum_query
    except StopIteration:
        pass
    return seq[sstart:sstop]

if len(sys.argv) < 4:
    logging.error("Missing arguments. USAGE: %s <READS.fa> <to_reaf.paf> <out_file.bin>", sys.argv[0])
    sys.exit(1)

reads_fname = sys.argv[1]
to_ref_fname = sys.argv[2]
out_fname = sys.argv[3]

reads = {}
with opener(reads_fname) as f:
    rid = None
    curr_read = []
    for line in f:
        if line.startswith('>'):
            if rid is not None:
                reads[rid] = "".join(curr_read)
                curr_read = []
            rid = line.lstrip('>').rstrip()
        else:
            curr_read.append(line.rstrip())
    if rid is not None:
        reads[rid] = "".join(curr_read)

to_ref_v = [(0, 0, 0, None)]
to_ref_m = {}
discarded_min_frac_mapped = 0
discarded_min_len = 0
discarded_min_ovlp_len = 0
n_overlaps = 0
n_bad_overlaps = 0
overlaps = {}
with opener(to_ref_fname) as f:
    for line in f:
        line = line.rstrip().split()
        logging.info("READ %s", line[0])
        if not any((x == 'tp:A:P' for x in line[12:])):
            logging.info("Skipped non-primary: %s", line)
            continue
        cg = list(filter(lambda x: x.startswith("cg:Z:"), line[12:]))
        assert len(cg) > 0
        cg = cg[0][5:]
        (len1, start1, stop1) = (line[1], line[2], line[3]) = (int(line[1]), int(line[2]), int(line[3]))
        if stop1 - start1 < MIN_FRAC_MAPPED * len1:
            discarded_min_frac_mapped += 1
            logging.info("Discarded read %s because has %f < %f mapped length", line[0], (stop1 - start1) / len1, MIN_FRAC_MAPPED)
            continue
        (len2, start2, stop2) = (line[6], line[7], line[8]) = (int(line[6]), int(line[7]), int(line[8]))
        #assert line[4] == '+'

        # Start and stop on genomic + hard clipped pref/suf
        start3 = start2 - start1
        stop3 = stop2 + (len1 - stop1)

        pos = len(to_ref_v)
        while pos > 0 and to_ref_v[pos-1][2] > start2:
            pos -= 1
            (oth_start2, oth_stop2, _, oth_id, oth_start1, oth_stop1, oth_len1) = to_ref_v[pos]
            
            # Start and stop on genomic + hard clipped pref/suf
            oth_start3 = oth_start2 - oth_start1
            oth_stop3 = oth_stop2 + (oth_len1 - oth_stop1)
            
            if (min(oth_stop2, stop2) - max(oth_start2, start2)) / (max(stop3, oth_stop3) - min(start3, oth_start3)) < MIN_OVLP_RATE:
                logging.debug("Discarded short overlap %d-%d -- %d-%d", start2, stop2, oth_start2, oth_stop2)
                discarded_min_ovlp_len += 1
                continue
            ovlp_start2 = max(start2, oth_start2)
            ovlp_stop2 = min(stop2, oth_stop2)
            logging.debug(
                "Overlap with read %s [map: %d-%d -- %d-%d] (ovlp: %d--%d)",
                oth_id, start2, stop2, oth_start2, oth_stop2,
                ovlp_start2, ovlp_stop2
            )
            ed = None
            read1 = reads[line[0]][start1:stop1]
            if line[4] == '-':
                read1 = revcompl(read1)
            reg1 = extract_region(read1, cg, ovlp_start2 - start2, ovlp_stop2 - start2)
            read2 = reads[oth_id][oth_start1:oth_stop1]
            if to_ref_m[oth_id][0][4] == '-':
                read2 = revcompl(read2)
            reg2 = extract_region(read2, to_ref_m[oth_id][1], ovlp_start2 - oth_start2, ovlp_stop2 - oth_start2)
            ed = edlib.align(reg1, reg2, k = int(ceil(max(len(reg1), len(reg2))*FRAC_BAD_OVERLAPS)))
            logging.debug(reg1[:110] + "..." + reg1[-30:])
            logging.debug(reg2[:110] + "..." + reg2[-30:])
            logging.debug(ed)
            if ed['editDistance'] == -1:
                n_bad_overlaps += 1
            if line[0] < oth_id:
                overlaps[(line[0], oth_id)] = (ed)
            else:
                overlaps[(oth_id, line[0])] = (ed)

            n_overlaps += 1

        to_ref_v.append((start2, stop2, max(stop2, to_ref_v[-1][2]), line[0], start1, stop1, len1))
        to_ref_m[line[0]] = (line, cg)
        

logging.info("Writing to file %s", out_fname)
out_dict = {
    'MIN_FRAC_MAPPED': MIN_FRAC_MAPPED,
    'MIN_OVLP_RATE': MIN_OVLP_RATE,
    'FRAC_BAD_OVERLAPS': FRAC_BAD_OVERLAPS,
    'reads': reads,
    'to_ref': to_ref_m,
    'overlaps': overlaps,
}
with bz2.open(out_fname, "wb") as f:
    pickle.dump(out_dict, f)


logging.info("No. of reads: %d", len(reads))
logging.info("Not min_frac_mapped (%.4f): %d", MIN_FRAC_MAPPED, discarded_min_frac_mapped)
logging.info("No. of discarded short overlaps (%f): %d", MIN_OVLP_RATE, discarded_min_ovlp_len)
logging.info("No. of overlaps: %d", n_overlaps)
logging.info("No. of bad overlaps (%.4f): %d", FRAC_BAD_OVERLAPS, n_bad_overlaps)