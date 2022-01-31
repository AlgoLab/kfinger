import bz2
import logging
import gzip
import pickle
import sys

MIN_COV = 0.5

logging.basicConfig(
    format='%(levelname)s:%(message)s',
    #level=logging.DEBUG,
    level=logging.INFO
)

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

def extract_region(cg, start, stop, offsetQ, offsetR):
    logging.debug("extract_region start=%d, stop=%d, offsetQ=%d, offsetR=%d", start, stop, offsetQ, offsetR)
    sstart = None
    sstop = None
    cum_query = offsetQ
    cum_ref = offsetR
    num = 0
    c = None
    cgops = cigar_ops(cg)
    if start < cum_query:
        sstart = cum_ref - (cum_query - start)
    else:
        # Find start
        try:
            while cum_query < start:
                (num, c) = next(cgops)
                assert c in ["M", "D", "I"]
                if c == "M" or c == "D":
                    cum_ref += num
                if c == "M" or c == "I":
                    cum_query += num
        except StopIteration:
            pass
        if cum_query == start or (cum_query > start and c != "M"):
            sstart = cum_ref
        elif cum_query > start:
            sstart = cum_ref - (cum_query - start)
    # Find stop
    try:
        while cum_query < stop:
            (num, c) = next(cgops)
            assert c in ["M", "D", "I"]
            if c == "M" or c == "D":
                cum_ref += num
            if c == "M" or c == "I":
                cum_query += num
    except StopIteration:
        pass
    if cum_query == stop or (cum_query > stop and c != "M"):
        sstop = cum_ref
    elif cum_query > stop:
        sstop = cum_ref - (cum_query - stop)
    elif cum_query < stop:
        sstop = cum_ref + (stop - cum_query)
    logging.debug("  -> %d--%d", sstart, sstop)
    return (sstart, sstop)

def ovlp_rate(i1, i2):
    return max(0, (min(i1[1], i2[1]) - max(i1[0], i2[0])) / (max(i1[1], i2[1]) - min(i1[0], i2[0])))

if len(sys.argv) < 3:
    logging.error("Missing arguments. USAGE: %s <index.bin> <ovlp.paf>", sys.argv[0])
    sys.exit(1)

index_fname = sys.argv[1]
ovlp_fname = sys.argv[2]

logging.info("Reading index from file %s", index_fname)
with bz2.open(index_fname, "rb") as f:
    in_dict = pickle.load(f)

MIN_FRAC_MAPPED = in_dict['MIN_FRAC_MAPPED']
MIN_OVLP_RATE = in_dict['MIN_OVLP_RATE']
FRAC_BAD_OVERLAPS = in_dict['FRAC_BAD_OVERLAPS']
reads = in_dict['reads']
to_ref = in_dict['to_ref']
overlaps = in_dict['overlaps']
del in_dict

logging.info("Reading overlaps from file %s", ovlp_fname)
found_overlaps = set()
best_perf = {}
with opener(ovlp_fname) as f:
    for line in f:
        line = line.rstrip().split()
        orig_line0 = line[0]
        orig_line5 = line[5]
        line[0] = line[0][:-2]
        line[5] = line[5][:-2]
        line[4] = "+" if orig_line0[-2:] == orig_line5[-2:] else "-"
        r1 = line[0]
        r2 = line[5]
        orig_r1 = orig_line0
        orig_r2 = orig_line5
        key = (min(r1, r2), max(r1, r2))
        if key in overlaps:
            (line[1], line[2], line[3]) = (int(line[1]), int(line[2]), int(line[3]))
            (line[6], line[7], line[8]) = (int(line[6]), int(line[7]), int(line[8]))
            (len1, start1, stop1, strand) = to_ref[r1][0][1:5]
            (start2, stop2) = to_ref[r1][0][7:9]
            (oth_len1, oth_start1, oth_stop1, oth_strand) = to_ref[r2][0][1:5]
            (oth_start2, oth_stop2) = to_ref[r2][0][7:9]

            true_ovlp_onG = (max(start2, oth_start2), min(stop2, oth_stop2))

            reg = extract_region(to_ref[r1][1], line[2], line[3], start1, start2)
            oth_reg = extract_region(to_ref[r2][1], line[7], line[8], oth_start1, oth_start2)

            pred_ovlp_onG = (max(reg[0], oth_reg[0]), min(reg[1], oth_reg[1]))

            ovlr = ovlp_rate(true_ovlp_onG, pred_ovlp_onG)

            if ((strand == oth_strand) != (line[4] == "+") or ovlr < MIN_COV):
                logging.info("Overlap %s-%s (%.3f)", orig_r1, orig_r2, ovlr)
                logging.error("OVERLAP IS NOT EQUAL: %d-%d   %d-%d", line[2], line[3], line[7], line[8])
                logging.info(to_ref[r1])
                logging.info(to_ref[r2])
                logging.info(line)
                logging.info(overlaps[key])
                logging.info("On genomic: True: %d-%d, Pred: %d-%d", *true_ovlp_onG, *pred_ovlp_onG)
            else:
                found_overlaps.add(key)
            if key not in best_perf or best_perf[key][0] < ovlr:
                best_perf[key] = (ovlr, (strand == oth_strand) == (line[4] == "+"), *pred_ovlp_onG)
        else:
            #logging.info("not found")
            pass

logging.info("Found overlaps: %d", len(found_overlaps))
logging.info("Not found overlaps: %d", len(overlaps.keys() - found_overlaps))

print("r1,r2,true_start_on_G,true_stop_on_G,true_ovlp_edit_dist,ovlp_ratio,strand_ok,pred_start_on_G,pred_stop_on_G")
for r1,r2 in overlaps.keys():
    (start2, stop2) = to_ref[r1][0][7:9]
    (oth_start2, oth_stop2) = to_ref[r2][0][7:9]
    ovlp_start2 = max(start2, oth_start2)
    ovlp_stop2 = min(stop2, oth_stop2)
    print(
        "\"{}\",\"{}\",{},{},{},{},{},{},{}".format(
            r1, r2, ovlp_start2, ovlp_stop2, overlaps[(r1, r2)]["editDistance"],
            *(best_perf[(r1,r2)] if (r1,r2) in best_perf else (-1,-1,-1,-1))
        )
    )

        



