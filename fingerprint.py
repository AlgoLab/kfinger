import argparse

import shutil
import os.path

from multiprocessing.pool import Pool
from functools import partial

from utils.fingerprint_utils import extract_reads,extract_long_reads,compute_fingerprint_by_list,compute_long_fingerprint_by_list,mapping_projection
from utils.factorizations import CFL, ICFL_recursive, CFL_icfl, CFL_icfl_cfl
from utils.factorizations_comb import d_cfl, d_icfl, d_cfl_icfl


# Create fingerprint files (args.step = 'fingerprint') #################################################################
def experiment_fingerprint_1f_np_step(args):

    # Input FASTA file containing transcripts
    input_fasta = args.fasta

    # Extract of reads (Format = ID GENE read)
    read_lines = extract_reads(name_file=input_fasta, filter=args.filter)

    if len(read_lines) == 0:
        print('No reads extracted!')
        exit(-1)

    print('\nCompute fingerprint by list (%s, %s) - start...' % (args.type_factorization, args.fact))

    fingerprint_file = open("%s" % args.path + "fingerprint_" + args.type_factorization + ".txt", 'w')
    fact_fingerprint_file = None
    if args.fact == 'create':
        # Create file containing factorizations
        fact_fingerprint_file = open("%s" % args.path + "fact_fingerprint_" + args.type_factorization + ".txt", 'w')

    # SPLIT for multiprocessing
    size = int(len(read_lines)/args.n)
    splitted_lines = [read_lines[i:i + size] for i in range(0, len(read_lines), size)]

    with Pool(args.n) as pool:

        type_factorization = args.type_factorization

        # Check type factorization
        factorization = None
        T = None
        #C = None
        if type_factorization == "CFL":
            factorization = CFL
        elif type_factorization == "ICFL":
            factorization = ICFL_recursive
        elif type_factorization == "CFL_ICFL-10":
            factorization = CFL_icfl
            T = 10
        elif type_factorization == "CFL_ICFL-20":
            factorization = CFL_icfl
            T = 20
        elif type_factorization == "CFL_ICFL-30":
            factorization = CFL_icfl
            T = 30
        elif type_factorization == "CFL_COMB":
            factorization = d_cfl
        elif type_factorization == "ICFL_COMB":
            factorization = d_icfl
        elif type_factorization == "CFL_ICFL_COMB-10":
            factorization = d_cfl_icfl
            T = 10
        elif type_factorization == "CFL_ICFL_COMB-20":
            factorization = d_cfl_icfl
            T = 20
        elif type_factorization == "CFL_ICFL_COMB-30":
            factorization = d_cfl_icfl
            T = 30


        func = partial(compute_fingerprint_by_list, args.fact, args.shift, factorization, T)

        fingerprint_lines = []
        fingerprint_fact_lines = []
        for res in pool.map(func, splitted_lines):

            fingerprint_lines = fingerprint_lines + res[0]
            fingerprint_fact_lines = fingerprint_fact_lines + res[1]

        fingerprint_file.writelines(fingerprint_lines)
        if args.fact == 'create':
            fact_fingerprint_file.writelines(fingerprint_fact_lines)

        fingerprint_file.close()

        if args.fact == 'create':
            fact_fingerprint_file.close()

        print('\nCompute fingerprint by list (%s, %s) - stop!' % (args.type_factorization, args.fact))

# Create fingerprint files (args.step = '1f_np') #######################################################################
def experiment_fingerprint_long_reads_step(args):

    # Input FASTA file containing transcripts
    input_fasta = args.fasta

    # Extract of long reads (Format = ID GENE read)
    read_lines = extract_long_reads(name_file=input_fasta)
    print("read_lines SIZE: ", len(read_lines))
    if len(read_lines) == 0:
        print('No reads extracted!')
        exit(-1)

    print('\nCompute fingerprint by list (%s, %s, %s) - start...' % (args.type_factorization, args.fact, args.s_fact))

    fingerprint_file = open("%s" % args.path + "fingerprint_" + args.type_factorization + "_s_" + args.s_fact + ".txt", 'w')
    fact_fingerprint_file = None
    if args.fact == 'create':
        # Create file containing factorizations
        fact_fingerprint_file = open("%s" % args.path + "fact_fingerprint_" + args.type_factorization + "_s_" + args.s_fact + ".txt", 'w')

    # SPLIT for multiprocessing
    size = int(len(read_lines)/args.n)
    splitted_lines = [read_lines[i:i + size] for i in range(0, len(read_lines), size)]

    with Pool(args.n) as pool:

        type_factorization = args.type_factorization

        # Check type factorization
        factorization = None
        T = None
        if type_factorization == "CFL":
            factorization = CFL
        if type_factorization == "ICFL":
            factorization = ICFL_recursive
        if type_factorization == "CFL_ICFL-10":
            factorization = CFL_icfl
            T = 10
        if type_factorization == "CFL_ICFL-20":
            factorization = CFL_icfl
            T = 20
        if type_factorization == "CFL_ICFL-30":
            factorization = CFL_icfl
            T = 30
        if type_factorization == "CFL_COMB":
            factorization = d_cfl
        if type_factorization == "ICFL_COMB":
            factorization = d_icfl
        if type_factorization == "CFL_ICFL_COMB-10":
            factorization = d_cfl_icfl
            T = 10
        if type_factorization == "CFL_ICFL_COMB-20":
            factorization = d_cfl_icfl
            T = 20
        if type_factorization == 'CFL_ICFL_COMB-30':
            factorization = d_cfl_icfl
            T = 30
        elif type_factorization == "CFL_ICFL_CFL":
            factorization = CFL_icfl_cfl
            T = 30

        func = partial(compute_long_fingerprint_by_list, args.fact, factorization, T, s_fact=int(args.s_fact))

        fingerprint_lines = []
        fingerprint_fact_lines = []
        for res in pool.map(func, splitted_lines):

            fingerprint_lines = fingerprint_lines + res[0]
            fingerprint_fact_lines = fingerprint_fact_lines + res[1]

        fingerprint_file.writelines(fingerprint_lines)
        if args.fact == 'create':
            fact_fingerprint_file.writelines(fingerprint_fact_lines)

        fingerprint_file.close()

        if args.fact == 'create':
            fact_fingerprint_file.close()

        print('\nCompute fingerprint by list (%s, %s, %s) - stop!' % (args.type_factorization, args.fact, args.s_fact))

########################################################################################################################
# Create fingerprint files (args.step = 'mapping') #######################################################################
def fingerprint_mapping(args):
    # Input FASTA file containing transcripts
    input_fasta = args.fingerprint

    mapped_file = open("%s" % args.path + "mapped_" + os.path.basename(args.fingerprint).strip('.txt') + ".fa" , 'w')
    mapped_lines = mapping_projection(input_fasta)
    mapped_file.writelines(mapped_lines)
    mapped_file.close()

    # Create a txt copy of the file
    #fing = args.fingerprint.strip('.txt')
    #shutil.copy("%s" % args.path + "mapped_" + fing + '.txt',
    #            "%s" % args.path + "mapped_" + fing + ".fa")


########################################################################################################################

##################################################### MAIN #############################################################
########################################################################################################################
if __name__ == '__main__':

    # Gestione argomenti ###############################################################################################
    parser = argparse.ArgumentParser()

    parser.add_argument('--type', dest='type', action='store', default='1f_np')
    parser.add_argument('--path', dest='path', action='store', default='training/')
    parser.add_argument('--type_factorization', dest='type_factorization', action='store',default='CFL')
    parser.add_argument('--fasta', dest='fasta', action='store', default='transcript_genes.fa')
    parser.add_argument('--fingerprint', dest='fingerprint', action='store', default='prova_fingerprint.txt')
    parser.add_argument('--fact', dest='fact', action='store', default='create')
    parser.add_argument('--shift', dest='shift', action='store', default='shift')
    parser.add_argument('--filter', dest='filter', action='store', default='list')
    parser.add_argument('-n', dest='n', action='store', default=1, type=int)
    parser.add_argument('--s_fact', dest='s_fact', action='store', default='list')

    args = parser.parse_args()

    if args.type == '1f_np':
        print('\nFingerprint Step: 1f_np...\n')
        experiment_fingerprint_1f_np_step(args)

    elif args.type == 'long_reads':
        print('\nFingerprint long reads...\n')
        experiment_fingerprint_long_reads_step(args)

    elif args.type == 'mapping':
        print('\nMapping projecyion of fingerprint files...\n')
        fingerprint_mapping(args)
