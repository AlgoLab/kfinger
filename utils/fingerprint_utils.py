import gzip

from .factorizations import CFL
from .factorizations import ICFL_recursive
from .factorizations import CFL_icfl
from .factorizations_comb import d_cfl
from .factorizations_comb import d_icfl
from .factorizations_comb import d_cfl_icfl
from .factorizations_comb import reverse_complement

# Given a list of lengths computes the list of k-fingers
# Given a k_finger, return the enriched string
def get_enrich_str(facts_list):
    enrich_str = None
    if len(facts_list) > 2:
        facts_list.pop(0)
        facts_list.pop()

        if len(facts_list) == 1:
            # Only 1 element
            if len(facts_list[0]) <= 20:
                enrich_str = reverse_complement(facts_list[0])
            else:
                enrich_str = (facts_list[0])[:10] + (facts_list[0])[(len(facts_list[0]) - 10):]
                enrich_str = reverse_complement(enrich_str)
        else:
            # More elements
            max = 0
            max_fact = ''
            for fact in facts_list[::-1]:
                if len(fact) > max:
                    max = len(fact)
                    max_fact = fact
            enrich_str = reverse_complement(max_fact)
            if len(enrich_str) <= 20:
                enrich_str = reverse_complement(enrich_str)
            else:
                enrich_str = (enrich_str)[:10] + (enrich_str)[(len(enrich_str) - 10):]
                enrich_str = reverse_complement(enrich_str)

    # Padding of length 20 for the enriched string
    if len(enrich_str) <= 20:
        for i in range(len(enrich_str), 20):
            enrich_str += 'N'
    return enrich_str


# Given a k_finger, returns the normalized version:
def normalize(k_finger):
    k_finger_inv = k_finger[::-1]

    for a, b in zip(k_finger, k_finger_inv):
        a = int(a)
        b = int(b)

        if a < b:
            return k_finger
        elif b < a:
            return k_finger_inv
        else:
            continue

    return k_finger


# Shift of size d of a string s
def shift_string(string='', size=100, shift='no_shift'):
    list_of_shift = []

    if shift == 'no_shift':
        list_of_shift.append(string)
    else:
        if len(string) < size:
            list_of_shift.append(string)
        else:
            for i in range(len(string)):
                fact = string[i:i + size]
                if i + size > len(string):
                    pref = string[:(i + size) - len(string)]
                    fact = fact + pref
                list_of_shift.append(fact)
    return list_of_shift


# Split long reads in subreads
def factors_string(string='', size=300):
    list_of_factors = []

    if len(string) < size:
        list_of_factors.append(string)
    else:
        # print('len(string): ', len(string))
        for i in range(0, len(string), size):
            # print("i: ", i)
            if i + size > len(string):
                fact = string[i:len(string)]
            else:
                fact = string[i:i + size]

            list_of_factors.append(fact)

    return list_of_factors


# Read file FQ containing reads
def read_fq(fasta_lines):
    lines = []
    read = ''
    step = ''

    i = 0
    while True:
        # ID_GENE
        l_1 = str(fasta_lines[i])
        l_1 = l_1.replace('@m54329U_', '')
        l_1 = l_1.replace('/ccs', '')
        id_gene = l_1
        read = read + id_gene + ' '

        # Read
        l_2 = str(fasta_lines[i + 1])
        read = read + l_2

        lines.append(read)
        read = ''

        i += 4
        if i == len(fasta_lines):
            break

    return lines


# Read file FQ containing reads
# For each read consider two sequences:
#   1) ID 1 original_sequence
#   2) ID 0 R&C_sequence
def read_fq_2_steps(fasta_lines):
    lines = []
    read_original = ''
    read_rc = ''
    step = ''

    i = 0
    while True:
        # ID_GENE
        l_1 = str(fasta_lines[i])
        l_1 = l_1.replace('@m54329U_', '')
        l_1 = l_1.replace('/ccs', '')
        l_1 = l_1.replace('\n', '')
        id_gene = l_1
        read_original = read_original + id_gene + '_1 '
        read_rc = read_rc + id_gene + '_0 '

        # Read
        l_2 = str(fasta_lines[i + 1])
        read_original = read_original + l_2
        read_rc = read_rc + reverse_complement(l_2.replace('\n', ''))

        lines.append(read_original)
        lines.append(read_rc)
        read_original = ''
        read_rc = ''

        i += 4
        if i == len(fasta_lines):
            break

    return lines


# Read file FQ containing reads
def read_long_fasta(fasta_lines):
    lines = []
    read = ''
    step = ''

    i = 0
    while True:
        # ID_GENE
        l_1 = str(fasta_lines[i])
        l_1 = l_1.replace('>', '')
        id_gene = l_1
        read = read + id_gene + ' '

        # Read
        l_2 = str(fasta_lines[i + 1])
        read = read + l_2

        lines.append(read)
        read = ''

        i += 2
        if i == len(fasta_lines):
            break

    return lines


# Read file FQ containing reads
# For each read consider two sequences:
#   1) ID 1 original_sequence
#   2) ID 0 R&C_sequence
def read_long_fasta_2_steps(fasta_lines):
    lines = []
    read_original = ''
    read_rc = ''
    step = ''

    i = 0
    while True:
        # ID_GENE
        l_1 = str(fasta_lines[i])
        l_1 = l_1.replace('>', '')
        l_1 = l_1.replace('\n', '')
        id_gene = l_1
        read_original = read_original + id_gene + '_1 '
        read_rc = read_rc + id_gene + '_0 '

        # Read
        l_2 = str(fasta_lines[i + 1])
        read_original = read_original + l_2
        read_rc = read_rc + reverse_complement(l_2.replace('\n', ''))

        lines.append(read_original)
        lines.append(read_rc)
        read_original = ''
        read_rc = ''

        i += 2
        if i == len(fasta_lines):
            break

    return lines


# Read file GZ containing reads
def read_gz(fasta_lines):
    lines = []
    read = ''
    step = ''

    i = 0
    while True:
        print(i)
        print(fasta_lines)
        # ID_GENE
        l_1 = str(fasta_lines[i])
        l_1 = l_1.replace('b\'', '')
        s_l1 = l_1.split()
        if len(s_l1) == 2:
            id_gene = s_l1[1]
            id_gene = id_gene.replace('\\n', '')
            id_gene = id_gene.replace('\'', '')
            read = read + id_gene + ' '
        else:
            s_l1 = l_1.split(',')
            read = read + s_l1[1] + ' '

        # Read
        l_2 = str(fasta_lines[i + 1])
        l_2 = l_2.replace('b\'', '')
        l_2 = l_2.replace('\\n', '')
        l_2 = l_2.replace('\'', '')
        read = read + l_2

        lines.append(read)
        read = ''

        i += 4
        if i == len(fasta_lines):
            break

    return lines


# Read file FASTA
def read_fasta(fasta_lines):
    lines = []
    read = ''
    step = ''
    for s in fasta_lines:
        if s[0] == '>':
            if read != '':
                lines.append(read)
                read = ''

            s = s.replace('\n', '')
            s_list = s.split()
            read = s_list[1] + ' '
        else:
            s = s.replace('\n', '')
            read += s

    return lines


# Given a FASTA file returns the list containing ONLY the reads, for each read, the corresponding fingerprint
def extract_reads(name_file, filter='list', n_for_genes=None):
    print('\nExtract reads - start...')

    # Creazione  dizionario per cpntare i geni trovati
    list_id_genes = None
    dict_n_for_genes = None
    if n_for_genes != None:
        file_experiment = open('list_experiment.txt')
        list_id_genes = file_experiment.readlines()
        list_id_genes = [s.replace('\n', '') for s in list_id_genes]
        dict_n_for_genes = {i: 0 for i in list_id_genes}
        file_experiment.close()

    file = None
    lines = []

    # Scrittura su file
    if name_file.endswith('.gz'):
        # GZ FILE
        file = gzip.open(name_file, 'rb')
        lines = read_gz(file.readlines())
    elif name_file.endswith('.fa') or name_file.endswith('.fasta') or name_file.endswith('.fastq'):
        # FASTA FILE
        file = open(name_file)
        lines = read_fasta(file.readlines())

    read_lines = []

    cont_genes = 0
    cont_ripper = 0

    for s in lines:

        new_line = ''
        new_fact_line = ''

        str_line = s.split()
        id_gene = str_line[0]
        id_gene = id_gene.replace('\n', '')

        # Create lines
        file_experiment = open('list_experiment.txt')
        list_id_genes = file_experiment.readlines()
        list_id_genes = [s.replace('\n', '') for s in list_id_genes]
        if filter == 'list':
            if id_gene in list_id_genes:

                ########################################################################################################
                if n_for_genes != None:
                    n_for_id_gene = dict_n_for_genes[id_gene]
                    if n_for_id_gene < n_for_genes:
                        lbl_id_gene = id_gene + ' '

                        # UPPER fingerprint
                        sequence = str_line[1]

                        sequence = sequence.upper()
                        new_line = lbl_id_gene + ' ' + sequence
                        new_line += '\n'
                        read_lines.append(new_line)

                        n_for_id_gene += 1
                        dict_n_for_genes[id_gene] = n_for_id_gene
                else:
                    ########################################################################################################
                    lbl_id_gene = id_gene + ' '

                    sequence = str_line[1]
                    new_line = lbl_id_gene + ' ' + sequence
                    new_line += '\n'
                    read_lines.append(new_line)

        elif filter == 'list_ripper':
            ID_GENE_RIPPER = 'RIPPER'

            if id_gene in list_id_genes:

                lbl_id_gene = id_gene + ' '

                sequence = str_line[1]
                new_line = lbl_id_gene + ' ' + sequence
                new_line += '\n'
                read_lines.append(new_line)
                cont_genes += 1
            else:
                if cont_ripper < cont_genes:
                    lbl_id_gene = ID_GENE_RIPPER + ' '

                    sequence = str_line[1]
                    new_line = lbl_id_gene + ' ' + sequence
                    new_line += '\n'
                    read_lines.append(new_line)
                    cont_ripper += 1

        else:
            lbl_id_gene = id_gene + ' '
            new_line = lbl_id_gene + str_line[1]
            read_lines.append(new_line)

        file_experiment.close()

    file.close()

    print('\nExtract reads - stop!')

    return read_lines

# Mapping projection of fingerprints in a file
def mapping_projection(fingerprint_file_path):
    file = open(fingerprint_file_path)
    lines = file.readlines()
    mapped_lines = []

    for line in lines:
        line = line.replace('|','')
        str_line = line.split()
        id_gene = str_line[0]
        fingerprint = [int(i) for i in str_line[1:]]
        mapped_fingerprint = fingerprint_projection(fingerprint)
        str_mapped_fingerprint = '>' + id_gene + '\n' + mapped_fingerprint + "\n"
        mapped_lines.append(str_mapped_fingerprint)
    return mapped_lines

# Mapping fingerprint of a string
def fingerprint_projection(fingerprint: str):
    #alf = '#0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZàèìòù%^-+:_çé!?£$§°*'
    alf = 'ABCDEEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz!"#$%&()*+,-./0123456789:;<=?@[\]^_`{|}~¡¢£¤¥¦§¨©ª«¬®¯°±²³´µ¶·¸¹º»¼½¾¿ÀÁÂÃÄÆÇÈÉËÌÍÎÏÐÑÒÓÔÕÖ×ØÙÚÛÜÝÞßàáâãäæçèéëìíîïðñòóôõö÷øùúûüýþÿ€ƒ„…†‡ˆ‰Š‹ŒŽ•–—˜™š›ƀƁƂƃƆƉƊƋƌƍƎƏƐƑƒƓƔƕƗƘƙƚƛƜƝƞƟƠơƢƣƤƥƦƧƨƩƪƫƬƭƮƯưƱƲƳƴƵƶ•‣‰※⁅⁆⁇⁈⁉⁊⁋⁌⁍⁎⁏⁐⁑⁒⁔⁕⁖⁗⁘⁙⁚⁛⁜⁞⁝₨₯₵₪₺₼₽₾₿'
    mapped_fingerprint: str = ''.join([alf[f] for f in fingerprint])
    return mapped_fingerprint


# Given a FASTA file returns the list containing ONLY the reads, for each read, the corresponding fingerprint
def extract_long_reads(name_file):
    print('\nExtract long reads - start...')

    lines = []

    file = open(name_file)
    # lines = read_fq_2_steps(file.readlines())
    # lines = read_fq(file.readlines())
    # lines = read_long_fasta(file.readlines())
    lines = read_long_fasta_2_steps(file.readlines())

    read_lines = []

    i = 0
    for s in lines:
        print(i)

        i += 1

        str_line = s.split()

        if len(str_line[1]) >= 0:
            id_gene = str_line[0]
            id_gene = id_gene.replace('\n', '')

            lbl_id_gene = id_gene + ' '

            # UPPER fingerprint
            sequence = str_line[1]

            sequence = sequence.upper()
            new_line = lbl_id_gene + ' ' + sequence
            new_line += '\n'
            read_lines.append(new_line)

    print('\nExtract long reads - stop!')

    return read_lines


# Given a list of reads and a factorization technique, compute the list containing, for each read, the corresponding fingerprint
def compute_fingerprint_by_list(fact_file='no_create', shift='no_shift', factorization=CFL, T=None, list_reads=[]):
    fingerprint_lines = []
    fingerprint_fact_lines = []

    id_gene = ''
    step = ''
    for s in list_reads:

        str_line = s.split()

        id_gene = str_line[0]
        read = str_line[1]

        list_of_shifts = shift_string(read, 100, shift)
        for sft in list_of_shifts:
            list_fact = factorization(sft, T)

            # Remove special characters
            if '>>' in list_fact:
                list_fact[:] = (value for value in list_fact if value != '>>')
            if '<<' in list_fact:
                list_fact[:] = (value for value in list_fact if value != '<<')

            # Create lines
            lbl_id_gene = id_gene + ' '
            new_line = lbl_id_gene + ' '.join(str(len(fact)) for fact in list_fact)
            new_line += '\n'
            fingerprint_lines.append(new_line)
            if fact_file == 'create':
                new_fact_line = lbl_id_gene + ' '.join(fact for fact in list_fact)
                new_fact_line += '\n'
                fingerprint_fact_lines.append(new_fact_line)

    return fingerprint_lines, fingerprint_fact_lines


# Given a list of reads and a factorization technique, compute the list containing, for each read, the corresponding fingerprint
def compute_long_fingerprint_by_list(fact_file='no_create', factorization=CFL, T=None, list_reads=[],s_fact=300, C=None):
    fingerprint_lines = []
    fingerprint_fact_lines = []

    id_gene = ''
    step = ''
    for s in list_reads:

        str_line = s.split()

        id_gene = str_line[0]
        read = str_line[1]

        # Create lines
        lbl_id_gene = id_gene + ' '
        new_line = lbl_id_gene + ' '
        new_fact_line = lbl_id_gene + ' '
        list_of_factors = factors_string(read,s_fact)
        for sft in list_of_factors:
            list_fact = factorization(sft,T)

            # Remove special characters
            if '>>' in list_fact:
                list_fact[:] = (value for value in list_fact if value != '>>')
            if '<<' in list_fact:
                list_fact[:] = (value for value in list_fact if value != '<<')

            new_line = new_line + ' '.join(str(len(fact)) for fact in list_fact)
            new_line += ' | '
            if fact_file == 'create':
                new_fact_line = new_fact_line + ' '.join(fact for fact in list_fact)
                new_fact_line += ' | '

        new_line += '\n'
        new_fact_line += '\n'
        fingerprint_lines.append(new_line)
        fingerprint_fact_lines.append(new_fact_line)

    return fingerprint_lines, fingerprint_fact_lines

