import argparse
import re
import sys
import textwrap
import math
from dataclasses import dataclass
from copy import deepcopy

import ensembl_rest
import pandas as pd
import yaml

# imports from dante files
from arguments.input import nonempty_file
from arguments import yaml_reader
from prefiltering.SimpleFilter import base_mapping


@dataclass
class Repetition:
    seq: str
    num: int

    def __repr__(self):
        return self.seq + '-' + str(self.num)


def get_regexp(seq: str):
    tmp = ''.join(base_mapping[ch] for ch in seq)
    return tmp


def load_arguments():
    """
    Loads and parses arguments.
    :return: config_dict - config file in dictionary format.
             table_df - motif table in dataframe format.
             args - parsed arguments.
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''Python program to convert tables and arguments to config.yaml file.
The way of typing repetitive sequence with 1 repetition into YAML file.
    1 (default)
      - seq: left_flank 
      - seq: XXX              (sequences with 1 repetition)                     
      - seq: N-XXX            (sequences with more than 1 repetition)
      - seq: right_flank
    2
      - seq: left_flank + sequences with 1 repetition
      - seq: N-XXX            (sequences with more than 1 repetition)
      - seq: right_flank
    3
      - seq: left_flank
      - seq: 1-XXX            (sequences with 1 repetition)
      - seq: N-XXX            (sequences with more than 1 repetition)
      - seq: right_flank '''))

    # add arguments
    parser.add_argument('table_file', type=nonempty_file, help='Table of motifs in tsv or csv format')
    parser.add_argument('config_file', type=nonempty_file, help='YAML configuration file')

    parser.add_argument('-o', '--output-file', type=str, help='Output file path. Default=config_file', default=None)

    parser.add_argument('--flank', type=int, help='Length of flanking sequence (default=30)', default=30)
    parser.add_argument('--seq', type=int,
                        help='The way of typing repetitive sequence with 1 repetition into YAML file. (default=1)',
                        default=1)

    parser.add_argument('--min-flank-post', '--pmf', type=int, help='Minimum length of flanking sequence to pass postfilter (default=5)',
                        default=5)
    parser.add_argument('--min-flank-post-complex', '--pmfc', type=int,
                        help='Minimum length of flanking sequence to pass postfilter for repetition pairs (default=3)',
                        default=3)
    parser.add_argument('--min-rep-post', '--pmr', type=int, help='Minimum number of repetitions to pass prefilter and postfilter (default=2)',
                        default=2)

    parser.add_argument('-m', '--min-mapq', type=int,
                        help='Minimum mapping quality to pass selection process, default=0', default=0)

    parser.add_argument('-e', '--max-errors', type=int,
                        help='Maximum number of errors in repetitive sequence from table compared to reference sequence (default=1)', default=1)

    parser.add_argument('-u', '--include-unmapped', help='Include unmapped reads', action='store_true')
    parser.add_argument('-r', help='Remove current motifs in YAML configuration file', action='store_true')
    parser.add_argument('--quiet', help='Does not print errors and warning messages', action='store_true')
    parser.add_argument('--postfilter', help='Options to add to postfilter, for example "max_errors: 0.2" for 20%% of errors allowed')
    parser.add_argument('--prefilter', help='Applies prefilter, default is to not use prefilter for BAM inputs', action='store_true')
    parser.add_argument('--nomenclature-column', help='Name of the column for nomenclature (default=nomenclature)', default='nomenclature')
    parser.add_argument('--skip-check', help='Skip checks for the repetition parts. Good if you need to supply non-perfect nomenclatures.',
                        action='store_true')

    args = parser.parse_args()

    # read config file
    config_dict = yaml_reader.load_arguments(args.config_file)

    # remove motifs
    if args.r and 'motifs' in config_dict:
        del config_dict['motifs']

    # read table file
    if args.table_file.endswith('.csv'):
        table_df = pd.read_csv(args.table_file)
    elif args.table_file.endswith('.tsv'):
        table_df = pd.read_csv(args.table_file, sep='\t')
    else:
        print('Unsupported format of table file. Supported formats are csv or tsv.')
        sys.exit()

    return config_dict, table_df, args


def motif_to_string(motif):
    """
    Get string representation of a motif.
    :param motif: dict - motif dictionary
    :return: str - string representation
    """
    ret = motif['full_name'] + ' ('
    ret += ','.join([module['seq'] for module in motif['modules']])
    ret += ')'
    return ret


def get_ref_sequence(chromosome, start_pos, end_pos):
    """
    Get sequence from reference genome.
    :param chromosome: str - chromosome
    :param start_pos: int - reference start position
    :param end_pos: int - reference end position
    :return: seq - str - sequence from reference genome
    """
    region = chromosome + ':' + str(start_pos) + '..' + str(end_pos)

    # use ensembl_rest api to get sequence from reference genome
    try:
        temp = ensembl_rest.sequence_region('human', region)
        seq = temp['seq']
    except ensembl_rest.HTTPError as e:
        print('HTTPError: ', e)
        print('Error occurred in region: ', region)
        return None

    return seq


def parse_nomenclature(nc, flank):
    """
    Parse nomenclature. Get sequence from reference genome, repetitive sequences from table and position in reference genome.
    :param nc: str - nomenclature in format 'NC_chromosome.x:g.start_pos_end_posSEQUENCE[N]'
    :param flank: int - length of flank sequence
    :return: ref_seq - str - sequence from reference genome
             repetitive - list - repetitive sequences from table of motifs
             chromosome - str - chromosome
             start_pos - int - reference start position
             end_pos - int - reference end position
    """

    def get_number(string):
        """
        Get first number from string and remove all characters from string until the end of this number.
        :param string: str - string containing numbers
        :return: string - str - string without characters, that were before the end of number
                 num - int - first number in string
        """
        temp = re.search(r'\d', string)
        # if string does not contain number return None
        if temp is None:
            return None, None

        start_num = temp.start()
        string = string[start_num:]

        temp = re.search(r'[^\d]', string)

        # if there are no other characters than numbers, whole string is number
        end_num = len(string)
        if temp is not None:
            end_num = temp.start()

        num = int(string[:end_num])
        string = string[end_num:]

        return string, num

    def parse_rep_sequence(string):
        """
        Parse repetitive sequence into list.
        :param string: str - string containing repetitive sequence
        :return: seq_list - list(Repetition) - list containing repetitive sequences and their number of occurrences
        """
        seq_list = []
        digit = False
        seq = num = ''

        for i in string:
            if i.isalpha():
                if digit:
                    digit = False
                    seq_list.append(Repetition(seq, int(num)))
                    seq = num = ''
                seq += i
            elif i.isdigit():
                digit = True
                num += i

        if seq != '':
            if num == '':
                num = '1'
            seq_list.append(Repetition(seq, int(num)))

        return seq_list

    # remove all spaces and commas
    nc = re.sub(r'[ ,]', '', nc)

    # split string according first occurrence of colon
    colon = nc.index(':')
    name = nc[:colon]
    reg_seq = nc[colon:]

    # get chromosome - number before dot
    try:
        dot = name.index('.')
    except ValueError:  # dot not found
        dot = None
    chromosome = re.sub(r'[^\d]', '', name[:dot])
    chromosome = str(int(chromosome))

    if chromosome == '23':
        chromosome = 'X'
    if chromosome == '24':
        chromosome = 'Y'

    # get start and end positions
    reg_seq, start_pos = get_number(reg_seq)
    seq_in_table, end_pos = get_number(reg_seq)

    # get repetitive sequence
    repetitive = None
    if len(seq_in_table) > 0:
        repetitive = parse_rep_sequence(seq_in_table)

    # get reference sequence
    start_pos -= flank
    end_pos += flank
    ref_seq = get_ref_sequence(chromosome, start_pos, end_pos)

    return ref_seq, repetitive, chromosome, name, start_pos, end_pos


def check_repetitions(sequence, seq_left, seq_right, list_rep, max_errors=1, verbose=False):
    """
    :param sequence: str - sequence from reference genome
    :param seq_left: str - left flank sequence
    :param seq_right: str - right flank sequence
    :param list_rep: list(Repetition) - list of Repetitions(seq, num)
    :param max_errors: int - maximum number of errors in repetition
    :param verbose: bool - be verbose
    :return: seq_right - str - updated right flank sequence
             rep_list - list(Repetition) - updated list of Repetitions(seq, num)
    """

    rep_list = deepcopy(list_rep)
    error = False
    n_errors = 0

    for i, repetition in enumerate(rep_list):
        rep_seq = repetition.seq

        # check for deletions
        for zac in range(len(rep_seq)):
            if re.match(get_regexp(rep_seq[zac:] + rep_seq[:zac]), sequence):
                rep_seq = rep_seq[zac:] + rep_seq[:zac]
                repetition.seq = rep_seq[zac:] + rep_seq[:zac]
                break

        # check
        if not re.match(get_regexp(rep_seq), sequence):
            if i > 0:
                if verbose:
                    print(f'Error: Repetitive sequence {rep_seq} in {rep_list} does not match sequence from reference genome')
                    print(f'       Attempt to repair by increasing number of previous repetitions')
                # e.g. repetitions: (GCA - 2, ATC - 2), ref_sequence: GCAGCAGCAGCAATCATC,
                # repetitions after repair: (GCA - 4, ATC - 2)
                prev = rep_list[i - 1]
                while re.match(get_regexp(prev.seq), sequence):
                    sequence = sequence[len(prev.seq):]
                    prev.num += 1
                if prev.num > list_rep[i - 1].num:
                    prev.num -= n_errors

            else:
                if verbose:
                    print(f'Error: Repetitive sequence {rep_seq} in {rep_list} does not match sequence from reference genome')

                if len(rep_list) > 1:
                    if verbose:
                        print(f'       Attempt to repair by skipping this repetition')
                    # check if sequence match sequence from second repetition
                    # e.g. repetitions: (GCA - 2, ATC - 2), ref_sequence ATCATC,
                    # repetitions after repair: (GCA - 0, ATC - 2)
                    second_rep = rep_list[1]
                    if re.match(get_regexp(second_rep.seq), sequence):
                        repetition.num = 0
                        continue

                if verbose:
                    print(f'       Attempt to repair by ignoring first nucleotides in repetitive sequence and adding them to flank')
                # check if sequence from first repetitions is in sequence
                # e.g. repetitions: (GCA - 2, ATC - 2), ref_sequence "TCT"GCAGCAATCATC,
                # after_repair: left_flank = left_flank + TCT
                temp = re.search(get_regexp(rep_seq), sequence)
                if temp is not None:
                    ind = temp.start()
                    seq_left = seq_left + sequence[:ind]
                    sequence = sequence[ind + len(rep_seq):]

        n_errors = 0
        for r in range(repetition.num):
            if re.match(get_regexp(rep_seq), sequence):
                sequence = sequence[len(rep_seq):]
            else:
                temp = re.search(get_regexp(rep_seq), sequence)
                if temp is not None and n_errors < max_errors:
                    if verbose:
                        print(f'Error: Repetitive sequence {rep_seq} in {rep_list} does not match sequence from reference genome')
                        print(f'       Attempt to repair by ignoring this error')
                    n_errors += 1
                    error = True
                    if temp.start() <= len(rep_seq):
                        sequence = sequence[temp.start():]
                    else:
                        sequence = sequence[len(rep_seq):]
                else:
                    if verbose:
                        print(f'Error: Repetitive sequence {rep_seq} in {rep_list} does not match sequence from reference genome')
                        print(f'       Attempt to repair by decreasing number of this repetition')
                    # e.g. repetitions: (GCA - 2, ATC - 2), ref_sequence GCAATCATC,
                    # repetitions after repair: (GCA - 1, ATC - 2)
                    repetition.num = r
                    break

    if sequence != '':
        if verbose:
            print(f'Error: Repetitive sequence {rep_list} does not match sequence from reference genome')
            print(f'       Attempt to repair by increasing number of previous repetitions')
        # e.g. repetitions: (GCA - 2, ATC - 2), ref_sequence GCAGCAATCATCATCATC,
        # repetitions after repair: (GCA - 2, ATC - 4)
        last = rep_list[-1]
        while re.match(get_regexp(last.seq), sequence):
            sequence = sequence[len(last.seq):]
            last.num += 1
        if last.num > list_rep[-1].num:
            last.num -= n_errors

    if sequence != '':
        if verbose:
            print(f'       Attempt to repair failed. Adding nucleotides to right flank sequence')
        # e.g. repetitions: (GCA - 2, ATC - 2), ref_sequence GCAGCAATCATC"CTG",
        # after repair: right_flank = CTG + right_flank
        seq_right = sequence + seq_right

    return seq_left, seq_right, rep_list, error


def check_left_flank(seq_left, list_rep, verbose=False):
    """
    Check if end of left flank sequence contains first repetitive sequence.
    :param seq_left: str - left flank sequence
    :param list_rep: list(Repetition) - list of Repetitions(seq, num)
    :param verbose: bool - be verbose
    :return: seq_left - str - updated left flank sequence
             rep_list - list(Repetition) - updated list of Repetitions(seq, num)
    """
    rep_list = deepcopy(list_rep)
    first = rep_list[0]

    if re.match(get_regexp(first.seq), seq_left[-len(first.seq):]):
        if verbose:
            print('Repetitive sequence find in left flank region. Adding this sequence into repetitions.')
        while re.match(get_regexp(first.seq), seq_left[-len(first.seq):]):
            # cut repetitive sequence from flank sequence and add it to the list of repetitive sequences
            seq_left = seq_left[:-len(first.seq)]
            first.num += 1
    return seq_left, rep_list


def check_right_flank(seq_right, list_rep, verbose=False):
    """
    Check if start of right flank sequence contains last repetitive sequence.
    :param seq_right: str - right flank sequence
    :param list_rep: list(Repetition) - list of Repetitions(seq, num)
    :param verbose: bool - be verbose
    :return: seq_right - str - updated right flank sequence
             rep_list - list(Repetition) - updated list of Repetitions(seq, num)
    """
    rep_list = deepcopy(list_rep)
    last = rep_list[-1]

    if re.match(get_regexp(last.seq), seq_right):
        if verbose:
            print('Repetitive sequence find in right flank region. Adding this sequence into repetitions.')
        while re.match(get_regexp(last.seq), seq_right):
            # cut repetitive sequence from flank sequence and add it to the list of repetitive sequences
            seq_right = seq_right[len(last.seq):]
            last.num += 1
    return seq_right, rep_list


def compare_sequences(ref_seq, rep_list, flank, max_errors, skip_check=False, verbose=False):
    """
    Compare (repair) sequence from repetitive area in the table and sequence from reference genome.
    :param ref_seq: str - sequence from reference genome
    :param rep_list: list(Repetition) - list of Repetitions(seq, num)
    :param flank: int - length of flank sequence
    :param max_errors: int - maximum number of errors in repetition
    :param skip_check: bool - skip check of repetitions
    :param verbose: bool - be verbose
    :return: rep_list - list(Repetition) - updated list of Repetitions(seq, num)
             seq_left - str - left flank sequence
             seq_right - str - right flank sequence
    """
    # cut flank regions from reference sequence
    seq_left = ref_seq[:flank]
    current_seq = ref_seq[flank:]
    seq_right = current_seq[-flank:]
    current_seq = current_seq[:-flank]

    error = False
    if not skip_check:
        seq_left, seq_right, rep_list, error = check_repetitions(current_seq, seq_left, seq_right,
                                                                 rep_list, max_errors, verbose)

    # move repetitive sequence from flank sequence to list of repetitive sequences,
    # if repetitive sequence is found in flank region
    seq_left, rep_list = check_left_flank(seq_left, rep_list, verbose)
    seq_right, rep_list = check_right_flank(seq_right, rep_list, verbose)

    return rep_list, seq_left, seq_right, error


def write_bases_repetitions(bases, repetitions, rep_list, min_flank, index, back):
    """
    Write bases and repetitions for postfilter.
    :param bases: list(int) - bases in postfilter - IN-PLACE CHANGE
    :param repetitions: list(int) - repetitions in postfilter - IN-PLACE CHANGE
    :param rep_list: list(Repetition) - sequences for motif
    :param min_flank: int - minimum length of flanking sequence to pass postfilter
    :param index: int - index of current bases/repetition
    :param back: int - move backwards in list?
    """

    j = index

    while min_flank > 0:
        temp = rep_list[j]
        temp_num = max(temp.num, 1)
        bases[j] = temp_num * len(temp.seq)

        if bases[j] > min_flank:
            bases[j] = min_flank

        repetitions[j] = min(math.ceil(min_flank / len(temp.seq)), temp_num)
        min_flank -= bases[j]

        if back:
            j -= 1
        else:
            j += 1


def single_filter(rep_list, min_flank, min_rep):
    """
    Make single postfilters (1 repetitive sequence) and prefilter.
    :param rep_list: list(Repetition) - sequences for motif
    :param min_flank: int - minimum length of flanking sequence to pass postfilter
    :param min_rep: int - minimum number of repetitions to pass postfilter and prefilter
    :return: postfilters - list(dict)
             prefilter - list(str)
    """

    postfilters = []
    prefilters = []

    for i, rep in enumerate(rep_list):
        seq = rep.seq
        num = rep.num
        if num > 0:
            repetitions = [0] * len(rep_list)
            bases = [0] * len(rep_list)

            repetitions[i] = min(min_rep, num)
            bases[i] = repetitions[i] * len(seq)

            # bases before
            write_bases_repetitions(bases, repetitions, rep_list, min_flank, i - 1, True)

            # bases after
            write_bases_repetitions(bases, repetitions, rep_list, min_flank, i + 1, False)

            postfilters.append({'bases': ','.join(map(str, bases)), 'repetitions': ','.join(map(str, repetitions)),
                                'index_rep': i + 1})
            prefilters.append(str(repetitions[i]) + '-' + seq)

    return prefilters, postfilters


def complex_filter(rep_list, min_flank, min_rep, indexes):
    """
    Make complex postfilters (2 repetitive sequences).
    :param rep_list: list(Repetition) - sequences for motif
    :param min_flank: int - minimum length of flanking sequence to pass postfilter
    :param min_rep: int - minimum number of repetitions to pass postfilter
    :param indexes: list(int) - indexes of repetitive sequences
    :return: postfilters - list(dict)
    """
    postfilters = []

    for first, second in zip(indexes[:-1], indexes[1:]):

        repetitions = [0] * len(rep_list)
        bases = [0] * len(rep_list)

        j = first
        # bases between 2 repetitions
        while j <= second:
            temp = rep_list[j]
            temp_num = max(temp.num, 1)
            temp_num = min(temp_num, min_rep)
            bases[j] = temp_num * len(temp.seq)

            repetitions[j] = temp_num
            j += 1

        # bases before
        write_bases_repetitions(bases, repetitions, rep_list, min_flank, first - 1, True)

        # bases after
        write_bases_repetitions(bases, repetitions, rep_list, min_flank, second + 1, False)

        postfilters.append({'bases': ','.join(map(str, bases)), 'repetitions': ','.join(map(str, repetitions)),
                            'index_rep': first + 1, 'index_rep2': second + 1})

    return postfilters


def list_for_filters(modules):
    """
    Make list of modules (repetitive sequences) in better format for making filters.
    :param modules: dict - sequences for motif in dict format
    :return: rep_list - list(Repetition)
    """
    rep_list = []
    for module in modules:
        rep = module['seq']
        if re.match(r'\d', rep) is not None:
            num = int(re.search(r'\d+', rep).group())
        else:
            # 0 means 1 repetition without generating postfilter
            num = 0
        seq = re.search(r'[A-za-z]+', rep).group()
        rep_list.append(Repetition(seq, num))

    return rep_list


def make_filters(modules, min_flank, min_flank_complex, min_rep, no_prefilter):
    """
    Make postfilter and prefilter for config file.
    :param modules: dict - sequences for motif in dict format
    :param min_flank: int - minimum length of flanking sequence to pass postfilter
    :param min_flank_complex: int - minimum length of flanking sequence to pass postfilter for repetition pairs
    :param min_rep: int - minimum number of repetitions to pass postfilter and prefilter
    :param no_prefilter: bool - disable prefilter
    :return: prefilter - dict
             postfilters - list of dict
    """
    rep_list = list_for_filters(modules)

    prefilters, postfilters1 = single_filter(rep_list, min_flank, min_rep)

    prefilter = {'type': 'SimpleFilter', 'seq': ','.join(prefilters)}
    if no_prefilter:
        prefilter = {'type': 'DummyFilter'}

    if len(postfilters1) < 2:
        return prefilter, postfilters1

    indexes = []
    for post in postfilters1:
        indexes.append(int(post['index_rep']) - 1)

    postfilters2 = complex_filter(rep_list, min_flank_complex, min_rep, indexes)

    postfilters = postfilters1 + postfilters2

    return prefilter, postfilters


def make_modules(seq_left, seq_right, rep_list, rep_type):
    """
    Make modules for config file.
    :param seq_left: str - left flank sequence
    :param seq_right: str - right flank sequence
    :param rep_list: list(Repetition) - list of Repetitions(seq, num)
    :param rep_type: int - the way of writing repetitive sequences with 1 repetition
    :return: modules - dictionary
    """
    modules = []
    seq = seq_left

    # write all sequences into modules in user specified format
    for repetition in rep_list:
        if repetition.num > 1 or (rep_type == 3 and repetition.num > 0):
            if seq:
                modules.append({'seq': seq})
                seq = ''
            modules.append({'seq': str(repetition.num) + '-' + repetition.seq})
        elif repetition.num == 1:
            if rep_type == 2:
                seq += repetition.seq
            else:
                if seq:
                    modules.append({'seq': seq})
                    seq = ''
                modules.append({'seq': repetition.seq})

    modules.append({'seq': seq + seq_right})

    return modules


def make_motif(desc, full_name, rep_list, seq_left, seq_right, chromosome, chromosome_version, start, end, args):
    """
    Make motif for config file.
    :param desc: str - motif description
    :param full_name: str - motif full name
    :param rep_list: list(Repetition) - list of Repetitions(seq, num)
    :param seq_left: str - left flank sequence
    :param seq_right: str - right flank sequence
    :param chromosome: str - chromosome
    :param start: int - reference start
    :param end: int - reference end
    :param args: dict/arguments - arguments passed to the program
    :return: motif - dictionary
    """

    motif = {'description': desc, 'full_name': full_name, 'chromosome': 'chr' + chromosome,
             'chromosome_version': chromosome_version, 'ref_start': start, 'ref_end': end,
             'include_unmapped': args.include_unmapped}

    modules = make_modules(seq_left, seq_right, rep_list, args.seq)

    motif['modules'] = modules

    motif['prefilter'], motif['postfilter'] = make_filters(modules, args.min_flank_post, args.min_flank_post_complex,
                                                           args.min_rep_post, not args.prefilter)

    if args.min_mapq is not None:
        motif['min_mapq'] = args.min_mapq

    return motif


def print_report(n_motifs, without_corrections, removed_repetitions, decreased_repetitions, increased_repetitions,
                 errors_in_bases, max_errors, errors_in_table, errors_empty_nom, length_flank, in_config, verbose):
    """
    Print final report.
    :param n_motifs: int - number of all motifs
    :param without_corrections: int - number of motifs without corrections
    :param removed_repetitions: list(Repetition) - list of removed Repetitions(seq, num)
    :param decreased_repetitions: int - number of decreased repetitions
    :param increased_repetitions: int - number of increased repetitions
    :param errors_in_bases: int - number of motifs where error in bases was found
    :param max_errors: int - maximum number of errors in repetitions
    :param errors_in_table: int - number of motifs with invalid parameters in table
    :param errors_empty_nom: int - number of motifs with empty nomenclature field
    :param length_flank: int - length of flank sequence
    :param verbose: bool - be verbose
    :param in_config: int - number of motifs in config file
    """
    print()
    print(f'From available {n_motifs} motifs:')
    print(f'    {n_motifs - len(removed_repetitions) - errors_in_table - errors_empty_nom} converted')
    print(f'        {without_corrections} converted without problems')
    print(f'        {n_motifs - without_corrections - len(removed_repetitions) - errors_in_table - errors_empty_nom} '
          f'converted after error correction')
    print(f'        Corrections:')
    print(f'            {decreased_repetitions} decreased repetitions')
    print(f'            {increased_repetitions} increased repetitions')
    print(f'            {errors_in_bases} matched with maximum {max_errors} errors')
    print(f'        There could have been more than 1 correction in 1 motif')
    print()
    print(f'    {len(removed_repetitions) + errors_empty_nom + errors_in_table} removed due to failed conversion')
    if verbose:
        for ref_seq, reps in removed_repetitions:
            print(f'        Reference sequence: {ref_seq[length_flank:-length_flank]}')
            print(f'        Repetitions in table: {reps}\n')
    print(f'        {errors_in_table} removed due to error in table (invalid reference genome or no repetitions)')
    print(f'        {errors_empty_nom} removed due to empty "nomenclature" column')
    print()
    print(f'{in_config} motifs written into config file -> {n_motifs - len(removed_repetitions) - errors_in_table - in_config - errors_empty_nom} '
          f'converted motifs were removed in deduplication process')


"""
Start of real code
"""
if __name__ == "__main__":
    # load arguments
    config, table, args = load_arguments()
    verbose = not args.quiet
    motifs = []

    without_problem = 0
    rep_dec = 0
    rep_inc = 0
    er_in_bases = 0
    er_in_table = 0
    er_empty_nom = 0
    rep_rem = []

    # iter over dataframe in motif table
    for i, (index, rows) in enumerate(table.iterrows()):
        if i % 100 == 0:
            print(f'Progress: {i:6d}/{len(table):6d}')
        nomenclature = rows[args.nomenclature_column]
        disease = rows['disease']
        gene = rows['gene' if 'gene' in rows else 'Gene']
        description = f'{rows["description"]} (gene: {gene})'

        # skip empty nomenclature
        if str(nomenclature) == 'nan':
            er_empty_nom += 1
            if verbose:
                print(f'Error: Skipping row due to empty "nomenclature" ({nomenclature})')
            continue

        # parse nomenclature (get information from nomenclature)
        ref_sequence, repetitions, chromosome, chromosome_version, ref_start, ref_end = parse_nomenclature(nomenclature, args.flank)

        # print motif name
        if verbose:
            print(f'Motif: {description}:')
            print(f'       ref_sequence={ref_sequence}')
            print(f'       nomenclature={nomenclature}')

        # skip motif if there is no sequence from reference genome or repetitive sequence from table
        if ref_sequence is None:
            if verbose:
                print(f'Error: Sequence did not find in reference genome. Ignoring region {chromosome}:{ref_start}..{ref_end}.')
            er_in_table += 1
            continue
        if repetitions is None:
            er_in_table += 1
            if verbose:
                print(f'Error: Repetitive sequences not in table. Ignoring region {chromosome}:{ref_start}..{ref_end}.')
            continue

        # check if reference sequence and repetitive sequence are same (repair repetitive sequence if they are not same)
        new_repetitions, left_flank, right_flank, seq_error = compare_sequences(ref_sequence, repetitions, args.flank,
                                                                                args.max_errors, args.skip_check,
                                                                                verbose)

        new_errors = False
        under_half = False
        for new_rep, old_rep in zip(new_repetitions, repetitions):
            if new_rep.num < int(old_rep.num / 2):
                rep_rem.append([ref_sequence, repetitions])
                under_half = True
                break
            elif new_rep.num < old_rep.num:
                rep_dec += 1
                new_errors = True
            elif new_rep.num > old_rep.num:
                rep_inc += 1
                new_errors = True

        if under_half:
            continue

        if seq_error:
            er_in_bases += 1
            new_errors = True

        if not new_errors:
            without_problem += 1

        if len(right_flank) < args.flank or len(left_flank) < args.flank:
            # flank sequences are shorter than they are supposed to be, update them
            # it is possible because repetitive sequence can be cut out from flank region
            ref_start = ref_start - args.flank + len(left_flank)
            ref_end = ref_end + args.flank - len(right_flank)

            seq_new = get_ref_sequence(chromosome, ref_start, ref_end)

            if seq_new is not None:
                left_flank = seq_new[:args.flank]
                right_flank = seq_new[-args.flank:]

        # make motif
        motif = make_motif(description, disease, new_repetitions, left_flank, right_flank,
                           chromosome, chromosome_version, ref_start, ref_end, args)

        # append postfilter values
        if args.postfilter is not None and ':' in args.postfilter:
            key, value = args.postfilter.strip().split(':')
            try:
                value = float(value.strip())
            except ValueError:
                pass
            for postfilter in motif['postfilter']:
                postfilter[key.strip()] = value

        # check duplicate modules
        new_motif = True
        for mot in motifs:
            if mot['modules'] == motif['modules']:
                new_motif = False
                if verbose:
                    print(f'Duplicated motif removed: {motif_to_string(motif)}')
                    print(f'          Original motif: {motif_to_string(mot)}')
                break

        if new_motif:
            motifs.append(motif)

    if 'motifs' in config:
        config['motifs'].extend(motifs)
    else:
        config['motifs'] = motifs

    # write config dictionary into .yaml file
    with open(args.output_file if args.output_file is not None else args.config_file, 'w') as file:
        yaml.dump(config, file)

    print_report(table.shape[0], without_problem, rep_rem, rep_dec, rep_inc, er_in_bases,
                 args.max_errors, er_in_table, er_empty_nom, args.flank, len(motifs), verbose)
