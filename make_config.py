# imports from dante files
from arguments.input import nonempty_file
from arguments import yaml_reader

# imports from global libraries
import argparse
import pandas as pd
import ensembl_rest
import re
import sys
import yaml
import textwrap
from dataclasses import dataclass


@dataclass
class Repetition:
    seq: str
    num: int

    def __repr__(self):
        return self.seq + '-' + str(self.num)


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

    parser.add_argument('--flank', type=int, help='Length of flanking sequence (default=30)', default=30)
    parser.add_argument('--seq', type=int, help='The way of typing repetitive sequence with 1 repetition into YAML file. (default=1)', default=1)
    parser.add_argument('-m', type=int, help='Min map quality - I do not think we use this parameter in dante. It is used in BamFilter but the line is in comments')

    parser.add_argument('-u', help='Include unmapped reads', action="store_true")
    parser.add_argument('-r', help='Remove current motifs in YAML configuration file', action="store_true")

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

    return config_dict, table_df, args.flank, args.seq, args.u, args.config_file, args.m


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
    dot = name.index('.')
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

    return ref_seq, repetitive, chromosome, start_pos, end_pos


def check_first(sequence, seq_left, rep_list):
    """
    Check if first sequence from Repetition list match sequence from reference genome. If not, try to repair.
    :param sequence: str - sequence from reference genome
    :param seq_left: str - left flank sequence
    :param rep_list: list(Repetition) - list of Repetitions(seq, num)
    :return: sequence - str - sequence from reference genome possibly without first letters
             seq_left - str - (updated) left flank sequence
             rep_list - list(Repetition) - updated list of Repetitions(seq, num)
    """
    first_rep = rep_list[0]
    first_rep_seq = re.sub('N', '.', first_rep.seq)

    # check if sequence match sequence from first repetition
    if not re.match(first_rep_seq, sequence):
        print(f'Error: Repetitive sequence {first_rep_seq} in {rep_list} does not match sequence from referenced genome')

        if len(rep_list) > 1:
            print('Attempt to repair by skipping this repetition')
            # check if sequence match sequence from second repetition
            # e.g. repetitions: (GCA - 2, ATC - 2), ref_sequence ATCATC, repetitions after repair: (GCA - 0, ATC - 2)
            second_rep = rep_list[1]
            second_rep_seq = re.sub('N', '.', second_rep.seq)
            if re.match(second_rep_seq, sequence):
                first_rep.num = 0
                return sequence, seq_left, rep_list

        print('Attempt to repair by ignoring first nucleotides in repetitive sequence and adding them to flank')
        # check if sequence from first repetitions is in sequence
        # e.g. repetitions: (GCA - 2, ATC - 2), ref_sequence "TCT"GCAGCAATCATC, after_repair: left_flank = left_flank + TCT
        temp = re.search(first_rep_seq, sequence)
        if temp is not None:
            ind = temp.start()
            seq_left = seq_left + sequence[:ind]
            sequence = sequence[ind + len(first_rep.seq):]

    return sequence, seq_left, rep_list


def check_repetitions(sequence, seq_right, rep_list):
    """
    Check if sequences from Repetition list match sequence from reference genome. If not, try to repair.
    :param sequence: str - sequence from reference genome
    :param seq_right: str - right flank sequence
    :param rep_list: list(Repetition) - list of Repetitions(seq, num)
    :return: seq_right - str - updated right flank sequence
             rep_list - list(Repetition) - updated list of Repetitions(seq, num)
    """
    for i, repetition in enumerate(rep_list):
        rep_seq = re.sub('N', '.', repetition.seq)
        if not re.match(rep_seq, sequence) and i > 0:
            print(f'Error: Repetitive sequence {rep_list} does not match sequence from reference genome')
            print('Attempt to repair by increasing number of previous repetitions')
            # e.g. repetitions: (GCA - 2, ATC - 2), ref_sequence: GCAGCAGCAGCAATCATC, repetitions after repair: (GCA - 4, ATC - 2)
            prev = rep_list[i - 1]
            prev_seq = re.sub('N', '.', prev.seq)
            while re.match(prev_seq, sequence):
                sequence = sequence[len(prev_seq):]
                prev.num += 1

        for r in range(repetition.num):
            if re.match(rep_seq, sequence):
                sequence = sequence[len(rep_seq):]
            else:
                print(f'Error: Repetitive sequence {rep_seq} in {rep_list} does not match sequence from reference genome')
                print('Attempt to repair by decreasing number of this repetition')
                # e.g. repetitions: (GCA - 2, ATC - 2), ref_sequence GCAATCATC, repetitions after repair: (GCA - 1, ATC - 2)
                repetition.num = r
                break

    if sequence != '':
        print(f'Error: Repetitive sequence {rep_list} does not match sequence from reference genome')
        print('Attempt to repair by increasing number of previous repetitions')
        # e.g. repetitions: (GCA - 2, ATC - 2), ref_sequence GCAGCAATCATCATCATC, repetitions after repair: (GCA - 2, ATC - 4)
        last = rep_list[-1]
        last_seq = re.sub('N', '.', last.seq)
        while re.match(last_seq, sequence):
            sequence = sequence[len(last_seq):]
            last.num += 1

    if sequence != '':
        print('Attempt to repair failed. Adding nucleotides to right flank sequence')
        # e.g. repetitions: (GCA - 2, ATC - 2), ref_sequence GCAGCAATCATC"CTG", after repair: right_flank = CTG + right_flank
        seq_right = sequence + seq_right

    return seq_right, rep_list


def check_left_flank(seq_left, rep_list):
    """
    Check if end of left flank sequence contains first repetitive sequence.
    :param seq_left: str - left flank sequence
    :param rep_list: list(Repetition) - list of Repetitions(seq, num)
    :return: seq_left - str - updated left flank sequence
             rep_list - list(Repetition) - updated list of Repetitions(seq, num)
    """
    first = rep_list[0]
    first_seq = re.sub('N', '.', first.seq)

    if re.match(first_seq, seq_left[-len(first_seq):]):
        print('Repetitive sequence find in flank region. Adding this sequence into repetitions.')
        while re.match(first_seq, seq_left[-len(first_seq):]):
            # cut repetitive sequence from flank sequence and add it to the list of repetitive sequences
            seq_left = seq_left[:-len(first_seq)]
            first.num += 1
    return seq_left, rep_list


def check_right_flank(seq_right, rep_list):
    """
    Check if start of right flank sequence contains last repetitive sequence.
    :param seq_right: str - right flank sequence
    :param rep_list: list(Repetition) - list of Repetitions(seq, num)
    :return: seq_right - str - updated right flank sequence
             rep_list - list(Repetition) - updated list of Repetitions(seq, num)
    """
    last = rep_list[-1]
    last_seq = re.sub('N', '.', last.seq)

    if re.match(last_seq, seq_right):
        print('Repetitive sequence find in flank region. Adding this sequence into repetitions.')
        while re.match(last_seq, seq_right):
            # cut repetitive sequence from flank sequence and add it to the list of repetitive sequences
            seq_right = seq_right[len(last_seq):]
            last.num += 1
    return seq_right, rep_list


def compare_sequences(ref_seq, rep_list, flank):
    """
    Compare (repair) sequence from repetitive area in the table and sequence from reference genome.
    :param ref_seq: str - sequence from reference genome
    :param rep_list: list(Repetition) - list of Repetitions(seq, num)
    :param flank: int - length of flank sequence
    :return: rep_list - list(Repetition) - updated list of Repetitions(seq, num)
             seq_left - str - left flank sequence
             seq_right - str - right flank sequence
    """
    # cut flank regions from reference sequence
    seq_left = ref_seq[:flank]
    current_seq = ref_seq[flank:]
    seq_right = current_seq[-flank:]
    current_seq = current_seq[:-flank]

    # check first repetitive sequence
    current_seq, seq_left, rep_list = check_first(current_seq, seq_left, rep_list)

    # check other repetitive sequences
    seq_right, rep_list = check_repetitions(current_seq, seq_right, rep_list)

    # move repetitive sequence from flank sequence to list of repetitive sequences, if repetitive sequence is found in flank region
    seq_left, rep_list = check_left_flank(seq_left, rep_list)
    seq_right, rep_list = check_right_flank(seq_right, rep_list)

    return rep_list, seq_left, seq_right


def make_motif(desc, full_name, rep_list, seq_left, seq_right, rep_type, chromosome, start, end, unmapped, mapq):
    """
    Make motif to config file.
    :param desc: str - motif description
    :param full_name: str - motif full name
    :param rep_list: list(Repetition) - list of Repetitions(seq, num)
    :param seq_left: str - left flank sequence
    :param seq_right: str - right flank sequence
    :param rep_type: int - the way of writing repetitive sequences with 1 repetition
    :param chromosome: str - chromosome
    :param start: int - reference start
    :param end: int - reference end
    :param unmapped - boolean - include unmapped?
    :param mapq: int - minimum mapping quality
    :return: motif - dictionary
    """

    motif = {'description': desc, 'full_name': full_name}

    # write modules
    modules = []

    still_flank = False
    first_right_flank = 0

    if rep_type != 2:
        modules.append({'seq': seq_left})
    else:
        # find out where is the first sequence which can be written into right flank sequence
        still_flank = True
        for i, repetition in reversed(list(enumerate(rep_list))):
            if repetition.num > 1:
                first_right_flank = i + 1
                break

    bases = '5'
    repetitions = '1'

    # write all sequences into modules in user specified format
    for i, repetition in enumerate(rep_list):
        if repetition.num > 1 or (rep_type == 3 and repetition.num > 0):
            if still_flank:
                modules.append({'seq': seq_left})
                still_flank = False
            repetitions += ',' + str(int(repetition.num / 3))
            bases += ',' + str((int(repetition.num / 3)) * len(repetition.seq))
            modules.append({'seq': str(repetition.num) + '-' + repetition.seq})
        elif repetition.num == 1:
            if still_flank:
                seq_left = seq_left + repetition.seq
            elif seq_type == 2 and i >= first_right_flank:
                seq_right = repetition.seq + seq_right
            else:
                repetitions += ',0'
                bases += ',0'
                modules.append({'seq': repetition.seq})

    if still_flank:
        modules.append({'seq': seq_left})

    modules.append({'seq': seq_right})

    bases += ',5'
    repetitions += ',1'

    motif['modules'] = modules

    # write prefilter
    prefilter = {'type': 'BamFilter', 'chromosome': 'chr' + chromosome, 'ref_start': start, 'ref_end': end, 'include_unmapped': unmapped}
    if mapq is not None:
        prefilter['min_mapq'] = mapq
    motif['prefilter'] = prefilter

    # write postfilter
    # postfilter needs to be edited
    motif['postfilter'] = [{'bases': bases, 'repetitions': repetitions}]

    return motif


if __name__ == "__main__":
    # load arguments
    config, table, flank_len, seq_type, unmapped, path, min_mapq = load_arguments()
    motifs = []

    # iter over dataframe in motif table
    for index, rows in table.iterrows():
        nomenclature = rows['nomenclature']
        disease = rows['disease']
        description = rows['description']

        # parse nomenclature (get information from nomenclature)
        ref_sequence, repetitions, chromosome, ref_start, ref_end = parse_nomenclature(nomenclature, flank_len)

        # skip motif if there is no sequence from reference genome or repetitive sequence from table
        if ref_sequence is None:
            print(f'Error: Sequence did not find in referenced genome. Ignoring region {chromosome}:{ref_start}..{ref_end}.')
            continue
        if repetitions is None:
            print(f'Error: Repetitive sequences not in table. Ignoring region {chromosome}:{ref_start}..{ref_end}.')
            continue

        # check if reference sequence and repetitive sequence are same (repair repetitive sequence if they are not same)
        repetitions, left_flank, right_flank = compare_sequences(ref_sequence, repetitions, flank_len)

        if len(right_flank) < flank_len or len(left_flank) < flank_len:
            # flank sequences are shorter than they are supposed to be, update them
            # it is possible because repetitive sequence can be cut out from flank region
            ref_start = ref_start - flank_len + len(left_flank)
            ref_end = ref_end + flank_len - len(right_flank)

            seq_new = get_ref_sequence(chromosome, ref_start, ref_end)

            if seq_new is not None:
                left_flank = seq_new[:flank_len]
                right_flank = seq_new[-flank_len:]

        # make motif
        motif = make_motif(disease, description, repetitions, left_flank, right_flank, seq_type, chromosome, ref_start, ref_end, unmapped, min_mapq)

        # check duplicate modules
        new_motif = True
        for mot in motifs:
            if mot['modules'] == motif['modules']:
                new_motif = False
                break

        if new_motif:
            motifs.append(motif)

    config['motifs'] = motifs

    # write config dictionary into .yaml file
    with open(path, 'w') as file:
        yaml.dump(config, file)
