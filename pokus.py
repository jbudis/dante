import argparse
from arguments.input import nonempty_file
from arguments import yaml_reader
import pandas as pd
import ensembl_rest
import re
import sys
import yaml
import textwrap


def load_arguments():
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

    parser.add_argument('table_file', type=nonempty_file, help='Table of motifs in tsv or csv format')
    parser.add_argument('config_file', type=nonempty_file, help='YAML configuration file')

    parser.add_argument('--flank', type=int, help='Length of flanking sequence')
    parser.add_argument('--seq', type=int, help='The way of typing repetitive sequence with 1 repetition into YAML file.')
    parser.add_argument('-m', type=int, help='Min map quality - I do not think we use this parameter in dante. It is used in BamFilter but the line is in comments')

    parser.add_argument('-u', help='Include unmapped reads', action="store_true")
    parser.add_argument('-r', help='Remove current motifs in YAML configuration file', action="store_true")

    args = parser.parse_args()

    config_dict = yaml_reader.load_arguments(args.config_file)

    if args.r and 'motifs' in config_dict:
        del config_dict['motifs']

    flank = 30
    if args.flank is not None:
        flank = args.flank

    rep_type = 1
    if args.seq is not None:
        rep_type = args.seq

    if args.table_file.endswith('.csv'):
        table_df = pd.read_csv(args.table_file)
    elif args.table_file.endswith('.tsv'):
        table_df = pd.read_csv(args.table_file, sep='\t')
    else:
        print('Unsupported format of table file. Supported formats are csv or tsv.')
        sys.exit()

    return config_dict, table_df, flank, rep_type, args.u, args.config_file, args.m


def parse_nomenclature(nc, flank):
    def get_number(string):
        temp = re.search(r'\d', string)
        if temp is None:
            return None, None

        start_num = temp.start()
        string = string[start_num:]

        temp = re.search(r'[^\d]', string)

        end_num = len(string)
        if temp is not None:
            end_num = temp.start()

        num = int(string[:end_num])
        string = string[end_num:]

        return string, num

    def parse_sequence(string):
        seq_list = []
        digit = False
        seq = ''
        num = ''

        for i in string:
            if i.isalpha():
                if digit:
                    digit = False
                    seq_list.append([seq, int(num)])
                    seq = ''
                    num = ''
                seq += i
            elif i.isdigit():
                digit = True
                num += i

        if seq != '':
            if num == '':
                num = '0'
            seq_list.append([seq, int(num)])

        return seq_list

    nc = re.sub(r'[ ,]', '', nc)
    colon = nc.index(':')
    name = nc[:colon]
    reg_seq = nc[colon:]

    dot = name.index('.')

    chr = re.sub(r'[^\d]', '', name[:dot])
    chr = str(int(chr))

    if chr == '23':
        chr = 'X' # Neviem kedy X a kedy Y

    reg_seq, start_pos = get_number(reg_seq)
    seq_in_table, end_pos = get_number(reg_seq)

    repetitive = None
    if len(seq_in_table) > 0:
        repetitive = parse_sequence(seq_in_table)

    start_pos -= flank
    end_pos += flank

    region = chr + ':' + str(start_pos) + '..' + str(end_pos)

    try:
        temp = ensembl_rest.sequence_region('human', region)
        seq = temp['seq']
    except ensembl_rest.HTTPError as e:
        print('HTTPError: ', e)
        print('Error occurred in region: ', region)
        seq = None

    return seq, repetitive, chr, start_pos, end_pos


def find_seq(string, seq):
    if 'N' not in seq:
        return string.find(seq)
    if seq.count('N') == 1:
        characters = ['A', 'C', 'T', 'G']
        min = len(string)
        for i in characters:
            temp = seq.replace('N', i)
            ind = string.find(temp)
            if ind > -1 and ind < min:
                min = ind

        if min == len(string):
            return -1
        else:
            return min
    return 0


def check_first(seq, left_flank, rep_list):
    if find_seq(seq, rep_list[0][0]) != 0:
        print('Error: Repetitive sequence ', rep_list, 'does not match sequence from referenced genome')
        print('Attempt to repair by skipping this repetition')
        if len(rep_list) > 1 and find_seq(seq, rep_list[1][0]) == 0:
            rep_list[0][1] = 0
        else:
            print('Attempt to repair by ignoring first nucleotides in repetitive sequence and adding them to flank')
            ind = find_seq(seq, rep_list[0][0])
            if ind > -1:
                left_flank = left_flank + seq[:ind]
                seq = seq[ind + len(rep_list[0][0]):]

    return seq, left_flank, rep_list


def check_repetitions(seq, rep_list):
    for i, (rep_seq, rep_num) in enumerate(rep_list):
        if find_seq(seq, rep_seq) != 0 and i > 0:
            print('Error: Repetitive sequence', rep_list, 'does not match sequence from referenced genome')
            print('Attempt to repair by increasing number of previous repetitions')
            prev = rep_list[i - 1]
            while find_seq(seq, prev[0]) == 0:
                seq = seq[len(prev[0]):]
                prev[1] += 1

        for r in range(rep_num):
            if find_seq(seq, rep_seq) == 0:
                seq = seq[len(rep_seq):]
            else:
                print('Error: Repetitive sequence', rep_list, ' does not match sequence from referenced genome')
                print('Attempt to repair by skipping this repetition')
                rep_list[i][1] = r
                break

    if seq != '':
        print('Error: Repetitive sequence', rep_list, ' does not match sequence from referenced genome')
        print('Attempt to repair by increasing number of previous repetitions')
        last = rep_list[-1]
        while find_seq(seq, last[0]) == 0:
            seq = seq[len(last[0]):]
            last[1] += 1

    return seq, rep_list


def check_flanks(left_flank, right_flank, rep_list):
    def check_left_flank(l_flank, rep_list):
        first_seq = rep_list[0][0]
        if l_flank[-len(first_seq):] == first_seq:
            print('Repetitive sequence find in flank region. Adding this sequence into repetitions.')
            while l_flank[-len(first_seq):] == first_seq:
                l_flank = l_flank[:-len(first_seq)]
                rep_list[0][1] += 1
        return l_flank, rep_list

    def check_right_flank(r_flank, rep_list):
        last_seq = rep_list[-1][0]
        if find_seq(r_flank, last_seq) == 0:
            print('Repetitive sequence find in flank region. Adding this sequence into repetitions.')
            while find_seq(r_flank, last_seq) == 0:
                r_flank = r_flank[len(last_seq):]
                rep_list[-1][1] += 1
        return r_flank, rep_list

    left_flank, rep_list = check_left_flank(left_flank, rep_list)
    right_flank, rep_list = check_right_flank(right_flank, rep_list)

    return left_flank, right_flank, rep_list


def make_motifs(desc, full_name, repeat, left_flank, right_flank, rep_type, chr, start, end, unmapped, mapq):
    motif = {'description': desc, 'full_name': full_name}

    still_flank = False
    modules = []

    if rep_type != 2:
        modules.append({'seq': left_flank})
    else:
        still_flank = True
        first_right_flank = 0
        for i, (rep_seq, rep_num) in reversed(list(enumerate(repeat))):
            if rep_num > 1:
                first_right_flank = i + 1
                break

    bases = '5'
    repetitions = '1'

    for i, (rep_seq, rep_num) in enumerate(repeat):
        if rep_num > 1 or (rep_type == 3 and rep_num > 0):
            if still_flank:
                modules.append({'seq': left_flank})
                still_flank = False
            repetitions += ',' + str(int(rep_num / 3))
            bases += ',' + str((int(rep_num / 3)) * len(rep_seq))
            modules.append({'seq': str(rep_num) + '-' + rep_seq})
        elif rep_num == 1:
            if still_flank:
                left_flank = left_flank + rep_seq
            elif seq_type == 2 and i >= first_right_flank:
                right_flank = rep_seq + right_flank
            else:
                modules.append({'seq': rep_seq})
                repetitions += ',0'
                bases += ',0'
    if still_flank:
        modules.append({'seq': left_flank})

    modules.append({'seq': right_flank})

    bases += ',5'
    repetitions += ',1'

    motif['modules'] = modules

    prefilter = {'type': 'BamFilter', 'chromosome': 'chr' + chr, 'ref_start': start, 'ref_end': end, 'include_unmapped': unmapped}

    if mapq is not None:
        prefilter['min_mapq'] = mapq

    motif['prefilter'] = prefilter

    motif['postfilter'] = {'bases': bases, 'repetitions': repetitions}

    return motif


config, table, flank_len, seq_type, unmapped, path, min_mapq = load_arguments()

motifs = []

for index, rows in table.iterrows():
    nomenclature = rows['nomenclature']
    disease = rows['disease']
    description = rows['description']

    ref_sequence, repetitions, chromosome, ref_start, ref_end = parse_nomenclature(nomenclature, flank_len)

    if ref_sequence is None:
        print('Error: Sequence did not find in referenced genome. Ignoring region %s:%d..%d.' % (chromosome, ref_start, ref_end))
        continue
    if repetitions is None:
        print('Error: Repetitive sequences not in table. Ignoring region  %s:%d..%d.' % (chromosome, ref_start, ref_end))
        continue

    sequence = ref_sequence
    s_left = sequence[:flank_len]
    sequence = sequence[flank_len:]
    s_right = sequence[-flank_len:]
    sequence = sequence[:-flank_len]

    sequence, s_left, repetitions = check_first(sequence, s_left, repetitions)

    sequence, repetitions = check_repetitions(sequence, repetitions)

    if sequence != '':
        print('Attempt to repair failed. Adding nucleotides to right flank sequence')
        s_right = sequence + s_right

    s_left, s_right, repetitions = check_flanks(s_left, s_right, repetitions)

    """if len(s_right) < flank_len or len(s_left) < flank_len:
        ref_start = ref_start - flank_len + len(s_left)
        ref_end = ref_end + flank_len - len(s_right)
        region = chromosome + ':' + str(ref_start) + '..' + str(ref_end)
        try:
            temp = ensembl_rest.sequence_region('human', region)
            s_right = temp['seq']
        except ensembl_rest.HTTPError as e:
            print('HTTPError: ', e)
            print('Error occurred in region: ', region) """

    motif = make_motifs(disease, description, repetitions, s_left, s_right, seq_type, chromosome, ref_start, ref_end, unmapped, min_mapq)
    motifs.append(motif)


if unmapped:
    config['general']['include_unmapped'] = True
else:
    config['general']['include_unmapped'] = False

config['motifs'] = motifs

with open(path, 'w') as file:
    yaml.dump(config, file)
