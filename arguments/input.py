import os
import argparse
import textwrap
import templates
from parser import ReadFile
import re
import sys
import report
import logging
from annotation import MOTIF_NUCLEOTIDES

import yaml_reader

FILTER_TYPES = ("DummyFilter", "RegexFilter", "LevenshteinFilter", "SimpleFilter")

# TODO test if all validators are correctly implemented

DEFAULT_MISMATCHES = 5
DEFAULT_LEVENSHTEIN = 4


def change_suffix(filename, suffix):
    """
    Change suffix into a new one.
    :param filename: str - filename
    :param suffix: str - string to append
    :return: str - new filename with new suffix
    """
    return filename[:filename.rfind('.')] + suffix


def load_arguments():
    """
    Loads and parses the arguments.
    :return: args - parsed arguments
    """
    description = templates.DANTE_DESCRIPTION
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent(description))

    required = parser.add_argument_group('Required')
    required.add_argument('config_file', type=nonempty_file, help="YAML configuration file.")

    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    script_dir = '/'.join(script_dir.split('/')[:-1])
    config = yaml_reader.load_arguments(args.config_file)
    defaults = yaml_reader.load_arguments('%s/input_data/default.yaml' % script_dir)

    if config is None:
        print('ERROR: Failed to load config file %s' % args.config_file)
        exit(-1)
    if defaults is None:
        print('ERROR: Failed to load default config file "%s/input_data/default.yaml"' % script_dir)
        exit(-1)

    try:
        if not os.path.exists(config['general']['output_dir']):
            os.makedirs(config['general']['output_dir'])
        #basename = args.config_file.split('/')[-1]
        basename = 'config.yaml'
        config_output = '%s/%s' % (config['general']['output_dir'], basename)
        yaml_reader.save_arguments(config, config_output)
        yaml_reader.add_defaults(config, defaults)
        yaml_reader.save_arguments(config, change_suffix(config_output, "_with_defaults.yaml"))
    except yaml_reader.ParameterException as e:
        print("ParameterException: '%s'" % e.message)
        exit(-1)

    return config


def load_arguments_old():
    """
    Loads and parses the arguments.
    :return: args - parsed arguments
    """
    description = templates.DANTE_DESCRIPTION
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent(description))

    required = parser.add_argument_group('Required')
    required.add_argument('output_prefix', type=writable_prefix, help="Prefix of output file names")
    required.add_argument('motif', type=str_motif, help="Motif sequence (e.g., ATCGGA,2-TTG,ATCGAT). THE REPEAT NUMBER MUST BE GREATER THAN 1")
    required.add_argument('read_file', type=nonempty_file, nargs='+', help="Read file(s) to annotate, if no file names provided, input is read from stdin and --file-type is required argument. ")

    prefilter = parser.add_argument_group('Prefilter')
    prefilter.add_argument('--filter-type', choices=FILTER_TYPES, default='SimpleFilter', help="Select the type of prefiltering, defaults to SimpleFilter")
    prefilter.add_argument('-f', '--filter', type=str_motif, default=None,
                           help="Filter sequence. To use multiple filters, add them separated by a comma (e.g., ATTATTATT,5-CTG,ATGCT), defaults to motif. Recommendations:\n"
                                "SimpleFilter -- use only the repeating element with less repetitions than in motif (e.g., 4-CTG).\n"
                                "other filters -- use the whole motif, set additional arguments properly.")
    prefilter.add_argument('--levensh-devs', type=multiple_positive_ints, default=None,
                           help="Allowed number of deviations in LevenshteinFilter for each motif. To set allowed number of deviations for each filter separately, "
                                "add them separated by a comma (e.g., 4,2,4), defaults to %d" % DEFAULT_LEVENSHTEIN)
    prefilter.add_argument('--mismatches', type=positive_int, default=None,
                           help="Number of allowed mismatches in regex filtering, defaults to %d (greatly slows down computation with higher (>7) values)" % DEFAULT_MISMATCHES)

    annotation = parser.add_argument_group('Annotation')
    annotation.add_argument('--delete-prob', type=probability, default=0.01, help="Base delete probability")
    annotation.add_argument('--insert-prob', type=probability, default=0.01, help="Base insert probability")
    annotation.add_argument('--max-delete-skips', type=positive_int, default=2, help="Max number of consequently deleted bases")
    annotation.add_argument('--motif-frequency', type=probability, default=0.01, help="Probability of leaving the start state to first pattern state")
    annotation.add_argument('--snp-chance', type=probability, default=0.02, help="Probability for base being a SNP")
    annotation.add_argument('-p', '--processes', type=positive_nonzero_int, default=1, help="Number of processes for annotation")

    postfilter = parser.add_argument_group('Postfilter')
    postfilter.add_argument('-r', '--required-repetitions', type=multiple_positive_ints, default=None,
                            help="Minimal repetition count for each module,e.g. 0,3,1 means that remains only read with at least 3 repetitions of the motif and one repetition of the last primer")
    postfilter.add_argument('-b', '--required-bases', type=multiple_positive_ints, default=None,
                            help="Minimal number of bases annotated by each module,e.g. 0,10,5 means that remains only read with at least 10 bases annotated as the motif and 5 as the last primer")

    other = parser.add_argument_group('Other')
    other.add_argument('-t', '--file-type', choices=ReadFile.SUPPORTED_FORMATS, help="Type of all of the read files. Required if reads comes directly from stdin.")
    other.add_argument('--stranded', choices=ReadFile.STRANDED_TYPES, default='both',
                       help='Specify, if data are from strand-specific assay. "yes" would examine sequences, '
                            'as defined in read file, "reverse" would examine reverse complemented sequences, "both" would '
                            'examine both of them. Defaults to "both"')
    other.add_argument('--max-reads', type=positive_int, default=None, help="Maximal number of reads to process. Default - process all in the input.")

    args = parser.parse_args()

    return args


def check_arguments(args):
    """
    Check if all arguments are valid and set defaults.
    :param args: argparse arguments
    :return: argparse arguments
    """

    # check unusable options:
    if args.filter_type != "RegexFilter" and args.mismatches is not None:
        warn = "WARNING: --mismatches ignored with %s filter!" % args.filter_type
        report.log_str(warn, priority=logging.WARNING)
    if args.filter_type != "LevenshteinFilter" and args.levensh_devs is not None:
        warn = "WARNING: --levensh-devs ignored with %s filter!" % args.filter_type
        report.log_str(warn, priority=logging.WARNING)
    if len(args.read_file) == 0 and args.file_type is None:
        error = "ERROR: --file-type is required argument when no filename for reads provided (reads are read from stdin)."
        report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)

    # set defaults:
    if args.filter is None:
        args.filter = args.motif
    if args.mismatches is None:
        args.mismatches = DEFAULT_MISMATCHES
    if args.levensh_devs is None:
        args.levensh_devs = [DEFAULT_LEVENSHTEIN]
    if len(args.read_file) == 0:
        args.read_file = [sys.stdin]

    # check if all the: motif, filter, levensh-devs, required-bases, and required-repetitions have the same number:
    length_motif = len(args.motif)
    first_motif = ",".join(map(lambda (s, r): s if r == 1 else "%d-%s" % (r, s), args.motif))
    vals_others = [args.filter, args.levensh_devs, args.required_bases, args.required_repetitions]
    names_others = ["--filter", "--levensh-devs", "--required-bases", "--required_repetitions"]
    for n, v in zip(names_others, vals_others):
        if v is not None and len(v) != 1 and len(v) != length_motif:
            error = "ERROR: %s option (%s) does not have the same number of motifs as motif sequence (%s)." % (n, ",".join(map(str, v)), first_motif)
            report.log_str(error, priority=logging.ERROR)
            raise argparse.ArgumentTypeError(error)

    return args


def multiple_positive_ints(value, max_limit=None):
    """
    Validate if value is composed of multiple positive values, separated by comma
    :param value: code provided by user
    :param max_limit: maximal allowed value of ints, skip validation if None
    :return: list of decoded percents
    """
    items = value.split(',')
    return map(lambda x: positive_int(x, max_limit), items)


def multiple_percents(value):
    """
    Validate if value is composed of multiple positive values from 0 to 100, separated by comma
    :param value: code provided by user
    :return: list of decoded percents
    """
    return multiple_positive_ints(value, 100)


def probability(value):
    """
    Validator for double value in interval <0, 1>
    :param value: String value to validate
    :return: double value or ArgumentTypeError
    """
    try:
        float_value = float(value)
    except ValueError:
        error = "Value %s is not float" % value
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)
    if 0 <= float_value <= 1:
        return float_value
    else:
        error = "Value %s is not in interval <0, 1>" % value
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)


def positive_int(value, max_limit=None):
    """
    Represents positive decimal number, 0 included
    :param value: string value to estimate
    :param max_limit: maximal allowed value, skip validation if None
    :return: integer value, if can be converted into positive int else ArgumentTypeError
    """
    try:
        int_value = int(value)
    except ValueError:
        error = "Value %s is not integer" % value
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)
    if int_value < 0:
        error = "Value %s is not positive integer" % value
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)
    if max_limit and max_limit < int_value:
        error = "Value %s must be lower or equal to %s" % (value, max_limit)
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)
    return int_value


def positive_nonzero_int(value):
    """
    Represents positive decimal number, 0 excluded
    :param value: string value to estimate
    :return: integer value, if can be converted into positive int else ArgumentTypeError
    """
    int_value = positive_int(value)
    if int_value == 0:
        error = "Value %s cannot be 0" % value
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)
    return int_value


def nonempty_file(file_path):
    """
    Checks if the filename in input is non-empty file.
    :param file_path: str - filename to check
    :return: str - filename
    """
    if not os.path.exists(file_path):
        error = "File %s does not exist" % file_path
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)
    if os.path.getsize(file_path) == 0:
        error = "File %s is empty" % file_path
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)
    return file_path


def writable_prefix(prefix):
    """
    Checks if this prefix is writable and exists.
    :param prefix: str - prefix to check
    :return: str - prefix
    """
    directory = os.path.dirname(prefix)
    if not os.path.exists(directory):
        error = "Output directory %s does not exist (%s)" % (directory, prefix)
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)
    if not os.access(directory, os.W_OK):
        error = "Output directory %s is not writable (%s)" % (directory, prefix)
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)
    return prefix


def str_motif(code):
    """
    Split motif code into individual parts, such as primers and STRs
    :param code: motif code, e.g ACACAGT,3-AGT,ACCAC
    :return: sequence and number of repetition for each primer and STR of the motif
    """

    def parse_word(sub_code):
        """
        Split individual parts of motif code into sequence and number of repetitions
        :param sub_code: part of motif code, either primer or STR
        :return: sequence of the part and number of its repetition (1 for primer)
        """
        is_primer = sub_code[0] in MOTIF_NUCLEOTIDES
        if is_primer:
            return sub_code, 1

        return sub_code[sub_code.find('-') + 1:], int(sub_code[:sub_code.find('-')])

    code_upper = code.upper()
    reg = re.compile('^((\d+-)?[{nucls}]+(,(\d+-)?[{nucls}]+)*)$'.format(nucls=''.join(MOTIF_NUCLEOTIDES)))
    if not reg.match(code_upper):
        error = "Illegal motif code: %s" % code
        # report.log_str(error, priority=logging.ERROR)
        raise argparse.ArgumentTypeError(error)

    return [parse_word(word) for word in code_upper.split(',')]


def infer_index_rep(motif, postfilter):
    """
    Infer index of a repeated element from postfilter and motif. Return 2 if not found.
    :param motif: list(tuple(seq, rep)) - motif object
    :param postfilter: list(int)/None - postfilter numbers
    :return: int - 1-based inferred index of a repeated element
    """
    started = False
    if postfilter is None:
        return 2
    for i, ((seq, rep), pf) in enumerate(zip(motif, postfilter)):
        if rep > 1 and pf > 0 and started:
            return i + 1
        if pf > 0:
            started = True
    return 2
