import os
import argparse
import textwrap
import templates
import sys
import report
import logging

import arguments.yaml_reader

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
    try:
        parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                         description=textwrap.dedent(description))

        required = parser.add_argument_group('Config files')
        required.add_argument('config_file', type=nonempty_file, help='YAML configuration file.')
        required.add_argument('--motif-file', type=nonempty_file, help='YAML configuration file with motifs.')

        options = parser.add_argument_group('Options')
        options.add_argument('--start-motif', type=int, help='Starting motif index (from 0). Default=0.', default=0)
        options.add_argument('--max-motifs', type=int, help='Maximal number of motifs to load. Default: All', default=None)
        options.add_argument('--cpu', type=int, help='Overwrite config number of CPUs. Default: as in config file', default=None)
        options.add_argument('--nomenclatures', type=int, help='Number of nomenclature strings to add to reports', default=5)
        # options.add_argument('--skip-annotation', action='store_true', help='Skip annotation and do only inference. Debug purposes only.')

        args = parser.parse_args()
    except argparse.ArgumentTypeError as e:
        print('ERROR: Argument parser error. ' + str(e.message))
        exit(-1)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    script_dir = '/'.join(script_dir.split('/')[:-1])
    config = arguments.yaml_reader.load_arguments(args.config_file)
    defaults = arguments.yaml_reader.load_arguments('%s/input_data/default.yaml' % script_dir)

    if config is None:
        print('ERROR: Failed to load config file %s' % args.config_file)
        exit(-1)
    if defaults is None:
        print('ERROR: Failed to load default config file "%s/input_data/default.yaml"' % script_dir)
        exit(-1)

    # load external motifs
    if args.motif_file is not None:
        motifs = arguments.yaml_reader.load_arguments(args.motif_file)
        if motifs is not None:
            config['motifs'] = motifs
    # change cpu levels
    if args.cpu is not None:
        config['general']['cpu'] = args.cpu
    # change nomenclature counts
    if args.nomenclatures is not None:
        config['general']['nomenclatures'] = args.nomenclatures
    # adjust number of motifs
    config['motifs'] = config['motifs'][args.start_motif:]
    if args.max_motifs is not None:
        num_motifs = len(config['motifs'])
        config['motifs'] = config['motifs'][num_motifs-args.max_motifs:]

    try:
        if not os.path.exists(config['general']['output_dir']):
            os.makedirs(config['general']['output_dir'])
        basename = 'config.yaml'
        config_output = '%s/%s' % (config['general']['output_dir'], basename)
        arguments.yaml_reader.save_arguments(config, config_output)
        arguments.yaml_reader.add_defaults(config, defaults)
        arguments.yaml_reader.save_arguments(config, change_suffix(config_output, "_with_defaults.yaml"))
    except arguments.yaml_reader.ParameterException as e:
        print("ParameterException: '%s'" % e.message)
        exit(-1)

    # skip annotations
    # config['skip_annotation'] = args.skip_annotation

    return config


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
    first_motif = ",".join(map(lambda s, r: s if r == 1 else "%d-%s" % (r, s), args.motif))
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
    return list(map(lambda x: positive_int(x, max_limit), items))


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
