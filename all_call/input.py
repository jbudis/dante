from __future__ import print_function

from argparse import ArgumentParser, RawDescriptionHelpFormatter
import train
import numpy as np
import json
import sys
import pandas as pd
import re

# default parameters for inference
DEFAULT_MODEL_PARAMS = (-0.0107736, 0.00244419, 0.0, 0.00440608)
DEFAULT_READ_DROP = (479.596, -21.4382)
DEFAULT_READ_DROP_REL = (1.18332, -0.0475454)
DEFAULT_FIT_FUNCTION = "linear"

# functions for training
fit_functions = {"const": train.const_rate, "linear": train.linear_rate, "n2": train.n2_rate, "exp": train.exp_rate}


def load_arguments():
    """
    Loads all arguments and sets default values.
    :return: argparse arguments
    """
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    training = parser.add_argument_group('Training')
    training.add_argument('-t', '--train', action='store_true', help="Train parameters of the model based on the true values. Do not predict values. Default: False")
    training.add_argument('--fit-function', choices=fit_functions.keys(), default="linear", help="Function to approximate deletion rate of STRs. Default: linear")
    training.add_argument('--model-fig', type=str, default=None, help="File to write .png file with comparison of models and train data. Suffix determines the type of image file.")

    inference = parser.add_argument_group('Inference')
    inference.add_argument('--minimal-prob', type=float, default=0.0001, help="Base chance to generate any STR. Counters noise in data at the cost of decreased confidence.")
    inference.add_argument('--minimal-weight', type=float, default=0.15,
                           help="Minimal weight of the model for second allele. Comes into use only when alleles are far from each other, for example (5, 20).")
    inference.add_argument('--pcolor', type=str, default=None,
                           help="File prefix to write image file(s) with pcolor image of confidence of each combination of alleles.")
    inference.add_argument('--relative', action='store_true', help="Use relative read count drop parameters.")
    inference.add_argument('-o', '--output', type=str, default=None, help="File prefix to write the predicted values of alleles. If not provided, output goes to standard output. ")
    inference.add_argument('--index-rep', type=int, default=None, help="Index of the repeating component. If not provided, do all of them.")

    input_args = parser.add_argument_group('Input/Output')
    input_args.add_argument('--profiles', type=str, required=True, help="TSV file or .npy file with one or more profiles. Required.")
    input_args.add_argument('--old-profiles', action='store_true', help="Expect profiles in old style (not Dante's 'repetitions.txt')")
    input_args.add_argument('--true-values', type=str, default=None, help="True values for training of profiles. Required if --train is set.")
    input_args.add_argument('--params-file', type=str, default=None, help="File with parameters of the model (output of training, input of inference).")
    input_args.add_argument('-v', '--verbosity-level', type=int, choices=range(3), default=1, help="Level of verbosity, default 1.")
    input_args.add_argument('--output-profile', type=str, default=None, help="File, where to append the profile, default: None")
    input_args.add_argument('-l', '--len_repeating', type=int, default=3, help="Length of the STR. Used for read drop modelling.")

    args = parser.parse_args()

    # check
    if args.train and args.true_values is None:
        print("ERROR: --train requires --true-values to be set")
        exit(-1)

    if not args.train and args.params_file is None:
        print("WARNING: --params-file is not specified, default parameters are used!!!")

    return args


def read_dante(filename):
    """
    Read profiles from CSV from Dante.
    :param filename: str - filename to read
    :return: Pandas.DataFrame with read profiles or None if no read occurred
    """
    # now try to load tsv file:
    name = filename.split('/')[-2]

    try:
        profiles = pd.read_csv(filename, sep="\t", header=None, index_col=None, parse_dates=True)
    except Exception:
        return None

    new_profiles = pd.DataFrame()

    max_str = max(profiles.max(0)[1:]) + 2

    if profiles is not None:
        for column in profiles.columns[1:]:
            vals = np.zeros(max_str, dtype=int)
            for i, c in enumerate(profiles[column]):
                vals[int(c)] += profiles.iloc[i][0]
            new_profiles['%s_%d' % (name, column - 1)] = vals

    if len(new_profiles.index) > 0:
        profiles = new_profiles.transpose()

    return profiles


def fix_profile_file(filename):
    """
    Fix profile file to be able to read as a tsv.
    :param filename: str - filename to fix
    """
    # read the file
    with open(filename) as f:
        lines = f.readlines()

    # find separator
    sep = '\t' if len(lines[0].split('\t')) >= len(lines[0].split(None)) else None

    # count the number of columns:
    cols = np.zeros_like(lines, dtype=int)
    for i, line in enumerate(lines):
        cols[i] = len(line.split(sep))

    # print with the highest number
    max_cols = max(cols)
    with open(filename, 'w') as f:
        for i, line in enumerate(lines):
            f.write(line.strip())
            # append enough zeros
            for _ in range(max_cols - cols[i]):
                f.write('\t0')
            f.write('\n')


def read_profiles(filename):
    """
    Read profiles from CSV or from .npy file.
    :param filename: str - filename to read
    :return: Pandas.DataFrame with read profiles or None if no read occurred
    """
    # first try to load numpy array
    try:
        profiles = np.load(filename)
    except IOError:
        profiles = None

    if profiles is not None:
        profiles = pd.DataFrame(data=profiles[np.newaxis], index=[int(filename.split('.')[0].split('/')[-1])])

    # now try to load tsv file:
    if profiles is None:
        try:
            fix_profile_file(filename)
            profiles = pd.read_csv(filename, sep='\t', header=None, index_col=0, parse_dates=True)
        except IOError:
            profiles = None

    return profiles


def read_true(filename):
    """
    Read true values from json file or from .true file.
    :param filename: str - json file to read
    :return: dict - values read from the json file or None if no read occurred
    """

    class WrongCountError(Exception):
        pass

    true_values = None
    try:
        with open(filename) as f:
            true_values = json.load(f)
    except Exception:
        pass

    if true_values is None:
        try:
            with open(filename) as f:
                true_values = {}
                for line in f:
                    split = line.split()
                    if len(split) == 3:
                        m = re.search(r'_\d+$', split[0])
                        name = split[0]
                        if m is None:
                            name += '_1'
                        true_values[name] = (int(split[1]), int(split[2]))
                    elif len(split) > 3:
                        raise WrongCountError("Wrong number of parsed elements (expected 3, got %d)" % len(split))

        except Exception as e:
            print('ERROR: ', e)
            return None

    return true_values


def read_params(filename):
    """
    Reads all parameters written with write_params(print_all=True)
    :param filename: str - filename to read parameters from, if None, load default params
    :return: 4-tuple, 2-tuple, function - parameters for model, read count drop, and error function for model distributions
    """
    if filename is None:
        return DEFAULT_MODEL_PARAMS, DEFAULT_READ_DROP, DEFAULT_READ_DROP_REL, DEFAULT_FIT_FUNCTION

    # read 2nd and last line of the file
    with open(filename) as f:
        lines = f.readlines()
        fit_function = lines[1].strip().split()[1]
        split = map(float, lines[-1].strip().split())

    if len(split) < 8:
        print("ERROR: parameters were not read successfully, using defaults!", file=sys.stderr)
        return DEFAULT_MODEL_PARAMS, DEFAULT_READ_DROP, DEFAULT_READ_DROP_REL, DEFAULT_FIT_FUNCTION

    # extract parameters from last line of file
    model_params = tuple(split[0:4])
    read_drop_params = tuple(split[4:6])
    read_drop_params_rel = tuple(split[6:8])

    return model_params, read_drop_params, read_drop_params_rel, fit_function
