from __future__ import print_function

from argparse import ArgumentParser, RawDescriptionHelpFormatter
import train
import numpy as np
import json
import sys
import pandas as pd
import re
import os
from glob import glob
from arguments import yaml_reader

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
    parser.add_argument('dir_structure', type=path_exists, help='Directory with multiple Dante results directories. '
                                                                'Each Dante directory has filled "all_profiles.txt" and "all_profiles.true" files.')

    # training.add_argument('--model-fig', type=str, default=None, help="File to write .png file with comparison of models and train data. Suffix determines the type of image file.")

    # parser.add_argument('--profiles', type=str, required=True, help="TSV file or .npy file with one or more profiles. Required.")
    parser.add_argument('--output-params', type=convert_to_absolute, default=None, help="File with parameters of the model to save to. Default: dir_structure/params.txt")
    parser.add_argument('--output-profile', type=convert_to_absolute, default=None, help="File, where to collect all the profiles. Default: dir_structure/all_profiles.txt")
    parser.add_argument('--output-true', type=convert_to_absolute, default=None, help="File, where to collect all the true values. Default: dir_structure/all_profiles.true")

    parser.add_argument('--input-true', type=convert_to_absolute, default=None, help="File, with all the true values. Default: collect from Dante predictions")

    parser.add_argument('--config-dir', type=path_exists, default=None, help="Directory, where to save new config files. Default: without saving")
    parser.add_argument('--fit-function', choices=fit_functions.keys(), default="linear", help="Function to approximate deletion rate of STRs. Default: linear")
    parser.add_argument('-v', '--verbosity-level', type=int, choices=range(3), default=1, help="Level of verbosity, default 1.")
    # input_args.add_argument('-l', '--len_repeating', type=int, default=3, help="Length of the STR. Used for read drop modelling.")

    args = parser.parse_args()

    # check
    if args.output_profile is None:
        args.output_profile = '%s/all_profiles.txt' % args.dir_structure
    if args.output_true is None:
        args.output_true = '%s/all_profiles.true' % args.dir_structure
    if args.output_params is None:
        args.output_params = '%s/params.txt' % args.dir_structure

    return args

def convert_to_absolute(path):
    """
    Converts to absolute path, do not check if exists.
    :param path: str - path
    :return: str - absolute path
    """
    return os.path.abspath(path)

def path_exists(path):
    """
    Checks if the supplied path exists.
    :param path: str - path to a file or dir
    :return: str - absolute path to a file or dir
    """
    try:
        path = convert_to_absolute(path)
    except Exception:
        print('ERROR: %s directory does not exists' % path)
        exit(-1)

    return path


def crawl_dante(dir_structure):
    """
    Crawl Dante dir and collect config, profile, and true_vals files
    :param dir_structure: str - directory above the Dante directory structures, here we start the crawl
    :return: list(str) x3 - list of paths to configs, profiles, and true values
    """

    # read all configs
    configs = glob('%s/*/config.yaml' % dir_structure)

    good_configs = []
    profiles = []
    true_vals = []

    # check if every config has its profiles and true_vals
    for config in configs:
        profile = '%s/all_profiles.txt' % os.path.dirname(config)
        if not os.path.exists(profile):
            print('WARNING: "%s" exists but "%s" does not!!' % (config, profile))
            continue
        true_val = '%s/all_profiles.true' % os.path.dirname(config)
        if not os.path.exists(true_val):
            print('WARNING: "%s" exists but "%s" does not!!' % (config, true_val))
            continue

        # all ok, write them:
        good_configs.append(config)
        profiles.append(profile)
        true_vals.append(true_val)

    return good_configs, profiles, true_vals


def get_name(path):
    """
    Get directory name from path to config/profile/...
    :param path: str - path
    :return: str - directory name without blanks
    """
    directory = path.split('/')[-2]
    directory = directory.replace(' ', '_')
    return directory


def update_config(config_path, save_dir, params_file):
    """
    Create new config file with inputs from the outputs of Dante.
    :param config_path: str - path to the config file
    :param save_dir: str - directory where to save the new config
    :param save_dir: str - directory where to save the new config
    :return: None
    """

    # gather inputs:
    directory = os.path.dirname(config_path)
    inputs = glob('%s/*/annotations*' % directory)
    inputs += glob('%s/*/filtered_primer*' % directory)

    # read the old config:
    config = yaml_reader.load_arguments(config_path)

    # update the config with new inputs
    config['inputs'] = []
    for input in inputs:
        config['inputs'].append({'path': input})

    # update the config with new params
    config['allcall']['param_file'] = params_file

    # add "_retrained" to output dirs
    config['general']['output_dir'] = '%s_retrained' % config['general']['output_dir']

    # write it
    name = get_name(config_path)
    config_name = '%s/%s_config.yaml' % (save_dir, name)
    yaml_reader.save_arguments(config, config_name)


def merge_profiles(profiles, output_file):
    """
    Merge all profiles according to the name of dirs and output them.
    :param profiles: list(str) - list of paths to profiles
    :param output_file: str - output file for merged file
    :return: pd.DataFrame - merged DataFrame with all data
    """
    if len(profiles) == 0:
        return None

    # create empty dataframe
    all_profiles = pd.DataFrame()

    # and fill it
    for profile in profiles:
        name = get_name(profile)

        # get the maximal number of columns:
        max_cols = 0
        with open(profile) as f:
            for line in f:
                max_cols = max(max_cols, line.count('\t'))

        # write to aggregated file:
        current = pd.read_csv(profile, sep='\t', header=None, names=['index'] + range(max_cols), index_col=0, parse_dates=True, engine='python')
        current.index = map(lambda x: '%s_%s' % (name, x), current.index)
        all_profiles = pd.concat([all_profiles, current])

    # fill not available data:
    all_profiles = all_profiles.fillna(0).astype(int)
    all_profiles.sort_index(inplace=True)

    # save it:
    all_profiles.to_csv(output_file, sep='\t')

    # return it
    return all_profiles


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
