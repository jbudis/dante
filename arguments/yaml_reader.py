from __future__ import division, print_function

import yaml


class ParameterException(Exception):
    pass


def load_arguments(yaml_file):
    """
    Reads the yaml file and creates arguments for the rest of the program.
    :param yaml_file: str - path to yaml parameters file
    :return: dict - dictionary with parameters
    """
    read_yaml = None

    try:
        with open(yaml_file, 'r') as stream:
            try:
                read_yaml = yaml.load(stream)
            except yaml.YAMLError as exc:
                print(exc)
    except IOError:
        return None

    return read_yaml


def save_arguments(config, yaml_file=None):
    """
    Saves the current arguments into a yaml file.
    :param config: dict - dictionary with parameters
    :param yaml_file: str/stream/None - path/stream to yaml parameters file
    :return: None
    """
    if isinstance(yaml_file, str):
        with open(yaml_file, 'w') as stream:
            return yaml.dump(config, stream, default_flow_style=False)
    else:
        return yaml.dump(config, yaml_file, default_flow_style=False)


def add_defaults(config, default):
    """
    Look through the parameters and adds the default ones recursively.
    :param config: dict - dictionary with parameters
    :param default: dict - dictionary with default parameters
    :return: None
    """

    for key in default:
        if isinstance(default[key], dict):
            # go deeper
            if key not in config:
                config[key] = {}
            add_defaults(config[key], default[key])
        elif isinstance(default[key], list):
            # add it to all the configs
            if key in config:
                for item in config[key]:
                    add_defaults(item, default[key][0])
        elif isinstance(default[key], str) and default[key] == 'required':
            # required - check if it is available
            if key not in config:
                raise ParameterException('Required argument not available - %s' % key)
            if config[key] is None or config[key].strip() == "":
                raise ParameterException('Required argument is empty - %s' % key)
        else:
            # add it
            if key not in config:
                config[key] = default[key]
