import configparser
import os


def _build_dict():
    """
    Turn config object into flat dictionary (doesn't keep sections).

    :return: Parsed flat config file as dictionary
    """
    _obj = {}
    for section in conf.sections():
        options = conf.options(section)
        for option in options:
            _obj[option] = conf.get(section, option)
    return _obj


# Read config file
conf = configparser.ConfigParser()
# Get parent dir, then append config.ini
conf.read(os.path.abspath(os.sep.join(__file__.split(os.sep)[:-2]) + os.sep + 'config.ini'))
_dict = _build_dict()


def get(key, fallback=None):
    """
    Retrieve value from config file, wrapper for config dictionary get()

    :param key: Key to get value for, will be converted to lowercase
    :param fallback: Optional fallback value if key doesn't exist
    :return: Value from config dict or fallback
    """
    return _dict.get(key.lower(), fallback)
