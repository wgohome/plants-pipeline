# Module Header ################################################################
import sys
from os.path import realpath, dirname
import re

if __name__ == '__main__':
    abspath = realpath(dirname(__file__))
    parent_module = re.search('^(.*plants-pipeline)', abspath).group()
    sys.path.insert(0, parent_module)
################################################################################

import datetime as dt
import os
import time

def get_timestamp():
    return dt.datetime.now().strftime('%Y%m%d-%H%M%S')

def write_log(to_write, log_path):
    with open(log_path, 'a') as f:
        f.write(to_write)

__all__ = ['get_timestamp', 'write_log']
