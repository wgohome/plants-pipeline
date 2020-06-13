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
import pandas as pd
# relative imports
from config.constants import DATA_PATH

def get_timestamp():
    return dt.datetime.now().strftime('%Y%m%d-%H%M%S')

def initiate_logfile(log_type, headers, spe=''):
    log_path = f"{DATA_PATH}/download/logs/{log_type}/{get_timestamp()}-{spe}{log_type}.log"
    with open(log_path, 'w') as f:
        f.write('\t'.join(headers) + '\n')
    return log_path

def write_log(to_write, log_path):
    with open(log_path, 'a') as f:
        f.write(to_write)

def build_runtable_path(spe):
    return f"{DATA_PATH}/preprocess/out/sra-runtables/{spe}_sra_runtable.txt"

def read_runtable(spe, runtable_path=None):
    if runtable_path == None:
        runtable_path = build_runtable_path(spe)
    return pd.read_csv(runtable_path, sep=',', header=0, index_col=False,
        dtype='string', usecols=['Run', 'Bytes'])

__all__ = ['get_timestamp', 'initiate_logfile', 'write_log', 'build_runtable_path', 'read_runtable']
