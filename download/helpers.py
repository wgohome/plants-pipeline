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
    return f"{DATA_PATH}/preprocess/sra-runtables/{spe}_sra_runtable.txt"

def read_runtable(spe, runtable_path=None):
    if runtable_path == None:
        runtable_path = build_runtable_path(spe)
    return pd.read_csv(runtable_path, sep=',', header=0, index_col=False,
        dtype='string', usecols=['Run', 'Bytes', 'LibraryLayout'])

def get_fastq_routes(runid):
    """Returns tuple of trailing path of fastq file in vol1/fastq/ server's directory, for paired and unpaired libraries,
    and also file names for paired and unpaired libraries"""
    p_file, up_file = f"{runid}_1.fastq.gz", f"{runid}.fastq.gz"
    dir2 = ""
    if 9 < len(runid) <= 12:
        dir2 = "0" * (12 - len(runid)) + runid[-(len(runid) - 9):] + "/"
    dirs = f"{runid[:6]}/{dir2}{runid}/"
    return f"{dirs}{p_file}", f"{dirs}{up_file}", p_file, up_file

__all__ = ['get_timestamp', 'initiate_logfile', 'write_log', 'build_runtable_path', 'read_runtable', 'get_fastq_routes']
