# Module Header ################################################################
import sys
from os.path import realpath, dirname
import re

if __name__ == '__main__':
    abspath = realpath(dirname(__file__))
    parent_module = re.search('^(.*plants-pipeline)', abspath).group()
    sys.path.insert(0, parent_module)
################################################################################

import argparse
import os
import shutil
import pandas as pd
import pdb
# Relative imports of CONSTANTS in config/constants.py
from config.constants import DATA_PATH
from download import helpers

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'This script checks and sweep kallisto outputs and fastq, then creating a timestamped progress log file in pipeline-data/download/logs/progress dir.', epilog = 'By Mutwil Lab')
    parser.add_argument('-t', '--taxid', metavar='taxid',
                        help='Enter the taxid of the species to be downloaded. For instance, Arabidopsis thaliana\'s taxid would be \'3702\'.',
                        dest='taxid', type=int, required=True)
    args = parser.parse_args()
    taxid = args.taxid

def kallisto_present(runid):
    return os.path.exists(f"{DATA_PATH}/download/kallisto-out/{helpers.get_dirs(runid)}{runid}.zip")

def log_progress(completed, incomplete):
    log_path = helpers.initiate_logfile('progress', ['runid', 'status'], spe=f"taxid{taxid}-")
    for runid in completed:
        to_write = f"{runid}\tcompleted\n"
        helpers.write_log(to_write, log_path)
    for runid in incomplete:
        to_write = f"{runid}\tincomplete\n"
        helpers.write_log(to_write, log_path)

def validate_latest_batch(taxid, to_log=False):
    runs_df = helpers.read_runtable(f"taxid{taxid}")
    completed = [runid for _, runid in runs_df['Run'].iteritems() if kallisto_present(runid)]
    incomplete = list(set(runs_df['Run']) - set(completed))
    if to_log:
        log_progress(completed, incomplete)
    completed_df = runs_df[runs_df['Run'].isin(completed)]
    incomplete_df = runs_df[runs_df['Run'].isin(incomplete)]
    return completed_df, incomplete_df

def update_runinfo_main(taxid):
    runinfo_main = f"{DATA_PATH}/download/runinfo-main/taxid{taxid}_runinfo_main.txt"
    runinfo_logs = sorted([file for file in os.listdir(f"{DATA_PATH}/download/logs/runinfo") if f"taxid{taxid}" in file])
    assert runinfo_logs, "No runinfo logs, you have not attempted any download for this species?"
    log_dfs = [pd.read_csv(f"{DATA_PATH}/download/logs/runinfo/{log}", sep='\t') for log in runinfo_logs]
    main_df = pd.concat(log_dfs, ignore_index=True)
    main_df.drop_duplicates('runid',keep='last',inplace=True,ignore_index=True)
    main_df.to_csv(runinfo_main, sep='\t', index=None, mode='w')
    return main_df

__all__ = ['update_runinfo_main', 'validate_latest_batch']
