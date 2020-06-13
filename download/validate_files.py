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
import helpers

parser = argparse.ArgumentParser(description = 'This script checks and sweep kallisto outputs and fastq, then creating a timestamped progress log file in pipeline-data/download/logs/progress dir.', epilog = 'By Mutwil Lab')
parser.add_argument('-s', '--spe', nargs=1, metavar='species_alias',
                    help='Enter the species alias to be downloaded. For instance, Arabidopsis thaliana\'s species alias would be \'Ath\'.',
                    dest='spe', type=str, required=True)
args = parser.parse_args()
spe = args.spe[0].capitalize()

def species_progress(spe):
    """Returns 2 lists for completed and incomplete Run IDs respectively, for the species alias, spe, specified.
    Pull out latest log for species and find out completed/incomplete Run IDs.
    If no log present, pull runids from runtable."""
    progress_log_dir = f"{DATA_PATH}/download/logs/progress"
    progress_logs = [file for file in os.listdir(progress_log_dir) if spe in file]
    if progress_logs:
        latest_log = sorted(progress_logs)[-1]
        progress_df = pd.read_csv(f"{progress_log_dir}/{latest_log}", sep='\t')
        completed = progress_df.loc[progress_df['status'] == 'completed']['runid']
        incomplete = progress_df.loc[progress_df['status'] == 'incomplete']['runid']
    else:
        runs_df = helpers.read_runtable(spe)
        incomplete = runs_df['Run']
        completed = pd.Series([], dtype=str)
    return completed.tolist(), incomplete.tolist()

def kallisto_processed(runid):
    """Check if runid has been successfully processed by kallisto, i.e. kallisto output present for runids"""
    run_info_path = f"{DATA_PATH}/download/kallisto-tmp/{runid}/run_info.json"
    if os.path.exists(run_info_path):
        return not '"n_processed": 0' in open(run_info_path,'r').read()
    return False

def build_route(runid):
    dir2 = ""
    if 9 < len(runid) <= 12:
        dir2 = "0" * (12 - len(runid)) + runid[-(len(runid) - 9):] + "/"
    dirs = f"{runid[:6]}/{dir2}/"
    return dirs

def sweep(runid):
    """If kallisto output present,
    move RunID kallisto file to kallisto-tmp, arranged in ENA dir structure;
    delete fastq from fastq-tmp, if applicable."""
    tmp_dir = f"{DATA_PATH}/download/kallisto-tmp/{runid}/"
    target_dir = f"{DATA_PATH}/download/kallisto-out/{build_route(runid)}"
    os.makedirs(target_dir, exist_ok=True)
    _ = shutil.move(tmp_dir, target_dir)
    p_fastq = f"{DATA_PATH}/download/fastq-tmp/{runid}_1.fastq.gz"
    up_fastq = f"{DATA_PATH}/download/fastq-tmp/{runid}.fastq.gz"
    if os.path.exists(p_fastq):
        os.remove(p_fastq)
    elif os.path.exists(up_fastq):
        os.remove(up_fastq)

def validate_latest_batch(spe):
    completed, incomplete = species_progress(spe)
    for runid in incomplete:
        if kallisto_processed(runid):
            sweep(runid)
            completed.append(runid)
            incomplete.remove(runid)
    log_path = helpers.initiate_logfile('progress', ['runid', 'status'], spe=f"{spe}-")
    for runid in completed:
        to_write = f"{runid}\tcompleted\n"
        helpers.write_log(to_write, log_path)
    for runid in incomplete:
        to_write = f"{runid}\tincomplete\n"
        helpers.write_log(to_write, log_path)

if __name__ == '__main__':
    validate_latest_batch(spe)

# Pull out latest log for species and find out completed/incomplete Run IDs
# If no log present, pull runids from runtable
# Check if kallisto output present for runids
# If kallisto output present, move RunID kallisto file to kallisto-tmp, arranged in ENA dir structure
# If kallisto output present, delete fastq from fastq-tmp
# Create log to update completed RunIDs and incomplete RunIDs
