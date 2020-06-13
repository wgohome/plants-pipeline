# Module Header ################################################################
import sys
from os.path import realpath, dirname
import re

if __name__ == '__main__':
    abspath = realpath(dirname(__file__))
    parent_module = re.search('^(.*plants-pipeline)', abspath).group()
    sys.path.insert(0, parent_module)
################################################################################

import os
import argparse
import pandas as pd
# Relative imports
from config.constants import DATA_PATH
import helpers
import process_fastq as proc

parser = argparse.ArgumentParser(description = 'This script despatches runids to download.py to run the ascp download and kallisto quantification for each Run ID in parallel.', epilog = 'By Mutwil Lab')
parser.add_argument('-i', '--input', nargs=1, metavar='species_alias',
                    help='Enter the species alias to be downloaded. For instance, Arabidopsis thaliana\'s species alias would be \'Ath\'.',
                    dest='spe', type=str, required=True)
parser.add_argument('-c', '--cds', nargs=1, metavar='cds_filename',
                    help='Enter the file name of the cds .fasta file to use. Do not include the full path, just the filename. It is expected that cds file is places in pipeline-data/download/cds directory.',
                    dest='cds_filename', type=str, required=True)
args = parser.parse_args()
spe = args.spe[0].capitalize()
cds_path = f"{DATA_PATH}/download/cds/{args.cds_filename[0]}"
idx_path = f"{DATA_PATH}/download/idx/{spe}.idx"
runtable_path = helpers.build_runtable_path(spe)
assert os.path.exists(runtable_path), f"The runtable for {spe} is not in pipeline-data/preprocess/out/sra-runtables."

# TODO: Make it robust to SRA inconsistent header names
runs_df = helpers.read_runtable(spe, runtable_path)
runids = runs_df['Run'].iloc[::300][:5]

if __name__ == '__main__':
    if not os.path.exists(idx_path):
        runtime, exit_code = proc.kallisto_index(idx_path=idx_path, cds_path=cds_path)
    proc.process_batch(runids, idx_path, spe, curl_stream=True)
