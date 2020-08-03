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
import pdb
# Relative imports
from config.constants import DATA_PATH
from download_functions import process_batch, kallisto_index
import helpers

parser = argparse.ArgumentParser(description = 'This script despatches runids to download.py to run the ascp download and kallisto quantification for each Run ID in parallel.', epilog = 'By Mutwil Lab')
parser.add_argument('-s', '--species', nargs=1, metavar='species_id',
                    help='Enter the species id to be downloaded. For instance, Arabidopsis thaliana\'s species id would be \'taxid3702\'.',
                    dest='spe_id', type=str, required=True)
parser.add_argument('-c', '--cds', nargs=1, metavar='cds_filename',
                    help='Enter the file name of the cds .fasta file to use. Do not include the full path, just the filename. It is expected that cds file is places in pipeline-data/download/cds directory.',
                    dest='cds_filename', type=str, required=True)
parser.add_argument('-m', '--method', nargs=1, metavar='download_method',
                    help="This script allows for download methods: 'ascp-bash', 'ascp-python' or 'curl'.",
                    dest='download_method', type=str, required=True)
parser.add_argument('-l', '--linearmode', action='store_true', default=False, required=False, dest='linearmode', help="Include this optional tag if download is to be in linear mode, else default will be parallel processes. If --method specified is ascp-bash")
parser.add_argument('-w', '--workers',nargs=1, metavar='num_workers', default=['8'], required=False, dest='workers', help="Optional. Specify the number of workers to spawn for multiple process if `-l` is not chosen. Otherwise, this argument will be ignored..")
parser.add_argument('-t', '--threads',nargs=1, metavar='num_threads', default=['2'], required=False, dest='threads', help="Optional. Specify the number of threads to use for each kallisto quantification call.")
args = parser.parse_args()
spe_id = args.spe_id[0]
cds_path = f"{DATA_PATH}/download/cds/{args.cds_filename[0]}"
idx_path = f"{DATA_PATH}/download/idx/{spe_id}.idx"
runtable_path = helpers.latest_runtable_path(spe_id)
download_method = args.download_method[0].lower()
linearmode = args.linearmode
workers = int(args.workers[0])
threads = int(args.threads[0])

assert os.path.exists(runtable_path), f"The runtable for {spe_id} is not in pipeline-data/preprocess/sra-runtables."

if download_method not in ['ascp-bash', 'ascp-python', 'curl']:
    raise Exception("--method specified is invalid. Only accepts 'ascp-bash', 'ascp-python' or 'curl'")

# TODO: Make it robust to SRA inconsistent header names

runs_df = helpers.read_runtable(spe_id, runtable_path).iloc[::300][:3] # TODO REMOVE
runs_df['Bytes'] = runs_df.loc[:,'Bytes'].fillna('0')
runs_df['Bytes'] = runs_df.loc[:,'Bytes'].astype(int)

if __name__ == '__main__':
    # Create index file for species if not present
    if not os.path.exists(idx_path):
        assert os.path.exists(cds_path), f"The CDS for {spe_id} is not in pipeline-data/download/cds."
        runtime, exit_code, _ = kallisto_index(idx_path=idx_path, cds_path=cds_path)
    # Run the batch
    process_batch(runs_df, idx_path, spe_id, download_method=download_method, linear=linearmode, workers=workers, threads=threads)
