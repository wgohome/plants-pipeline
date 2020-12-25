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
from download import download_functions, checkfiles
import helpers

parser = argparse.ArgumentParser(description = 'This script despatches runids to download.py to run the ascp download and kallisto quantification for each Run ID in parallel.', epilog = 'By Mutwil Lab')
parser.add_argument('-s', '--taxid', metavar='taxid',
                        help='Enter the taxanomic id of the species to be downloaded, in quotes. For example, "3702" for Arabidopsis thaliana.',
                        dest='taxid', type=int, required=False)
parser.add_argument('-m', '--method', nargs=1, metavar='download_method',
                    help="This script allows for download methods: 'ascp-bash', 'ascp-python' or 'curl'.", choices=['ascp-bash', 'ascp-python', 'curl'],
                    dest='download_method', type=str, required=True)
parser.add_argument('-l', '--linearmode', action='store_true', default=False, required=False, dest='linearmode', help="Include this optional tag if download is to be in linear mode, else default will be parallel processes. If --method specified is ascp-bash")
parser.add_argument('-w', '--workers',nargs=1, metavar='num_workers', default=['8'], required=False, dest='workers', help="Optional. Specify the number of workers to spawn for multiple process if `-l` is not chosen. Otherwise, this argument will be ignored..")
parser.add_argument('-t', '--threads',nargs=1, metavar='num_threads', default=['2'], required=False, dest='threads', help="Optional. Specify the number of threads to use for each kallisto quantification call.")
args = parser.parse_args()
taxid = args.taxid
download_method = args.download_method[0].lower()
linearmode = args.linearmode
workers = int(args.workers[0])
threads = int(args.threads[0])

idx_path = f"{DATA_PATH}/download/idx/taxid{taxid}.idx"
cds_path = f"{DATA_PATH}/download/cds/taxid{taxid}.cds.fasta"
runtable_exists = lambda taxid: len([file for file in os.listdir(f"{DATA_PATH}/preprocess/sra-runtables/") if f"taxid{taxid}" in file])

# Check that both files exists
if not os.path.exists(idx_path):
    print(f"kallisto index for taxid{taxid} is not present.")
    if os.path.exists(cds_path):
        runtime, exit_code, _ = download_functions.kallisto_index(idx_path=idx_path, cds_path=cds_path)
        print(f"kallisto index for taxid{taxid} has been generated.")
    else:
        print(f"CDS for taxid{taxid} is also not present.")
elif not runtable_exists(taxid):
    print(f"Runtable for taxid{taxid} is not present.")
else:
    completed_df, incomplete_df = checkfiles.validate_latest_batch(taxid, to_log=False)
    print(f"Attempting to download {incomplete_df.shape[0]} undownloaded Run IDs out of a total of {completed_df.shape[0] + incomplete_df.shape[0]} Run IDs for taxid{taxid} ...")
    download_functions.process_batch(incomplete_df, idx_path=idx_path, spe_id=f"taxid{taxid}", download_method=download_method, workers=workers, threads=threads, linear=linearmode)
    completed_df, incomplete_df = checkfiles.validate_latest_batch(taxid, to_log=True)
    checkfiles.update_runinfo_main(taxid)
