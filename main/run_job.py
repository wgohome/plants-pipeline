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
import numpy as np
import pandas as pd
import os
import wget
import re
from time import time
import pdb
# Relative imports
from config.constants import DATA_PATH
from download import helpers, download_functions, checkfiles
from preprocess.despatch import fetch_runtable

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'This script despatches all species in job-list to run the ascp download and kallisto quantification for each Run ID in parallel.', epilog = 'By Mutwil Lab')
    parser.add_argument('-m', '--method', nargs=1, metavar='download_method',
                        help="This script allows for download methods: 'ascp-bash', 'ascp-python' or 'curl'.", choices=['ascp-bash', 'ascp-python', 'curl'],
                        dest='download_method', type=str, required=True)
    parser.add_argument('-w', '--workers',nargs=1, metavar='num_workers', default=['8'], required=True, dest='workers', help="Specify the number of workers to spawn for multiple process if `-l` is not chosen.")
    parser.add_argument('-t', '--threads',nargs=1, metavar='num_threads', default=['2'], required=True, dest='threads', help="Specify the number of threads to use for each kallisto quantification call.")
    args = parser.parse_args()
    download_method = args.download_method[0].lower()
    workers = int(args.workers[0])
    threads = int(args.threads[0])

def get_valid_jobs():
    # Go through latest job-list
    job_list_dir = f"{DATA_PATH}/preprocess/job-list/"
    latest_job_file = max(os.listdir(job_list_dir))
    raw_job_df = pd.read_csv(f"{job_list_dir}{latest_job_file}", sep='\t')
    job_df = raw_job_df.loc[raw_job_df['cds_link'].notnull(),:]
    taxids = job_df['taxid'].tolist()
    idx_exists = lambda taxid: os.path.exists(f"{DATA_PATH}/download/idx/taxid{taxid}.idx")
    runtable_exists = lambda taxid: len([file for file in os.listdir(f"{DATA_PATH}/preprocess/sra-runtables/") if f"taxid{taxid}" in file])
    ready_taxids = [taxid for taxid in taxids if (idx_exists(taxid) and runtable_exists(taxid))]
    return ready_taxids

def log_status(to_write):
    with open(STATUS_LOG_PATH, 'a') as f:
        f.write(to_write)

if __name__ == '__main__':
    taxids = get_valid_jobs()
    STATUS_LOG_PATH = f"{DATA_PATH}/download/logs/status/{helpers.get_timestamp()}_job.log"
    taxids_str = ', '.join([f"{i}:{taxid}" for i, taxid in enumerate(taxids, 1)])
    log_status(f"# Starting job for {len(taxids)} species.\n# Taxanomic IDs: {taxids_str}\n")
    log_status("index\ttaxid\tattempt\ttotal_downloaded\ttotal_not_downloaded\tfreshly_downloaded\ttime\n")
    for i, taxid in enumerate(taxids, 1):
        # Validate runtable headers first
        idx_path = f"{DATA_PATH}/download/idx/taxid{taxid}.idx"
        for attempt in range(1,4):
            start = time()
            completed_df, incomplete_df = checkfiles.validate_latest_batch(taxid, to_log=False)
            x0 = completed_df.shape[0]
            if incomplete_df.shape[0] == 0:
                continue
            download_functions.process_batch(incomplete_df, idx_path=idx_path, spe_id=f"taxid{taxid}", download_method=download_method, workers=workers, threads=threads)
            completed_df, incomplete_df = checkfiles.validate_latest_batch(taxid, to_log=True)
            checkfiles.update_runinfo_main(taxid)
            log_status(f"{i}\t{taxid}\t{attempt}\t{completed_df.shape[0]}\t{incomplete_df.shape[0]}\t{completed_df.shape[0]-x0}\t{time()-start}\n")
            # log_status(f"Completed Attempt {attempt} for taxid{taxid}, Successful: {completed_df.shape[0]}, Unsuccessful: {incomplete_df.shape[0]}\n")
    log_status("# Completed all valid jobs.\n")
