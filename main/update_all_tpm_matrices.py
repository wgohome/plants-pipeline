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
import pandas as pd
import concurrent.futures
import re
import pdb
# Relative imports
from config.constants import DATA_PATH
from download import helpers
from postprocess import pull_tpm_matrix
from main import run_job

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'This script gets all downloaded experiments for all species, and creates a tpm matrices for it in data-pipeline/post-process/tpm-matrices.', epilog = 'By Mutwil Lab')
    parser.add_argument('-w', '--workers', metavar='workers',
                        help='Enter the number of concurrent worker processes to be used.',
                        dest='workers', type=int, required=True)
    parser.add_argument('-n', '--newonly', action='store_true', default=False,
                        help='Include this optional tag if you only want to process species whose matrices does not yet exist.',
                        dest='newonly', required=False)
    args = parser.parse_args()
    workers = args.workers
    newonly = args.newonly

taxids = run_job.get_valid_jobs()
if newonly:
    files = os.listdir(f"{DATA_PATH}/postprocess/tpm-matrices")
    existing_runids = [int(re.findall("taxid(\d{4,6})", file)[0]) for file in files]
    taxids = list(set(taxids) - set(existing_runids))

if workers == 0:
    for taxid in taxids:
        pull_tpm_matrix.write_tpm_matrix(taxid)
elif workers > 0:
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(pull_tpm_matrix.write_tpm_matrix, taxid) for taxid in taxids]
        results = []
        for f in concurrent.futures.as_completed(futures):
            results.append(f.result())
else:
    raise Exception("Number of workers invalid")
