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
    parser.add_argument('-t', '--threads', metavar='threads',
                        help='Enter the number of concurrent threads to be used.',
                        dest='threads', type=int, required=True)
    parser.add_argument('-n', '--newonly', action='store_true', default=False,
                        help='Include this optional tag if you only want to process species whose matrices does not yet exist.',
                        dest='newonly', required=False)
    args = parser.parse_args()
    threads = args.threads
    newonly = args.newonly

taxids = run_job.get_valid_jobs()
if newonly:
    files = os.listdir(f"{DATA_PATH}/postprocess/tpm-matrices")
    existing_taxids = [int(re.findall("taxid(\d{4,6})", file)[0]) for file in files]
    taxids = list(set(taxids) - set(existing_taxids))

if threads == 0:
    for taxid in taxids:
        pull_tpm_matrix.write_tpm_matrix(taxid)
elif threads > 0:
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(pull_tpm_matrix.write_tpm_matrix, taxid) for taxid in taxids]
        results = []
        for f in concurrent.futures.as_completed(futures):
            results.append(f.result())
else:
    raise Exception("Number of threads invalid")
