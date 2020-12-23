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
import concurrent.futures
import warnings
import pdb
# Relative imports
from config.constants import DATA_PATH
from download import download_functions, helpers
from preprocess.despatch import fetch_runtable

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'This script gets all CDS and runtables for all species.', epilog = 'By Mutwil Lab')
    parser.add_argument('-t', '--threads', metavar='threads',
                        help='Enter the number of concurrent threads to be used.',
                        dest='threads', type=int, required=True)
    parser.add_argument('-n', '--newonly', action='store_true', default=False,
                        help='Include this optional tag if you only want to download species whose runtables do not yet exist.',
                        dest='newonly', required=False)
    args = parser.parse_args()
    threads = args.threads
    newonly = args.newonly

# Go through latest job-list
job_list_dir = f"{DATA_PATH}/preprocess/job-list/"
latest_job_file = max(os.listdir(job_list_dir))
raw_job_df = pd.read_csv(f"{job_list_dir}{latest_job_file}", sep='\t')
job_df = raw_job_df.loc[raw_job_df['cds_link'].notnull(),:]

def process_cds(taxid, cds_link):
    idx_path = f"{DATA_PATH}/download/idx/taxid{taxid}.idx"
    cds_path = f"{DATA_PATH}/download/cds/taxid{taxid}.cds.fasta"
    if os.path.exists(idx_path):
        print(f"😊 kallisto index already exists for taxid{taxid}")
        return 1, idx_path
    elif not os.path.exists(cds_path):
        wget.download(cds_link, f"{DATA_PATH}/download/cds/taxid{taxid}.cds.fasta")
        if os.path.exists(cds_path):
            print(f"🙌🏻 Downloaded CDS for taxid{taxid}")
        else:
            print(f"😰 Failed to download CDS for taxid{taxid}")
            return -1, idx_path
    runtime, exit_code, _ = download_functions.kallisto_index(idx_path=idx_path, cds_path=cds_path)
    os.remove(cds_path)
    print(f"✓ Created kallisto index for taxid{taxid}, removed CDS")
    return 0, idx_path

def download_job(taxid, cds_link):
    status, idx_path = process_cds(taxid, cds_link)
    if newonly:
        # Redownload runtable only if it doesnt exist for now
        if not helpers.latest_runtable_path(f"taxid{taxid}"):
            try:
                fetch_runtable(taxid=taxid, db='sra')
                # fetch_runtable(taxid=taxid, db='ena')
            except:
                warnings.warn("Failed to download runtable.")
    else:
        try:
            fetch_runtable(taxid=taxid, db='sra')
            # fetch_runtable(taxid=taxid, db='ena')
        except:
            warnings.warn("Failed to download runtable.")
    return status

if threads == 0:
    for row in job_df.itertuples():
        _ = download_job(row.taxid, row.cds_link)
elif threads > 0:
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(download_job, row.taxid, row.cds_link) for row in job_df.itertuples()]
        results = []
        for f in concurrent.futures.as_completed(futures):
            results.append(f.result())
else:
    raise Exception("Number of threads invalid")

taxids = set(job_df['taxid'].tolist())

files = os.listdir(f"{DATA_PATH}/download/idx")
idxs = [int(re.findall("taxid(\d{4,6})", file)[0]) for file in files]
idxs = set(idxs)

files = os.listdir(f"{DATA_PATH}/preprocess/sra-runtables")
runtables = [int(re.findall("taxid(\d{4,6})", file)[0]) for file in files]
runtables = set(runtables)

print(f"Succesfully have {len(idxs)}/{len(taxids)} species with kallisto index processed.")
print(f"Successfully have {len(runtables)}/{len(taxids)} species with runtables downlaoded at least once.")
