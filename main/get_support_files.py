# Module Header ################################################################
import sys
from os.path import realpath, dirname
import re

if __name__ == '__main__':
    abspath = realpath(dirname(__file__))
    parent_module = re.search('^(.*plants-pipeline)', abspath).group()
    sys.path.insert(0, parent_module)
################################################################################

import numpy as np
import pandas as pd
import os
import wget
import pdb
# Relative imports
from config.constants import DATA_PATH
from download import download_functions
from preprocess.despatch import fetch_runtable

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

for row in job_df.itertuples():
    status, idx_path = process_cds(row.taxid, row.cds_link)
    # Redownload runtable anyway to get the latest
    fetch_runtable(taxid=row.taxid, db='sra')
    fetch_runtable(taxid=row.taxid, db='ena')
