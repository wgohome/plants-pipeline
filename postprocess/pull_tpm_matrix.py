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
from zipfile import ZipFile
import pdb
# Relative imports of CONSTANTS in config/constants.py
from config.constants import DATA_PATH
from download import helpers, checkfiles

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'This script gets all downloaded experiments for a given species, and creates a tpm matrix for it in data-pipeline/post-process/tpm-matrices.', epilog = 'By Mutwil Lab')
    parser.add_argument('-t', '--taxid', metavar='taxid',
                        help='Enter the taxid of the species for which the TPM matrix is to be created. For instance, Arabidopsis thaliana\'s taxid would be \'3702\'.',
                        dest='taxid', type=int, required=True)
    args = parser.parse_args()
    taxid = args.taxid

def zipped_path(runid):
    return f"{DATA_PATH}/download/kallisto-out/{helpers.get_dirs(runid)}{runid}.zip"
def unzipped_path(runid):
    return f"{DATA_PATH}/download/kallisto-out/{helpers.get_dirs(runid)}{runid}"

def target_path(runid):
    return f"{DATA_PATH}/download/kallisto-out/{helpers.get_dirs(runid)}"

def get_tpm_matrix(taxid):
    completed_df, incomplete_df = checkfiles.validate_latest_batch(taxid, to_log=False) # Taxid runtable have to exist! Should not be a problem if calling from job_list
    if (completed_df.empty) and (incomplete_df.empty):
        print(f"taxid{taxid} is invalid. Check that its runtable is downloaded or if it is even a valid taxid.")
        return pd.DataFrame()
    series = []
    for runid in completed_df['Run']:
        # unzip kallisto file
        ZipFile(zipped_path(runid)).extractall(target_path(runid))
        # Get TPM values
        tpm_series = pd.read_csv(f"{unzipped_path(runid)}/abundance.tsv", sep='\t', index_col='target_id')['tpm']
        tpm_series.name = runid
        series.append(tpm_series)
        # Delete the unzipped kallisto file
        shutil.rmtree(unzipped_path(runid))
    if series:
        tpm_matrix = pd.concat(series, axis=1)
        return tpm_matrix
    else:
        return pd.DataFrame()

def tpm_matrix_name(taxid):
    return f"{DATA_PATH}/postprocess/tpm-matrices/{helpers.get_timestamp()}-taxid{taxid}_tpm_matrix.tsv"

def write_tpm_matrix(taxid):
    tpm_matrix = get_tpm_matrix(taxid)
    if tpm_matrix.empty:
        print(f"taxid{taxid} does not have any downloaded experiments.")
    else:
        tpm_matrix.to_csv(tpm_matrix_name(taxid), sep='\t')
        print(f"TPM matrix written/updated for taxid{taxid}")
    return tpm_matrix.shape

if __name__ == '__main__':
    write_tpm_matrix(taxid)

__all__ = ['write_tpm_matrix']
