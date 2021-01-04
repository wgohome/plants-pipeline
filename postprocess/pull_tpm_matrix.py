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
import gc
import fileinput
import pdb
# Relative imports of CONSTANTS in config/constants.py
from config.constants import DATA_PATH
from download import helpers, checkfiles

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'This script gets all downloaded experiments for a given species, and creates a tpm matrix for it in data-pipeline/post-process/tpm-matrices.', epilog = 'By Mutwil Lab')
    parser.add_argument('-t', '--taxid', metavar='taxid',
                        help='Enter the taxid of the species for which the TPM matrix is to be created. For instance, Arabidopsis thaliana\'s taxid would be \'3702\'.',
                        dest='taxid', type=int, required=True)
    parser.add_argument('-m', '--method', default=3,
                        help='Include this optional tag to indicate which method should be used to extract the TPM values. By default, method 3 will be used, which is the least RAM intensive',
                        dest='method', type=int, required=False)
    args = parser.parse_args()
    taxid = args.taxid
    method = args.method

# METHOD 1 #################################################################

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
        tpm_matrix.index = tpm_matrix.index.str.upper()
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

# END OF METHOD 1 ##########################################################

# METHOD 2 #################################################################
# Write into file line by line, cannot be parallelized
def chunk_runids(completed_df):
    total = completed_df['Run'].size
    n_full_chunks = total // 1000
    remainder_chunk = total % 1000
    print(f"{n_full_chunks} chunks of 1000s and {remainder_chunk} remainder RunIDs in the last chunk.")
    return n_full_chunks, remainder_chunk

def extract_tpm_values(runid):
    ZipFile(zipped_path(runid)).extractall(target_path(runid))
    tpm_series = pd.read_csv(f"{unzipped_path(runid)}/abundance.tsv", sep='\t', index_col='target_id')['tpm']
    tpm_series.name = runid
    shutil.rmtree(unzipped_path(runid))
    return tpm_series

def tpm_matrix_name2(taxid, chunk=None):
    if chunk:
        return f"{DATA_PATH}/postprocess/tmp-tpm/{helpers.get_timestamp()}-taxid{taxid}_tpm_chunk_{chunk}.tsv"
    else:
        return f"{DATA_PATH}/postprocess/tpm-matrices/{helpers.get_timestamp()}-taxid{taxid}_tpm_matrix.tsv"

def delete_files(paths):
    for path in paths:
        os.remove(path)

def write_tpm_matrix2(taxid):
    # Must be performed sequentially! Cannot be parallelized
    completed_df, incomplete_df = checkfiles.validate_latest_batch(taxid, to_log=False) # Taxid runtable have to exist! Should not be a problem if calling from job_list
    if (completed_df.empty) and (incomplete_df.empty):
        print(f"taxid{taxid} is invalid. Check that its runtable is downloaded or if it is even a valid taxid.")
        return 0
    print(f"\nExtracting TPM matrix for taxid{taxid} ...")
    # Chunk up the data in smaller files
    n_full_chunks, remainder_chunk = chunk_runids(completed_df)
    chunk_paths = []
    # Process full chunks of 1000, except the last remainder chunk
    for i in range(n_full_chunks):
        series = []
        for runid in completed_df['Run'][i:i+1000]:
            series.append(extract_tpm_values(runid))
        chunk_matrix = pd.concat(series, axis=1)
        filepath = tpm_matrix_name2(taxid, chunk=i+1)
        chunk_paths.append(filepath)
        chunk_matrix.to_csv(filepath, sep='\t')
        print(f"Written file for chunk {i+1}")
    # Process remainder chunk
    series = []
    for runid in completed_df['Run'][-remainder_chunk:]:
        series.append(extract_tpm_values(runid))
    chunk_matrix = pd.concat(series, axis=1)
    filepath = tpm_matrix_name2(taxid, chunk='r')
    chunk_paths.append(filepath)
    chunk_matrix.to_csv(filepath, sep='\t')
    print(f"Written file for remainder chunk")
    # Stitch all the chunks tgt
    dfs = []
    for path in chunk_paths:
        dfs.append(pd.read_csv(path, sep='\t', index_col='target_id'))
    master_df = pd.concat(dfs, axis=1, join='outer')
    delete_files(chunk_paths)
    del dfs, series, chunk_matrix
    gc.collect()
    master_df.to_csv(tpm_matrix_name2(taxid), sep='\t')
    print(f"Combined TPM matrix for taxid{taxid} written!")
    return 1

# END OF METHOD 2 ##########################################################

# METHOD 3 #################################################################

# Notes
# Avoid accumlating TPM values in RAM, write immediately
# Avoid having to transpose table (RAM-intesive for big tables), append to each line for every runid added instead.

# List of genes
# Initiate TPM Write a row of headers for genes
# For each runid
#   Create a list of TPM values matching the list of genes' indexes
#   If there are new genes not in header, append at the back and add the genes to the header
#   Write the tpm row into the file
def write_header(genes, filepath):
    to_write = "target_id\n" + '\n'.join(genes) + '\n'
    with open(filepath, 'w+') as f:
        f.write(to_write)

def update_header(new_genes, filepath):
    with open(filepath, 'r') as f:
        headers = f.readline().strip().split('\t')
    runs = len(headers) - 1
    dummy_tpms = ['0']*runs
    extension = '\t' + '\t'.join(dummy_tpms) + '\n'
    with open(filepath, 'a+') as f:
        for gene in new_genes:
            f.write(f"{gene}{extension}")

def write_row(tpms, runid, filepath):
    tpms = [str(tpm) for tpm in tpms]
    for line in fileinput.input(filepath, inplace=True):
        if fileinput.filelineno() == 1:
            header_extension = f"\t{runid}\n"
            print(line.replace('\n', header_extension), end='')
        else:
            tpm_extension = f"\t{tpms.pop(0)}\n"
            print(line.replace('\n', tpm_extension), end='')

def any_new_genes(master_genes, series):
    # Order TPM based on master_genes list
    # return list of new genes to be appended to the end of master_genes
    dic = series.to_dict()
    tpms = []
    for gene in master_genes:
        tpm = dic.get(gene)
        if tpm == None:
            tpms.append(0)
        else:
            tpms.append(tpm)
            dic.pop(gene)
    new_genes = []
    if dic: # not empty
        for gene, tpm in dic.items():
            new_genes.append(gene)
            tpms.append(tpm)
    return tpms, new_genes

def write_tpm_matrix3(taxid):
    completed_df, incomplete_df = checkfiles.validate_latest_batch(taxid, to_log=False) # Taxid runtable have to exist!
    if (completed_df.empty) and (incomplete_df.empty):
        print(f"taxid{taxid} is invalid. Check that its runtable is downloaded or if it is even a valid taxid.")
        return 0
    print(f"\nExtracting TPM matrix for taxid{taxid} ...")
    filepath = tpm_matrix_name(taxid)
    master_genes = None
    for runid in completed_df['Run']:
        series = extract_tpm_values(runid)
        if master_genes == None:
            master_genes = series.index.tolist()
            write_header(series.index, filepath)
            write_row(series.values, runid, filepath)
        else:
            tpms, new_genes = any_new_genes(master_genes, series)
            if new_genes: # There are genes not in master_genes
                update_header(new_genes, filepath)
                master_genes.extend(new_genes)
            write_row(tpms, runid, filepath)

# END OF METHOD 3 ##########################################################

if __name__ == '__main__':
    if method == 1:
        write_tpm_matrix(taxid)
    elif method == 2:
        write_tpm_matrix2(taxid)
    elif method == 3:
        write_tpm_matrix3(taxid)
    else:
        print("Method {method} is invalid.")

__all__ = ['write_tpm_matrix', 'write_tpm_matrix2', 'write_tpm_matrix3']
