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
import numpy as np
import pandas as pd
import json
import pdb
# Relative imports of CONSTANTS in config/constants.py
from config.constants import DATA_PATH
from download import helpers, checkfiles

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'This script calculates the pcc between every genes for a species from all filtered Run IDs, and creates a pcc matrix for it in data-pipeline/postprocess/pcc-matrices.', epilog = 'By Mutwil Lab')
    parser.add_argument('-t', '--taxids', metavar='taxids',
                        help='Enter the taxid of the species for which the PCC matrix is to be created. For instance, Arabidopsis thaliana\'s taxid would be \'3702\'. A list of taxids can be accepted separated by space.',
                        dest='taxids', type=int, required=True, nargs='+')
    args = parser.parse_args()
    taxids = args.taxids

def calc_species(path):
    df = pd.read_csv(path, sep='\t', index_col=0, header=0)
    # Filter from qc-out
    filt_genes = get_filtered_genes(taxid)
    df = df.loc[:, filt_genes]
    genes = df.index
    npdata = df.to_numpy().astype('float64')
    npdata = np.nan_to_num(npdata)
    gaps = npdata - npdata.mean(axis=-1).reshape(-1, 1)
    gaps_sq = ((gaps ** 2).sum(axis = -1))**0.5
    return gaps, gaps_sq, genes

def get_filtered_genes(taxid):
    path = latest_qc_out()
    with open(path, 'r') as f:
        qc_out = json.load(f)
    return qc_out.get(str(taxid))

def pcc_neighbours(gaps, gaps_sq, x):
    pcc_vec = np.dot(gaps, gaps[x])/(gaps_sq[x] * gaps_sq)
    return pcc_vec

def log(path, to_write):
    with open(path, 'a') as f:
        f.write(to_write)

def latest_tpm_matrix(taxid):
    files = sorted([file for file in os.listdir(f"{DATA_PATH}postprocess/tpm-matrices") if f"taxid{taxid}" in file])
    if files == []:
        return None
    else:
        return f"{DATA_PATH}postprocess/tpm-matrices/{files[-1]}"

def latest_qc_out():
    files = sorted([file for file in os.listdir(f"{DATA_PATH}postprocess/qc-out")])
    if files == []:
        return None
    else:
        return f"{DATA_PATH}postprocess/qc-out/{files[-1]}"

def process_species(taxid):
    tpm_path = latest_tpm_matrix(taxid)
    if not tpm_path:
        print(f"taxid{taxid} does not have a tpm matrix yet!")
        return None
    out_path = f"{DATA_PATH}postprocess/pcc-matrices/{helpers.get_timestamp()}-taxid{taxid}.pcc.txt"
    print(f"Calculating for taxid{taxid} ...")
    gaps, gaps_sq, genes = calc_species(tpm_path)
    # write headers
    to_write = '\t' + '\t'.join(genes.tolist()) + '\n'
    log(out_path, to_write)
    for x in range(gaps.shape[0]):
        pcc_vec = np.dot(gaps, gaps[x])/(gaps_sq[x] * gaps_sq)
        pcc_vec = np.nan_to_num(pcc_vec)
        # write row
        pccs = [str(pcc) for pcc in pcc_vec]
        to_write = genes[x] + '\t' + '\t'.join(pccs) + '\n'
        log(out_path, to_write)
    print(f"Completed for taxid{taxid}!")
    return None

if __name__ == '__main__':
    for taxid in taxids:
        process_species(taxid)
