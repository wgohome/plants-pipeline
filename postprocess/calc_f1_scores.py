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
import warnings
import pdb
# Relative imports of CONSTANTS in config/constants.py
from config.constants import DATA_PATH
from download import helpers, checkfiles

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'This script finds the PCC between each gene and every other genes, tests with PCC threshold of 0.1 to 0.9, find the percentage of ribosomal genes among the filtered neighbours. This is returned in a matrix of percentages.', epilog = 'By Mutwil Lab')
    parser.add_argument('-t', '--taxids', metavar='taxids',
                        help='Enter the taxid of the species for which the PCC matrix is to be created. For instance, Arabidopsis thaliana\'s taxid would be \'3702\'. A list of taxids can be accepted separated by space.',
                        dest='taxids', type=int, required=True, nargs='+')
    args = parser.parse_args()
    taxids = args.taxids

def calc_species(path):
    df = pd.read_csv(path, sep='\t', index_col=0, header=0)
    # Filter from qc-out
    filt_runids = get_filtered_runids(taxid)
    # filt_runids must exist in df, is any key doesn't, KeyError will be raised
    runids = set(df.columns) & set(filt_runids)
    if (set(filt_runids) - runids) != set():
        print(f"Not in tpm table: {', '.join((set(filt_runids) - runids))}")
    df = df.loc[:, runids]
    genes = df.index
    npdata = df.to_numpy().astype('float64')
    npdata = np.nan_to_num(npdata)
    gaps = npdata - npdata.mean(axis=-1).reshape(-1, 1)
    gaps_sq = ((gaps ** 2).sum(axis = -1))**0.5
    return gaps, gaps_sq, genes

def get_filtered_runids(taxid):
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

def get_genes_set(taxid, bincodes=["17.1.2.1", "17.1.3.1"]):
    path = f"{DATA_PATH}postprocess/gene-classifications/taxid{taxid}_mercator.txt"
    if not os.path.exists(path):
        print(f"taxid{taxid} gene annotations is not found in pipeline-data/postprocess/gene-classifications/. Make sure it is labelled as taxidXXXX_mercator.txt")
        return []
    df = pd.read_csv(path, sep='\t')
    df['BINCODE'] = df['BINCODE'].str.strip("'")
    df['IDENTIFIER'] = df['IDENTIFIER'].str.strip("'")
    submasks = [df['BINCODE'].str.startswith(bincode) for bincode in bincodes]
    mask = pd.concat(submasks, axis=1).sum(axis=1) > 0
    ribosomal_series = df[mask]['IDENTIFIER']
    ribosomal_genes = ribosomal_series[ribosomal_series != ''].str.upper().tolist()
    return ribosomal_genes

def process_species(taxid, bincodes=["17.1.2.1", "17.1.3.1"]):
    tpm_path = latest_tpm_matrix(taxid)
    if not tpm_path:
        print(f"taxid{taxid} does not have a tpm matrix yet!")
        return None
    print(f"Calculating for taxid{taxid} ...")
    # List of ribosomal genes
    ribosomal_genes = get_genes_set(taxid, bincodes=bincodes)
    if ribosomal_genes == []:
        warnings.warn("Check if gene annotations are available or correct in the pipeline-data/postprocess/gene-classifications directory!")
        return None
    # Calculate PCC components
    gaps, gaps_sq, genes = calc_species(tpm_path)
    percentages = {}
    # for x in range(genes.size):
    for x in range(10):
        percentages[genes[x]] = {}
        pcc_vec = np.dot(gaps, gaps[x])/(gaps_sq[x] * gaps_sq)
        pcc_vec = np.nan_to_num(pcc_vec)
        for pcc_cutoff in [i/10 for i in range(1,10)]:
            neighbors = genes[pcc_vec >= pcc_cutoff]
            neighbors = neighbors.str.upper()
            ribo_neighbors = set(neighbors) & set(ribosomal_genes)
            if neighbors.empty:
                percentages[genes[x]][pcc_cutoff] = 0
            else:
                percentages[genes[x]][pcc_cutoff] = len(ribo_neighbors)/len(neighbors)
    df = pd.DataFrame(percentages).T
    df.index = df.index.str.upper()
    df.to_csv(f"{DATA_PATH}postprocess/percentage-matrices/{helpers.get_timestamp()}-taxid{taxid}_percentages.txt", sep='\t')
    print(f"Completed for taxid{taxid}!")
    return None

if __name__ == '__main__':
    for taxid in taxids:
        process_species(taxid, ["17.1.2.1", "17.1.3.1"])

# "17.1.2.1": 'Protein biosynthesis.ribosome biogenesis.large ribosomal subunit (LSU).LSU proteome'
# "17.1.3.1": 'Protein biosynthesis.ribosome biogenesis.small ribosomal subunit (SSU).SSU proteome'
