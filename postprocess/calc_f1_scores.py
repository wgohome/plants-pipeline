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
from postprocess import coexpression1

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'This script returns the f1 scores for each species.', epilog = 'By Mutwil Lab')
    parser.add_argument('-t', '--taxids', metavar='taxids',
                        help='Enter the taxid of the species for which the f1 scores are to be calculated. For instance, Arabidopsis thaliana\'s taxid would be \'3702\'. A list of taxids can be accepted separated by space.',
                        dest='taxids', type=int, required=True, nargs='+')
    args = parser.parse_args()
    taxids = args.taxids

# t3218a, t3197, t3702, t3218b = os.listdir(f"{DATA_PATH}postprocess/percentage-matrices")
# df3218 = pd.read_csv(f"{DATA_PATH}postprocess/percentage-matrices/{t3218b}", sep='\t', index_col=0)
# df3197 = pd.read_csv(f"{DATA_PATH}postprocess/percentage-matrices/{t3197}", sep='\t', index_col=0)
# df3702 = pd.read_csv(f"{DATA_PATH}postprocess/percentage-matrices/{t3702}", sep='\t', index_col=0)
# taxid = 3702
# df = df3702

def read_percentage_matrix(taxid):
    files = sorted([file for file in os.listdir(f"{DATA_PATH}postprocess/percentage-matrices") if f"taxid{taxid}" in file])
    if files == []:
        return pd.DataFrame()
    df = pd.read_csv(f"{DATA_PATH}postprocess/percentage-matrices/{files[-1]}", sep='\t', index_col=0)
    df.index = df.index.str.upper()
    return df

def get_f1_stats(df, ribosomal_genes):
    df['ribosomal'] = df.index.isin(ribosomal_genes)
    f1_stats = {}
    for PCC_cutoff in df.columns:
        if PCC_cutoff == 'ribosomal':
            continue
        # f1_stats[PCC_cutoff] = {}
        for p_cutoff in [x/10 for x in range(1,10)]:
            tp = ((df['ribosomal'] == True) & (df[PCC_cutoff] >= p_cutoff)).sum()
            fp = ((df['ribosomal'] == False) & (df[PCC_cutoff] >= p_cutoff)).sum()
            fn = ((df['ribosomal'] == True) & (df[PCC_cutoff] < p_cutoff)).sum()
            tn = ((df['ribosomal'] == False) & (df[PCC_cutoff] < p_cutoff)).sum()
            precision = tp/(tp+fp)
            recall = tp/(tp+fn)
            f1 = (2*precision*recall)/(precision+recall)
            f1_stats[(PCC_cutoff, p_cutoff)] = {
                'tp': int(tp),
                'fp': int(fp),
                'fn': int(fn),
                'tn': int(tn),
                'precision': 0 if np.isnan(precision) else float(precision),
                'recall': 0 if np.isnan(recall) else float(recall),
                'f1': 0 if np.isnan(f1) else float(f1)
            }
    # print(json.dumps(f1_stats, indent=4))
    return f1_stats

def write_f1_stats(f1_stats, taxid, to_write=False):
    stats_df = pd.DataFrame(f1_stats).T
    stats_df.index.names = ('pcc_cutoff', 'per_cutoff')
    if to_write:
        stats_df.to_csv(f"{DATA_PATH}postprocess/f1-stats/{helpers.get_timestamp()}-taxid{taxid}_f1_stats.txt", sep='\t')
    return stats_df

def get_species_stats(taxid, bincodes=["17.1.2.1", "17.1.3.1"]):
    df = read_percentage_matrix(taxid)
    ribosomal_genes = coexpression1.get_genes_set(taxid, bincodes=bincodes)
    if df.empty:
        print(f"taxid{taxid} does not have a percentage matrix yet!")
        return None
    f1_stats = get_f1_stats(df, ribosomal_genes)
    stats_df = write_f1_stats(f1_stats, taxid, to_write=True)

if __name__ == '__main__':
    for taxid in taxids:
        get_species_stats(taxid, bincodes=["17.1.2.1", "17.1.3.1"])
