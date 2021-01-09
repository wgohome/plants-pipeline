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
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import json
import pdb
# Relative imports of CONSTANTS in config/constants.py
from config.constants import DATA_PATH
from download import helpers, checkfiles

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'This script generates qc statistics from all experiments of given speices.', epilog = 'By Mutwil Lab')
    parser.add_argument('-t', '--taxids', metavar='taxids',
                        help='Enter the taxids of the species for which the TPM matrix is to be created. For instance, Arabidopsis thaliana\'s taxid would be \'3702\'. A list of taxids can be accepted separated by space.',
                        dest='taxids', type=int, required=False, nargs='+')
    parser.add_argument('-u', '--update', action='store_true', default=False,
                        help='Include this optional tag if you only want to update the QC cutoffs based on the latest-timestamped qc-summary table you editted.',
                        dest='update', required=False)
    args = parser.parse_args()
    taxids = args.taxids
    update = args.update

###########################################################################

def get_qc_stats(taxid):
    df = read_runinfo(taxid)
    if df.empty:
        return None, None, None
    qc_df = get_qc_matrix(df)
    df, logn_cut, p_cut = set_threshold(df, qc_df)
    qc_df.to_csv(f"{DATA_PATH}postprocess/qc-matrices/{helpers.get_timestamp()}-taxid{taxid}_qc_matrix.txt", sep='\t')
    plot_qc(df, taxid)
    return logn_cut, p_cut, df['pass'].sum()

def update_qc_stats(taxid, logn_cut, p_cut):
    df = read_runinfo(taxid)
    if df.empty:
        return None, None, None
    df, logn_cut, p_cut = set_threshold(df, logn_cut=logn_cut, p_cut=p_cut)
    # Did not write qc-matrices
    plot_qc(df, taxid)
    return logn_cut, p_cut, df

def read_runinfo(taxid):
    path = f"{DATA_PATH}download/runinfo-main/taxid{taxid}_runinfo_main.txt"
    if not os.path.exists(path):
        print(f"{path} does not exist!")
        return pd.DataFrame()
    df = pd.read_csv(path, sep='\t', index_col='runid', usecols=['n_processed', 'p_pseudoaligned', 'runid'])
    df['log10(processed)'] = df['n_processed'].apply(math.log10)
    return df

def get_qc_matrix(df):
    max_logn = np.arange(6,8.2,0.2) # rows
    max_p = [p for p in range(10,100,10)] # columns
    qc_matrix = np.full(len(max_logn)*len(max_p), 0).reshape(len(max_logn),len(max_p))
    for i, p in enumerate(max_p): # columns
        for j, l in enumerate(max_logn): # rows
            qc_matrix[j, i] = ((df['p_pseudoaligned'] > p) & (df['log10(processed)'] > l)).sum()
    qc_df = pd.DataFrame(qc_matrix, index=max_logn, columns=max_p)
    qc_df.index.name = "log_processed|p_pseudoaligned"
    return qc_df

def set_threshold(df, qc_df=None, logn_cut=None, p_cut=None):
    if logn_cut == None or p_cut == None:
        logn_cut, p_cut = suggest_threshold(df, qc_df)
    df['pass'] = (df['p_pseudoaligned'] > p_cut) & (df['log10(processed)'] > logn_cut)
    df['pass'] = df['pass'].astype(int)
    return df, logn_cut, p_cut

def suggest_threshold(df, qc_df):
    # Make sure it is such that >=1000 experiments or >=60% of all experiments accoutned for
    row = qc_df.index.get_loc(7.0, method='nearest')
    col = qc_df.columns.get_loc(70, method='nearest')
    n = qc_df.iloc[row,col]
    if n >= 1000:
        return 7.0, 70
    elif (n/df.shape[0]) >= 0.6:
        return 7.0, 70
    else:
        minimum = df.shape[0]*0.5
        vals = qc_df[qc_df >= minimum].values.flatten()
        # Find biggest treshold to hit 60% of all samples
        try:
            n = vals[~np.isnan(vals)].min()
        except:
            return 6.5, 0.5
        n_position = qc_df.where(qc_df == n).dropna(axis=0, how='all').dropna(axis=1, how='all')
        # row = qc_df[(qc_df == n).sum(axis=1) == 1].where(qc_df == n).index[0]
        row = n_position.index[0]
        col = n_position.columns[0]
        return row, col

def plot_qc(df, taxid):
    sns.set(style="darkgrid")
    g = sns.jointplot('p_pseudoaligned', 'log10(processed)', data = df, kind="reg", xlim=(0, 100), color='b', scatter=False,
                      height=10, marginal_kws=dict(bins=20, rug=True))
    g.ax_joint.cla()
    # c_dic = {1: 'b', 0: 'r'}
    # g.ax_joint.scatter(df['p_pseudoaligned'], df['log10(processed)'], c=df['pass'].apply(lambda x: c_dic[x]), label=['passed', 'failed'], marker='x',)
    sc0 = g.ax_joint.scatter(df[df['pass'] == 0]['p_pseudoaligned'], df[df['pass'] == 0]['log10(processed)'], c='r', label='failed', marker='+', alpha=0.4)
    sc1 = g.ax_joint.scatter(df[df['pass'] == 1]['p_pseudoaligned'], df[df['pass'] == 1]['log10(processed)'], c='b', label='passed', marker='+', alpha=0.4)
    g.ax_marg_x.set_title(f"taxid{taxid}", fontsize=16, pad=20)
    g.ax_joint.set_xlabel("% reads pseudoaligned", fontsize=12)
    g.ax_joint.set_ylabel("log10(reads processed)", fontsize=12)
    g.ax_joint.set(xlim=(-4,100), ylim=(0,9))
    # g.fig.set_size_inches(12,10)
    g.ax_joint.legend(loc="lower center")
    plt.tight_layout()
    g.fig.savefig(f"{DATA_PATH}postprocess/qc-jointplots/taxid{taxid}_qc_jointplot.png")

def get_cutoffs(taxids):
    cutoffs = {}
    for taxid in taxids:
        print(f"Processing QC for taxid{taxid} ...")
        logn_cut, p_cut, n = get_qc_stats(taxid) # writes qc-matrices, plots qc-jointplots for taxid
        cutoffs[taxid] = {
            'logn_cut': logn_cut,
            'p_cut': p_cut,
            'n': n
        }
    cutoff_df = pd.DataFrame(cutoffs).T
    cutoff_df = cutoff_df[['logn_cut', 'p_cut', 'n']] # Make sure columns in the right order to write
    # writes qc summary
    cutoff_df.to_csv(f"{DATA_PATH}postprocess/qc-summary/{helpers.get_timestamp()}-cutoffs.txt", sep='\t', header=['log10(processed) cutoff', '% pseudoaligned cutoff', 'Number of Run IDs passed'])
    return None

def update_cutoffs():
    file = sorted(os.listdir(f"{DATA_PATH}postprocess/qc-summary/"))
    if file == []:
        print(f"Does not have a qc-matrix yet!")
        return None
    path = f"{DATA_PATH}postprocess/qc-summary/{file[-1]}"
    summary_df = pd.read_csv(path, sep='\t', header=0, names=['logn_cut', 'p_cut', 'n'], index_col=0)
    filtered_runids = {}
    for row in summary_df.itertuples():
        logn_cut, p_cut, df = update_qc_stats(row.Index, row.logn_cut, row.p_cut) # plots qc-jointplot
        summary_df.loc[row.Index, 'n'] = df['pass'].sum()
        filtered_runids[row.Index] = df[df['pass'] == 1].index.tolist()
    # Update qc-summary
    summary_df.to_csv(f"{DATA_PATH}postprocess/qc-summary/{helpers.get_timestamp()}-cutoffs.txt", sep='\t', header=['log10(processed) cutoff', '% pseudoaligned cutoff', 'Number of Run IDs passed'])
    # write to qc-out
    with open(f"{DATA_PATH}postprocess/qc-out/{helpers.get_timestamp()}-qc_out.json", 'w') as f:
        json.dump(filtered_runids, f)
    os.remove(path)
    return None

if __name__ == "__main__":
    if update:
        # if not taxids:
        #     print("Do specify taxids to process!")
        # else:
        update_cutoffs()
    else:
        if not taxids:
            print("Do specify taxids to process!")
        else:
            get_cutoffs(taxids)
