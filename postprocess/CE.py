'''
Loose approach: when a gene is assigned > 1 bin, consider both bins
Network of functional annotations
-> To find Mapman processes that are enriched in each HCCA cluster, connect these processes.
-> Generate a network of Mapman processes.

Legend:
fn: function; ls: list, dc: dictionary, tp: tuple, set: set,
ar: np array, sr: pd series, df: pd dataframe
'''

import argparse
import numpy as np
import pandas as pd
import random # To seed for testing only
from random import shuffle
from statsmodels.stats.multitest import multipletests
import pdb

################################################################################
# # Need to add argv command line interface

# # Parse arguments from command line interface
# parser = argparse.ArgumentParser(description = 'This script generates networks of MapmanBins that are over or under represented within their clusters',
#                                  epilog = 'By Mutwil Lab')
# parser.add_argument('-c', '--clusters', nargs=1, metavar='hcca_clusters_filepath',
#                     help='Path for a HCCA clusters file to be processed',
#                     dest='hcca_clusters_filepath', type=str, required=True)
# parser.add_argument('-m', '--mapman', nargs=1, metavar='mapman_filepath',
#                     help='Path for Mapman Bins assignment file from Mercator',
#                     dest='mapman_filepath', type=str, required=True)
# parser.add_argument('-t', '--tpm', nargs=1, metavar='tpm_filepath',
#                     help='Expression matix file, to get gene indexes',
#                     dest='tpm_filepath', type=str, required=True)
# parser.add_argument('-o', '--overrep', nargs=1, metavar='overrep_bins_filepath',
#                     help='Name and relative path of the expected output file for over-represented Bins',
#                     dest='overrep_bins_filepath', type=str, required=True)
# parser.add_argument('-u', '--underrep', nargs=1, metavar='underrep_bins_filepath',
#                     help='Name and relative path of the expected output file for under-represented Bins',
#                     dest='underrep_bins_filepath', type=str, required=True)
# parser.add_argument('-a', '--alpha', nargs=1, metavar='alpha',
#                     help='Familywise Error Rate for Benjamini-Hochberg Procedure in each cluster',
#                     dest='alpha', type=float, required=True)
# args = parser.parse_args()

# mapman_filepath = args.hcca_clusters_filepath[0]
# clusters_filepath = args.mapman_filepath[0]
# tpm_filepath = args.tpm_filepath[0]
# overrep_bins_filepath = args.overrep_bins_filepath[0]
# underrep_bins_filepath = args.underrep_bins_filepath[0]
# FDR = args.alpha[0]

# species = "Cpa"
# irt = f"../plants/{species}/"
# ort = f"CE_test/more_tests/{species}/"
# mapman_path = irt + f"{species}.mercator.txt"
# hcca_path = irt + f"{species}.hcca"
# tpm_path = irt + f"{species}_matrix.txt"
# overrep_path = ort + f"{species}_overrep_bins" + "_{M}.txt"
# underrep_path = ort + f"{species}_underrep_bins" + "_{M}.txt"
# FDR = 0.05
# METHOD = 2

################################################################################

def main():
    # random.seed(0)
    ref_df, master_df, clus_dc, samples_ar =\
        process_input(hcca_path, mapman_path, tpm_path, level=2)
    all_or_pairs, all_ur_pairs = get_all_bin_pairs(master_df, clus_dc, samples_ar)
    write_results(all_or_pairs, overrep_path.format(M = METHOD), ref_df)
    write_results(all_ur_pairs, underrep_path.format(M = METHOD), ref_df)

def process_input(hcca_path, mapman_path, tpm_path, level=2):
    mapman_df = read_mapman(mapman_path)
    genes_df = read_tpm(tpm_path)
    ref_df = merge_df(mapman_df, genes_df, level)
    master_df = ref_df.copy()
    if METHOD == 2: # For conservative method, collapse bins by gene
        master_df = collapse_by_gene(ref_df)
    clus_dc = read_hcca(hcca_path)
    samples_ar = generate_samples(master_df) # arrays of lists in METHOD 3
    return ref_df, master_df, clus_dc, samples_ar

def read_mapman(mapman_path):
    mapman = pd.read_csv(mapman_path, sep='\t', header=0,\
                             index_col=False, quotechar="'")
    mapman.dropna(axis=0, subset=['IDENTIFIER'] , inplace=True) # Remove bins with no associated gene IDs
    mapman['IDENTIFIER'] = mapman['IDENTIFIER'].str.lower()
    return mapman

def read_tpm(tpm_path):
    genes = pd.read_csv(tpm_path, sep='\t', header=0, index_col=False,\
                        usecols=[0])
    genes['genes_index'] = genes.index
    genes['GeneID'] = genes['GeneID'].str.lower()
    return genes

def merge_df(mapman, genes, level):
    mapman = mapman[mapman['IDENTIFIER'].isin(genes['GeneID'])] # Exclude gene IDs in Mapman that is not in tpm
    merged = genes.merge(mapman, left_on='GeneID', right_on='IDENTIFIER')\
                  .drop('GeneID', axis='columns')
    merged['BINCODE'] = merged['BINCODE'].apply(serialise) # Serialise all BINCODES to allow sorting and hashing
    merged['BINCODE'] = merged['BINCODE'].apply(set_level, args=(level,)) # Truncate bincode to the desired level
    merged['NAME'] = merged['NAME'].apply(set_level, args=(level,))
    return merged

def collapse_by_gene(master_df):
    """ If a gene has more than one bin, collapse the rows of the gene and
    represent the BINCODE as tuples and NAME/DESCRIPTION as lists """
    collapsed_df = master_df.groupby('genes_index')\
                         .aggregate({'IDENTIFIER':max, 'BINCODE':list,
                                     'NAME':list, 'DESCRIPTION':list,
                                     'genes_index':'first'}) # genes_index is not the index
    return collapsed_df

def serialise(bincode):
    """ Assumes each level is at most 2 digits """
    serialised = []
    for l in bincode.split('.'):
        if int(l) < 10:
            l = '0' + l
        serialised.append(l)
    return '.'.join(serialised)

def set_level(bincode, level):
    return '.'.join(bincode.split('.')[:level])

def read_hcca(hcca_path):
    clus = pd.read_csv(hcca_path, sep='\t', names=['gene_id', 'cluster_id'],
                           index_col=0)
    valid_clus = clus[clus['cluster_id'].str.isnumeric()].astype('int32')
    clus_dc = valid_clus.groupby('cluster_id').groups
    return clus_dc

def generate_samples(master_df):
    """ Return list samples of 1000 inner lists of:
    [Method 1] all l2 bins of every gene in tpm, in the corresponsing frequency
    [Method 2] if genes > 1 bin, bins combined into a hashable tuple"""
    all_bins = master_df['BINCODE']
    samples = np.repeat([all_bins], 1000, axis = 0)
    np.array(list(map(np.random.shuffle, samples))) # Shuffles samples in place
    return samples

def get_all_bin_pairs(master_df, clus_dc, samples_ar):
    all_or_pairs, all_ur_pairs = [], []
    for clus_id, clus_genes in clus_dc.items():
        clus_df = master_df[master_df['genes_index'].isin(clus_genes)]
        clus_samp_ar = samples_ar[:, :(clus_df.shape[0])] # Sample based on number of genes in cluster
        or_pairs, ur_pairs = clus_bin_pairs(clus_df, clus_samp_ar)
        all_or_pairs.extend(or_pairs)
        all_ur_pairs.extend(ur_pairs)
    all_or_pairs = list(set(all_or_pairs)) #Reduce redundancy
    all_ur_pairs = list(set(all_ur_pairs))
    return all_or_pairs, all_ur_pairs

def clus_bin_pairs(clus_df, clus_samp_ar):
    clus_df = clus_df.explode('BINCODE') # Need to unpack across axis=-1, list becomes strings
    clus_count = clus_df.groupby('BINCODE').aggregate({'IDENTIFIER':'count'})
    clus_bins = clus_count.index.to_numpy().reshape(-1,1,1)
    samp_count = ar_count(clus_bins, clus_samp_ar).sum(axis=-1)
    or_pairs = get_sig_pairs(clus_count, samp_count, goal="overrep")
    ur_pairs = get_sig_pairs(clus_count, samp_count, goal="underrep")
    return or_pairs, ur_pairs

def ar_count(clus_bins, clus_samp_ar):
    def count(clus_bins, clus_samp_ar):
        return True if clus_bins in clus_samp_ar else False
    vcount = np.vectorize(count)
    return vcount(clus_bins, clus_samp_ar)

def get_sig_pairs(clus_count, samp_count, goal="overrep"):
    if goal == "underrep":
        p_ar = ((samp_count <= clus_count.to_numpy()).sum(axis=-1)+1)/1000
    else:
        p_ar = ((samp_count >= clus_count.to_numpy()).sum(axis=-1)+1)/1000
    mask_bool_ar, corrected_p, _, _ = multipletests(p_ar, alpha=FDR,
        method="fdr_bh", is_sorted=False, returnsorted=False)
    sig_df = clus_count[mask_bool_ar]
    sig_pairs = pair_perm(list(sig_df.index))
    return sig_pairs

def pair_perm(sig_ls):
    pairs = [] # List of lists of pairs
    for i in range(len(sig_ls) - 1):
        for j in range(i + 1, len(sig_ls)):
            pair_ls = sorted([sig_ls[i], sig_ls[j]]) # Sort order within pair
            pair_tp = tuple(pair_ls) # Hashable
            pairs.append(pair_tp)
    return pairs

def write_results(all_pairs, path, master_df):
    to_write = '\t'.join(["Bin 1", "Bin 1 Name", "Bin 2", "Bin 2 Name"]) + "\n"
    for pair in all_pairs:
        name1, name2 = get_bin_names(pair, master_df)
        entries = list(map(str, [pair[0], name1, pair[1], name2]))
        to_write += '\t'.join(entries) + '\n'
    with open(path, 'w') as f:
        f.writelines(to_write)
    print("Created {path}".format(path=path))

def get_bin_names(pair, master_df):
    names = []
    for i in range(2):
        names.append(master_df[master_df['BINCODE'] == pair[i]].iloc[0]['NAME'])
    return names[0], names[1]
