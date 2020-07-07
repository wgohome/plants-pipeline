"""
2020-03-25 14:00
HCCA Version 3 (with Numpy, matrix PCC calc of covariance)
-> Generate HRR of every gene in a species based on its TPM from
its expression matrix.
"""

import os
import argparse
import numpy as np
import pandas as pd
from time import time
import pdb
import timeit

################################################################################
# Parse arguments from CLI
def parse():
    parser = argparse.ArgumentParser(description = 'This script generates the Highest Reciprocal Rank (HRR) between every pair of genes in a given expression matrix.',
                                     epilog = 'By Mutwil Lab')
    parser.add_argument('-i', '--input', nargs=1, metavar='tpm_matrix_filepath',
                        help='Path for a tpm matrix file to be processed, eg ./expression_matrices/filename.txt',
                        dest='tpm_matrix_path', type=str, required=True)
    parser.add_argument('-o', '--output', nargs=1, metavar='hrr_matrix_filepath',
                        help='File path for a hrr matrix file to be created, eg ./hrr_matrices/filename.txt',
                        dest='hrr_matrix_path', type=str, required=True)
    parser.add_argument('-c', '--cutoff', nargs=1, metavar='rank_cutoff',
                        help='Cutoff for rank of gene by PCC of gene pair in consideration, default cutoff = 100',
                        dest='CUTOFF', default=[100], type=int, required=False)
    args = parser.parse_args()
    global tpm_matrix_path, hrr_matrix_path, CUTOFF
    tpm_matrix_path = args.tpm_matrix_path[0]
    hrr_matrix_path = args.hrr_matrix_path[0] # To be created
    CUTOFF = args.CUTOFF[0]

################################################################################
# Main

def main():
    top100_matrix = get_pcc_matrix(tpm_matrix_path)
    write_hrr_matrix(top100_matrix)

################################################################################
# Functions

def get_pcc_matrix(tpm_matrix_path):
    """ Return top100_matrix of indexes of top100 neighbours for each gene """
    pddata = pd.read_csv(tpm_matrix_path, sep='\t', header=0, index_col=0)
    npdata = pddata.to_numpy().astype('float64')
    gaps = npdata - npdata.mean(axis = -1).reshape(npdata.shape[0], 1)
    gaps_sq = np.sqrt((gaps ** 2).sum(axis = -1)) # .reshape(npdata.shape[0], 1)
    top100_matrix = np.arange(CUTOFF)
    for x in range(gaps.shape[0]):
        pcc_vector = np.dot(gaps, gaps[x])/(gaps_sq[x] * gaps_sq)
        top100 = np.argsort(-pcc_vector)[:CUTOFF] # nan will be sorted as last by default
        top100_matrix = np.vstack((top100_matrix, top100))
    return top100_matrix[1:]

def write_hrr_matrix(top100_matrix):
    """ Write line by line to file """
    with open(hrr_matrix_path, "w") as hrr_matrix:
        hrr_matrix.write("")
    for gene1_index in range(top100_matrix.shape[0]):
        line = f"{gene1_index}\t\t\t\t"
        for gene2_rank in range(CUTOFF):
            line = get_hrr(gene1_index, gene2_rank, line, top100_matrix)
        with open(hrr_matrix_path, "a+") as hrr_matrix:
            hrr_matrix.write(line + '\n')

def get_hrr(gene1_index, gene2_rank, line, top100_matrix):
    """ Return line, to which line is appended with 'gene2_index+HRR'
    if gene1 is ranked above CUTOFF amongst gene2's neighbours,
    after finding the HRR between the pair gene1 and gene2"""
    gene2_index = top100_matrix[gene1_index, gene2_rank]
    search_array = np.where(top100_matrix[gene2_index] == gene1_index)[0]
    if search_array.size == 0:
        hrr = -1
    else:
        gene1_rank = int(search_array[0])
        hrr = max(gene1_rank, gene2_rank) + 1
        # TO_EDIT: Need to +1 here to compensate for starting from 0?
        # Depends on whether we consider the gene itself in the ranking
        line += f"\t{gene2_index}+{hrr}"
        # Only add gene_index to the hrr_matrix if hrr > 0, i.e. both genes are reciprocally ranked lower number than CUTOFF
    return line

################################################################################

if __name__ == '__main__':
    parse()
    main()
