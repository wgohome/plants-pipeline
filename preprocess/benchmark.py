# Module Header ################################################################
import sys
from os.path import realpath, dirname
import re

if __name__ == '__main__':
    abspath = realpath(dirname(__file__))
    parent_module = re.search('^(.*plants-pipeline)', abspath).group()
    sys.path.insert(0, parent_module)

################################################################################

# Packages
import os
import pandas as pd
import numpy as np
import pdb
# Relative imports
from config.constants import DATA_PATH
from preprocess.po_parser.po import po # po is an Ontology class instance
from preprocess import iohelper
from download import helpers

annotation1_files = os.listdir(f"{DATA_PATH}preprocess/sra-annotations")
annotation2_files = os.listdir(f"{DATA_PATH}preprocess/sra-annotation2s")
evorepro_filepath = "/Users/wirriamm/OneDrive - Nanyang Technological University/marek_lab/annotation-benchmark/evorepro_s1.csv"
evorepro_df = pd.read_csv(evorepro_filepath)
evorepro_df = evorepro_df[evorepro_df['Organ name'] != '-']

po_mapper = {
    'Spore': ['plant spore'],
    'Female': ['plant egg cell', 'plant ovary', 'plant ovule'],
    'Leaf': ['leaf'],
    'Stem': ['stem'],
    'Male': ['pollen', 'microspore'],
    'Flower': ['flower'],
    'Seeds': ['seed'], # is a leaf node. No sublasses.
    'Root': ['root'], # does not have any meristem as children
    'Root meristem': ['root meristem'], # path does not coincide with root po term
    'Apical meristem': ['apical meristem']
}
evorepro_to_taxid = {
    'MAIZE': 4577,
    'ORYSA': 4530,
    'SOLLC': 4081,
    'ARATH': 3702,
    'MARPO': 3197,
    'PHYPA': 3218,
    'PICAB': 3329
}

# Exploratory
# def name_search(search):
#     relevant_terms = [term for term in po.terms if search in term.name[0].lower()]
#     names = [term.name[0] for term in relevant_terms]
#     return names

def merge_annotation_to_evorepro(evorepro_df, species, taxid, annotation_method=2):
    if annotation_method == 1:
        file = [file for file in annotation1_files if f"taxid{taxid}" in file][-1]
        annotation_path = f"{DATA_PATH}/preprocess/sra-annotations/{file}"
    else:
        file = [file for file in annotation2_files if f"taxid{taxid}" in file][-1]
        annotation_path = f"{DATA_PATH}/preprocess/sra-annotation2s/{file}"
    standard = evorepro_df[evorepro_df['Species'] == species][['ID', 'Organ name']]
    df = pd.read_csv(annotation_path, sep='\t', header=None, names=['runid', 'po_term'])
    df.dropna(inplace=True)
    merged = standard.merge(df, how='inner', left_on='ID', right_on='runid')
    merged.set_index('runid', inplace=True)
    merged.drop('ID', axis=1, inplace=True)
    merged.rename(columns={'Organ name': 'organ'}, inplace=True)
    return merged, standard.shape[0], df.shape[0]

def compare_annotation(series, po_mapper=po_mapper):
    ref_organ, po_term = series['organ'], series['po_term']
    # naive matching
    if ref_organ.lower() in po_term.lower():
        return True
    # tree traversal matching
    query_term = po.find_by_name(po_term)
    ref_terms = [po.find_by_name(name) for name in po_mapper[ref_organ]]
    for ref_term in ref_terms:
        if ref_term in [*query_term.superclasses()][0]:
            return True
    return False

def count_stats(species, taxid, method=2):
    merged, total_evorepro, total_annotated = merge_annotation_to_evorepro(evorepro_df, species, taxid, annotation_method=method)
    merged['matched'] = merged.apply(compare_annotation, axis=1)
    tp = merged['matched'].sum()
    total = merged.shape[0]
    stats = {
        'matched': tp,
        'total_compared': total,
        'p_matched': tp*100/total,
        'total_evorepro': total_evorepro, # Assumed truth
        'total_annotated': total_annotated # Not relevant
    }
    return stats

matches1 = {}
matches2 = {}
for species, taxid in evorepro_to_taxid.items():
    matches1[taxid] = count_stats(species, taxid, method=1)
    matches2[taxid] = count_stats(species, taxid, method=2)
method1 = pd.DataFrame(matches1).T
method2 = pd.DataFrame(matches2).T
method1 = method1.astype({'matched': 'int32',
                'total_compared': 'int32',
                'total_evorepro': 'int32',
                'total_annotated': 'int32'
                })
method2 = method2.astype({'matched': 'int32',
                'total_compared': 'int32',
                'total_evorepro': 'int32',
                'total_annotated': 'int32'
                })
# method1.to_csv("../annotation-benchmark/method1.txt",sep='\t', index_label='taxid')
# method2.to_csv("../annotation-benchmark/method2.txt",sep='\t', index_label='taxid')
