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
import re
import pandas as pd
import numpy as np
import inflect
from functools import reduce
import pdb
# Relative imports
from preprocess.po_parser.po import po # po is an Ontology class instance
from preprocess import iohelper

# Constants
INFL = inflect.engine() # For pluralization

def annotate(df, species, db='ena'):
    # species refer to full species name
    global HEADERS
    HEADERS = df.columns.tolist()
    out_path = iohelper.create_annotation_file(species, db)
    cols1, cols2, cols3 = get_cols(db) # cols4 will be the remaining columns
    for index, series in df.iterrows():
        match = collect_matches(series, cols1, cols2, cols3)
        index = series.run_accession if db == 'ena' else series.Run
        # to_write = index + '\t' + '\t'.join(po_names) + '\n'
        po_name = match.name[0] if match != "" else ""
        to_write = f"{index}\t{po_name}\n"
        iohelper.write_annotation(to_write, out_path)

def get_cols(db='ena'):
    # Values here to be manually editted if desired
    if db == 'ena':
        cols1 = ['tissue_type', 'tissue_lib']
        cols2 = ['sample_title', 'study_title', 'dev_stage']
        cols3 = ['sample_description']
    elif db == 'sra':
        cols1 = ['tissue', 'sample_name', 'tissue_type', 'sample_type', 'sample_title', 'plant_structure', 'Cell_type', 'organism_part', 'Sample_Tissue', 'tissue-type']
        cols2 = ['title', 'Sample Name', 'Library Name', 'dev_stage', 'Age', 'source_name', 'Ecotype', 'Developmental_stage', 'plant_developmental_stage', 'growth_stage', 'stage', 'tissue_position', 'Development_stage', 'develomental_stage', 'sample_comment', 'developmental_stge']
        cols3 = ['Description']
    return cols1, cols2, cols3

def collect_matches(series, cols1, cols2, cols3):
    # Attempt matching cols1 first, then cols2
    for col in [cols1, cols2]:
        matches = iterate_cols_matching(series, col)
        if matches:
            return find_best_match(matches)
    # Finally match to long entries in cols3
    matches = iterate_cols_matching_long(series, cols3)
    return find_best_match(matches) if matches else ""

def iterate_cols_matching(series, cols):
    matches = []
    for col_raw in cols:
        for col in validate_col(col_raw, series):
            col_match = find_matches(series[col])
            if col_match: matches.append(col_match)
    return matches

def iterate_cols_matching_long(series, cols3):
    matches = []
    for col_raw in cols3:
        for col in validate_col(col_raw, series):
            query = series.loc[col]
            if pd.isnull(query):
                continue
            if len(query) > 100:
                rg = r"([\w\s(\. )]*sample[\w\s(\. )]*)"
                queries = re.findall(rg, query)
            else:
                queries = [query]
            for q in queries:
                match = find_matches(q)
                if match: matches.append(match)
    return matches

def validate_col(col, series):
    # Returns a list of col headers that are valid for this series/runtable
    perms = [col.lower(), col.upper(), col.capitalize(), col.title()]
    perms = list(filter(lambda x: (x in HEADERS), perms))
    return perms

def find_matches(query):
    # Match a single query (in a cell) against all PO terms
    # Choose the most specific term that matches this entry
    matches = []
    if pd.isnull(query):
        return matches
    query = query.lower()
    for term in po.terms:
        po_name = re.sub(r"(?i)\s*plant\s*", " ", term.name[0]).strip().lower()
        regex = re.compile(rf"\b{po_name}\b", re.I)
        regex_plural = re.compile(rf"\b{INFL.plural(po_name)}\b", re.I)
        match = re.findall(regex, query)
        match_plural = re.findall(regex_plural, query)
        if match or match_plural:
            matches.append(term)
    if matches:
        return choose_most_specific_term(matches)
    else:
        return None

def choose_most_specific_term(matches):
    """Return most specific term out of matches"""
    path_length = lambda path: len(path)
    # Take longest path of match term to represent the specificity of that term
    longest_path = lambda match: sorted([*match.superclasses()], key=path_length)[-1]
    # Take match term with the longest representaive path
    most_specific_term = lambda match: len(longest_path(match))
    return sorted(matches, key=most_specific_term)[-1]

def choose_least_specific_term(matches):
    path_length = lambda path: len(path)
    # Take longest path of match term to represent the specificity of that term
    longest_path = lambda match: sorted([*match.superclasses()], key=path_length)[0]
    # Take match term with the longest representaive path
    most_specific_term = lambda match: len(longest_path(match))
    return sorted(matches, key=most_specific_term)[0]

find_best_match = choose_least_specific_term

# def find_best_match(matches):
#     # Find the nodes closest to the leaf nodes.
#     # If there are matches from different lines, choose the common last node.
#     if matches == []: # reduce() cannot take an empty iterable as argument
#         return ""
#     matches_treepaths = [find_parents(match) for match in matches]
#     best_common_treepath = reduce(compare_treepath, matches_treepaths)
#     return best_common_treepath[-1] if best_common_treepath else ""

# def find_parents(match):
#     treepath = [*match.superclasses()]
#     treepath.reverse()
#     return treepath

# def compare_treepath(path1, path2):
#     merged_treepath = []
#     for node1, node2 in zip(path1, path2):
#         if node1 == node2:
#             merged_treepath.append(node1)
#         else:
#             break
#     else:
#         # Only exceuted if break is not triggered, i.e. path1, path2 ⊆ same path
#         # Returns the more specific between path1 and path2
#         return path1 if (len(path1) > len(path2)) else path2
#     # path1 and path2 are separate paths, truncate up to the last common node
#     return merged_treepath
