# Module Header ################################################################
import sys
from os.path import realpath, dirname
import re

if __name__ == '__main__':
    abspath = realpath(dirname(__file__))
    parent_module = re.search('^(.*pipeline)', abspath).group()
    sys.path.insert(0, parent_module)

################################################################################

# Packages
import re
import pandas as pd
import numpy as np
from functools import reduce
import pdb
# Relative imports
from datasource.po import po_terms, infl, po
from dataprocess import iohelper

def annotate(df, species, study, out_path, db='ena'):
    global HEADERS
    HEADERS = df.columns.tolist()
    # Spe = iohelper.species_shortform(species)
    with open(out_path, 'w') as f:
        f.write("")
    cols1, cols2, cols3 = get_cols(db) # cols4 will be the remaining columns
    for index, series in df.iterrows():
        # print(index)
        row_name = series.run_accession if db == 'ena' else series.Run
        to_write_list = [row_name]
        super_matches = collect_matches(series, cols1, cols2, cols3)
        for i, match_tuple in enumerate(super_matches):
            # print("cols" + str(i + 1))
            po_name = match_tuple[0].name if match_tuple[0] else ""
            origin_text = match_tuple[1]
            # print(po_name)
            # print(origin_text)
            to_write_list.extend([po_name, origin_text])
        to_write = '\t'.join(to_write_list) + '\n'
        # print()
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
    super_matches = []
    # Attempt matching cols1 first, then cols2
    for col in [cols1, cols2]:
        matches = iterate_cols_matching(series, col)
        to_append = find_best_match(matches) if matches else ('','')
        super_matches.append(to_append)
    # Finally match to long entries in cols3
    matches = iterate_cols_matching_long(series, cols3)
    to_append = find_best_match(matches) if matches else ('','')
    super_matches.append(to_append)
    return super_matches

def iterate_cols_matching(series, cols):
    matches = []
    for col_raw in cols:
        for col in validate_col(col_raw, series):
            col_match = find_matches(series[col])
            if col_match: matches.append((col_match, series[col]))
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
                if match: matches.append((match, series[col]))
    return matches

def validate_col(col, series):
    # Returns a list of col headers that are valid for this series/runtable
    perms = [col.lower(), col.upper(), col.capitalize(), col.title()]
    perms = list(filter(lambda x: (x in HEADERS), perms))
    return perms

def find_matches(query):
    # Match a single query (in a cell) against all PO terms
    matches = []
    if pd.isnull(query):
        return matches
    query = query.lower()
    for term in po_terms:
        po_name = re.sub(r"(?i)\s*plant\s*", " ", term.name).strip().lower()
        regex = re.compile(rf"\b{po_name}\b", re.I)
        regex_plural = re.compile(rf"\b{infl.plural(po_name)}\b", re.I)
        match = re.findall(regex, query)
        match_plural = re.findall(regex_plural, query)
        if match or match_plural:
            matches.append(term)
    if matches:
        return choose_more_specific_term(matches)
    else:
        return None

def choose_more_specific_term(matches):
    path_length = lambda match: len([*match.superclasses()])
    matches.sort(key=path_length, reverse=True)
    return matches[0]

def find_best_match(matches):
    # Find the nodes closest to the leaf nodes.
    # If there are matches from different lines, choose the common last node.
    if matches == []: # reduce() cannot take an empty iterable as argument
        return ""
    elif len(matches) == 1:
        return matches[0]
    matches_treepaths = [(find_parents(match[0]), match[1]) for match in matches]
    best_common_treepath = reduce(compare_treepath, matches_treepaths)
    if best_common_treepath[0]:
        best_match_tuple = (best_common_treepath[0][-1], best_common_treepath[1])
    else:
        best_match_tuple = ('', '')
    return best_match_tuple

def find_parents(match):
    treepath = [*match.superclasses()]
    treepath.reverse()
    return treepath

def compare_treepath(tuple1, tuple2):
    path1 = tuple1[0]
    path2 = tuple2[0]
    merged_treepath = []
    for node1, node2 in zip(path1, path2):
        if node1 == node2:
            merged_treepath.append(node1)
        else:
            break
    else:
        # Only exceuted if break is not triggered, i.e. path1, path2 ⊆ same path
        # Returns the more specific between path1 and path2
        return (path1, tuple1[1]) if (len(path1) > len(path2)) else (path2, tuple2[1])
    # path1 and path2 are separate paths, truncate up to the last common node
    return (merged_treepath, tuple1[1]) if merged_treepath else ('', '')
