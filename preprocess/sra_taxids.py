# Module Header ################################################################
import sys
from os.path import realpath, dirname
import re

if __name__ == '__main__':
    abspath = realpath(dirname(__file__))
    parent_module = re.search('^(.*pipeline)', abspath).group()
    sys.path.insert(0, parent_module)

################################################################################
"""
Connects to NCBI Taxanomy eUtils to:
-
"""
import requests
import pandas as pd
import datetime as dt
import xml.etree.ElementTree as et
import pdb

import ena

# from os.path import realpath, dirname
# import re
# abspath = realpath(dirname(__file__))
# parent_module = re.search('^(.*pipeline)', abspath).group()

study = 'study1'
timestamp = str(dt.datetime.now()).replace(' ', '_')
write_path = f"{parent_module}/../data/{study}/in/{timestamp}sra_taxids.txt"

def query_genome():
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=genome&retmode=json&term=viridiplantae[Organism]+AND+complete[Status]+AND+genome[Properties]&retmax=100000&usehistory=y"
    response = requests.get(url)
    if response.status_code != 200:
        raise ValueError("search terms are invalid, GET request failed")
    webenv = response.json()["esearchresult"]["webenv"]
    return webenv

def get_taxids(webenv):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=genome&retmode=json&query_key=1&WebEnv={webenv}"
    response = requests.get(url)
    if response.status_code != 200:
        raise ValueError("search terms are invalid, GET request failed")
    results = response.json()["result"]
    uids_list = results.pop('uids')
    taxids = [summ['taxid'] for summ in results.values()]
    names = [summ['organism_name'] for summ in results.values()]
    return taxids, names

def get_sra_nums(taxid):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&retmode=json&term={taxid}[orgn:__txid{taxid}] +AND+biomol%20rna[properties]+AND+platform%20Illumina[properties]&retmax=100000&usehistory=y"
    response = requests.get(url)
    if response.status_code != 200:
        raise ValueError("search terms are invalid, GET request failed")
    num = response.json()["esearchresult"]["count"]
    return num

# def query_genome_all():
#     retstart = 0
#     webenvs = []
#     while True:
#         print('retstart', retstart)
#         url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=genome&retmode=json&term=viridiplantae[Organism]+AND+genome[Properties]&retstart={retstart}&retmax=500&usehistory=y"
#         response = requests.get(url)
#         if response.status_code != 200:
#             raise ValueError("search terms are invalid, GET request failed")
#         results = response.json()
#         webenv = results["esearchresult"]["webenv"]
#         webenvs.append(webenv)
#         retstart += 500
#         results_count = int(results['esearchresult']['retmax'])
#         # total_count = int(results['esearchresult']['count'])
#         print('results', results_count)
#         # print('total', total_count)
#         print()
#         if results_count < 500:
#             break
#     print('done')
#     return webenvs

def query_genome_all():
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=genome&retmode=json&term=viridiplantae[Organism]+AND+genome[Properties]&retmax=100000&usehistory=y"
    response = requests.get(url)
    if response.status_code != 200:
        raise ValueError("search terms are invalid, GET request failed")
    webenv = response.json()["esearchresult"]["webenv"]
    return webenv

def get_taxids_all(webenvs):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=genome&query_key=1&WebEnv={webenv}"
    response = requests.get(url)
    if response.status_code != 200:
        raise ValueError("search terms are invalid, GET request failed")
    # results = response.json()["result"]
    # uids_list = results.pop('uids')
    # taxids = [summ['taxid'] for summ in results.values()]
    # names = [summ['organism_name'] for summ in results.values()]
    # taxids_ls.append(taxids)
    # names_ls.append(names)
    root = et.fromstring(response.content)
    names = [child[1].text for child in root]
    return names

# webenv = query_genome()
# taxids, names = get_taxids(webenv)
# sra_nums = [get_sra_nums(taxid) for taxid in taxids]
# df = pd.DataFrame({'taxids': taxids, 'species_names': names, 'sra_exp_numbers': sra_nums})
# df.to_csv("test.txt", sep='\t', index=False)

webenv = query_genome_all()
names = get_taxids_all(webenv)
# to_write = [name + '\n' for name in names]
# with open('all_plants.txt','w') as f:
#     f.writelines(to_write)
hybrids = [name for name in names if len(name.split())>2]
single_names = [name for name in names if len(name.split()) < 2]
trimmed_names = list(set(names) - set(hybrids) - set(single_names))
taxids = [ena.get_taxanomic_id(name) for name in trimmed_names]
sra_nums = [get_sra_nums(taxid) for taxid in taxids]
df = pd.DataFrame({'taxids': taxids, 'species_names': trimmed_names, 'sra_exp_numbers': sra_nums})
df.to_csv("all_plants_nofilter.txt", sep='\t', index=False)
