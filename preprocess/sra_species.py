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
Functions that connect to NCBI Taxanomy eUtils to:
- Get species list for all viridiplantae
- Query the number of sra entires for each species
"""
import os
import datetime as dt
from dotenv import load_dotenv
import pandas as pd
import requests
from time import sleep
import pdb
# Relative imports
from preprocess import iohelper
from config.constants import DATA_PATH

# Get NCBI API KEY
abspath = realpath(dirname(__file__))
parent_module = re.search('^(.*pipeline)', abspath).group()
dotenv_path = os.path.join(parent_module, '.env')
load_dotenv(dotenv_path)
NCBI_API_KEY = os.getenv('NCBI_API_KEY')

def esearch_genome(name, taxid, retstart=0, retmax=500):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=genome&retmode=json&term={name}[orgn:__txid{taxid}]&retstart={retstart}&retmax={retmax}&usehistory=y&api_key={NCBI_API_KEY}"
    response = requests.get(url)
    if response.status_code != 200:
        raise ValueError("search terms are invalid, GET request failed")
    webenv = response.json()["esearchresult"]["webenv"]
    count = int(response.json()["esearchresult"]["count"]) # total available records in db
    retmax = int(response.json()["esearchresult"]["retmax"]) # total record returned in this query
    return webenv, count, retmax

def esummary_webenv(webenv, retstart=0, retmax=500):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=genome&retmode=json&retstart={retstart}&retmax={retmax}&query_key=1&WebEnv={webenv}&api_key={NCBI_API_KEY}"
    response = requests.get(url)
    if response.status_code != 200:
        raise ValueError("search terms are invalid, GET request failed")
    results = response.json()["result"]
    return [{'taxid': attrs['taxid'], 'species': attrs['organism_name']} for key, attrs in results.items() if key != 'uids']

def query_species_list(name, taxid):
    master_list = []
    i = 0
    retmax = 500
    while retmax == 500:
        webenv, count, retmax = esearch_genome(name, taxid, retstart=(i * 500))
        taxid_to_species = esummary_webenv(webenv, retstart=(i * 500))
        master_list.extend(taxid_to_species)
        i += 1
        sleep(0.2)
    return master_list

def sra_exp_numbers(taxid, species):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&retmode=json&term={species}%5Borgn:__txid{taxid}%5D+AND+biomol rna[Properties]&usehistory=y&api_key={NCBI_API_KEY}"
    response = requests.get(url)
    if response.status_code != 200:
        raise ValueError("search terms are invalid, GET request failed")
    rna_count = int(response.json()['esearchresult']['count'])
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&retmode=json&term={species}%5Borgn:__txid{taxid}%5D+AND+biomol rna[Properties]+AND+platform illumina[Properties]&usehistory=y&api_key={NCBI_API_KEY}"
    response = requests.get(url)
    if response.status_code != 200:
        raise ValueError("search terms are invalid, GET request failed")
    illumina_rna_count = int(response.json()['esearchresult']['count'])
    sleep(0.2)
    return rna_count, illumina_rna_count

def make_species_report(name, taxid):
    master_list = query_species_list(name, taxid)
    print("Species list obtained")
    # 1. Check if last run of species report is terminated halfway
    species_list_files = sorted([file for file in os.listdir(f"{DATA_PATH}/preprocess/species-list/") if f"taxid{taxid}" in file])
    if species_list_files:
        species_list_path = f"{DATA_PATH}/preprocess/species-list/{species_list_files[-1]}"
        latest_df = pd.read_csv(species_list_path, sep='\t')
        start = latest_df.shape[0]
    else:
        start = 0
        species_list_path = f"{DATA_PATH}/preprocess/species-list/{iohelper.get_timestamp()}-taxid{taxid}.txt"
    # 2. Query numbers for each species
    for i in range(start, len(master_list)):
        rna_count, illumina_rna_count = sra_exp_numbers(master_list[i]['taxid'], master_list[i]['species'])
        master_list[i]['rna_count'] = rna_count
        master_list[i]['illumina_rna_count'] = illumina_rna_count
        df_tmp = pd.DataFrame([master_df[i]])
        header = True if i == 0 else False
        df_tmp.to_csv(species_list_path, sep='\t', index=False, mode='a', header=header)
        print(f"Queried for #{i+1}/{len(master_list)}. {master_list[i]['species']}")
    # 3. If all species in master_list queried, sort through by number of RNA experiments available
    master_df = pd.from_csv(species_list_path, sep='\t')
    if len(master_list) == master_df.shape[0]:
        master_df.sort_values(by=['illumina_rna_count', 'rna_count'], ascending=False, inplace=True)
        master_df['cds_link'] = None
        master_df.to_csv(f"{DATA_PATH}/preprocess/job-list/{iohelper.get_timestamp()}-taxid{taxid}.txt", sep='\t', index=False)
        print("TO EDIT: file saved in pipeline-data/preprocess/job-list, add cds links for species to be processed.")
    return master_df

__all__ = ['make_species_report']
