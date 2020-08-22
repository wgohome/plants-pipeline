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
# abspath = realpath(dirname(__file__))
# parent_module = re.search('^(.*pipeline)', abspath).group()
dotenv_path = os.path.join(parent_module, '.env')
load_dotenv(dotenv_path)
NCBI_API_KEY = os.getenv('NCBI_API_KEY')

def esearch_genome(retstart=0, retmax=500):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=genome&retmode=json&term=viridiplantae[orgn:__txid33090]&retstart={retstart}&retmax={retmax}&usehistory=y&api_key={NCBI_API_KEY}"
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

def query_species_list():
    master_list = []
    i = 0
    retmax = 500
    while retmax == 500:
        webenv, count, retmax = esearch_genome(retstart=(i * 500))
        taxid_to_species = esummary_webenv(webenv, retstart=(i * 500))
        master_list.extend(taxid_to_species)
        i += 1
        sleep(0.1)
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

def make_species_report():
    master_list = query_species_list()
    print("Species list obtained")
    for i in range(len(master_list)):
    # master_list = master_list[:3]
    # for i in range(3):
        rna_count, illumina_rna_count = sra_exp_numbers(master_list[i]['taxid'], master_list[i]['species'])
        master_list[i]['rna_count'] = rna_count
        master_list[i]['illumina_rna_count'] = illumina_rna_count
        print(f"Done for #{i+1}/{len(master_list)}. {master_list[i]['species']}")
    master_df = pd.DataFrame(master_list)
    master_df.to_csv(f"{DATA_PATH}/preprocess/species-list/{iohelper.get_timestamp()}-taxid33090.txt", sep='\t', index=False)
    print("File saved in pipeline-data/preprocess/species-list")
    return master_df

master_df = make_species_report()
