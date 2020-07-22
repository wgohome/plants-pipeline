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
Connects to SRA API to:
- Download RunTable from RunInfoSelector with Selenium

Options: (stored as an attributes dict)
- returned fields (complications with SRA inconsistent headers!)
- output format

Hardcoded:
- Query filters
- Max of 100 000 runIDs, need to make the next sra query to continue getting
- Id of download button targetted by selenium, might change if sra modifies it
"""

import pandas as pd
import re
import requests
import os
import shutil
import time
from selenium.webdriver import Chrome
from selenium.webdriver.chrome.options import Options
import warnings
import pdb
# Relative imports
from config.constants import DATA_PATH
from preprocess import iohelper
from preprocess import ena

# CONSTANTS
BASE_PATH = os.path.dirname(os.path.realpath(__file__)) + '/'
SRA_WAITING_ROOM_PATH = f"{BASE_PATH}dlwaitingroom/" # Waiting room as default dir for chrome to download to, before moved to OUT_PATH
OUT_PATH = f"{DATA_PATH}/preprocess/sra-runtables/"

# Local variables
default_attributes = {}

# Exported functions
"""
Usage notes:
1. get_sra_runtable => Make sra query based on filters, download runtable csv
2. read_sra_runtable => Read runtable csv into pd df
3. process_sra_runtable => Both step 1 and 2 together, returning the pd df
* tax_id can be obtained with `ena.get_taxanomic_id('Genus species')`
"""
get_taxanomic_id = ena.get_taxanomic_id

def process_sra_runtable(tax_id, species, query_attributes=default_attributes):
    runtable_path = get_sra_runtable(tax_id, species, query_attributes)
    if runtable_path == None:
        return None
    sra_df = read_sra_runtable(runtable_path, query_attributes)
    return sra_df

def get_sra_runtable(tax_id, species, query_attributes=default_attributes):
    sra_runtable_path = f"{DATA_PATH}/preprocess/sra-runtables/{iohelper.species_shortform(species)}_sra_runtable.txt"
    if os.path.exists(sra_runtable_path):
        return sra_runtable_path
    webenv = make_sra_query(tax_id, species, query_attributes)
    if webenv == None:
        return None
    runtable_path = dl_sra_runtable(tax_id, species, webenv)
    return runtable_path

def read_sra_runtable(runtable_path, query_attributes=default_attributes):
    sra_df = pd.read_csv(runtable_path, sep=',', header=0, index_col=False, dtype='string')
    return sra_df

# Local functions
def make_sra_query(tax_id, species, query_attributes=default_attributes):
    species = species.replace(' ', '%20')
    output_format = query_attributes.get('output_format') or 'json'
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&retmode={output_format}&term=\"{species}\"[orgn:__txid{tax_id}]+AND+\"biomol\%20rna\"[properties]+AND+\"platform\%20Illumina\"[properties]&retmax=100000&usehistory=y"
    response = requests.get(url)
    if response.status_code != 200:
        raise ValueError("query invalid")
    if response.text == 'No results.':
        warnings.warn(f"No results. Tax id {tax_id} ({species}) not found on ENA")
        return None
    sra_json = pd.DataFrame(response.json())
    webenv = sra_json["esearchresult"]['webenv']
    return webenv

def dl_sra_runtable(tax_id, species, webenv):
    # dl_path = os.path.dirname(os.path.realpath(__file__))
    # webenv = make_sra_query(tax_id, species)
    urlRunSelector = f"https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=1&WebEnv={webenv}"
    runtable_path = download_from_webpage(species, urlRunSelector)
    return runtable_path

def download_from_webpage(species, urlRunSelector):
    opts = Options()
    opts.add_argument('headless')
    # opts.add_argument('--window-size=1920x1080')
    # assert opts.headless
    prefs = {'download.default_directory': SRA_WAITING_ROOM_PATH}
    opts.add_experimental_option('prefs', prefs)
    browser = Chrome(options=opts, executable_path='chromedriver')
    browser.get(urlRunSelector)
    button = check_page_load(browser)
    button.click()
    runtable_path = check_download_status(species)
    browser.close()
    return runtable_path

def check_page_load(browser):
    counter = 0
    while counter < 24: #Limit page load time to 2 mins
        time.sleep(5)
        counter += 1
        try:
            button = browser.find_element_by_id("t-rit-all")
            return button
        except:
            continue
    raise TimeoutError("Page load tiemout")

def check_download_status(species):
    # Rename file after download is complete, move file to out/sra_runtables
    dirs = []
    attempt = 0
    while not len(dirs):
        time.sleep(5)
        dirs = [dir for dir in os.listdir(SRA_WAITING_ROOM_PATH) if dir[-4:] == '.txt'] # Check that download is complete
        attempt += 1
        if attempt > 24: # Limit dl time to 2 minutes
            raise TimeoutError("Download attempt timeout")
    runtable_path = move_downloaded_file(species, dirs)
    return runtable_path

def move_downloaded_file(species, dirs):
    origin_name = f"{SRA_WAITING_ROOM_PATH}{dirs[0]}"
    Spe = iohelper.species_shortform(species)
    new_name = f"{SRA_WAITING_ROOM_PATH}{Spe}_sra_runtable.txt"
    target_path = f"{OUT_PATH}{Spe}_sra_runtable.txt"
    os.rename(origin_name, new_name)
    shutil.move(new_name, target_path) # change target to full path instead of dir path to allow overwrite
    return target_path

__all__ = ['get_sra_runtable', 'read_sra_runtable', 'process_sra_runtable', 'get_taxanomic_id']
