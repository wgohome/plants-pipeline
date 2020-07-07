"""
Connects to ENA API to:
- Query taxanomic id based on species name
- Get RunIDs list for a given taxanomic id
- Stream Runtable (tsv .txt) for a given taxanomic id,
    returing data as pd df

Options: (stored as an attributes dict)
- returned fields
- output format

Hardcoded:
- Query filters
- Query to portal endpoint instead of browser endpoint
    (browser will download the data as xml)
"""

import pandas as pd
import requests
import warnings
import pdb

default_attributes = {'fields': 'run_accession, tissue_type,tissue_lib,sample_title,study_title,dev_stage,sample_description'}

# Exported functions
def get_taxanomic_id(scientific_name):
    """ Returns taxanomic id (str) queried from ENA based on scientific name """
    genus, species = scientific_name.lower().split()
    url = f"https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/{genus}%20{species}"
    response = requests.get(url)
    if response.status_code != 200:
        raise ValueError("scientific_name is invalid")
    if response.text == 'No results.':
        warnings.warn(f"No results. Tax id for {scientific_name} not found on ENA")
        return None
    tax_id = response.json()[0]['taxId']
    return tax_id

def get_ena_runtable(tax_id, query_attributes=default_attributes):
    response = query_ena(tax_id, query_attributes)
    ena_df = pd.DataFrame(response.json()) # Assumes json, will need another function for tsv
    return ena_df

def get_ena_runids(tax_id, query_attributes=default_attributes):
    ena_df = get_ena_runtable(tax_id, query_attributes)
    runids = ena_df['run_accession'].tolist()
    return runids

# Local functions
def query_ena(tax_id, query_attributes=default_attributes):
    output_format = query_attributes.get('output_format') or "json" # OR tsv. To get xml, use ena browser instead of portal.
    fields = query_attributes.get('fields') or "tissue_lib,tissue_type"
    url = f"https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=tax_tree({tax_id})%20AND%20instrument_platform=%22ILLUMINA%22&fields={fields}&format={output_format}"
    response = requests.get(url)
    if response.status_code != 200:
        raise ValueError("query invalid")
    return response

__all__ = ['get_taxanomic_id', 'get_ena_runids', 'get_ena_runtable']
