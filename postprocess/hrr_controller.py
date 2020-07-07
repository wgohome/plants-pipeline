"""
controller.py

HCCA worflow main script
-> Takes in species input file
-> Create directories for respective output
->
"""

import os
import datetime as dt

# Add argv

# species_job_list = ["Cpa", "Cre", "Osa", "Pab", "Sly", "Smo", "Vvi", "Zma", "Ath"]
species_job_list = ["Egu", "Nta"]

for species in species_job_list:
    start = dt.datetime.now()
    print('######################################################################')
    print(f"Running HRR script for {species}...")
    print('######################################################################')
    # mapman_path = f"../plants/{species}/{species}.mercator.txt"
    tpm_path = f"../plants/{species}/{species}_matrix.txt"
    hrr_path = f"../plants/{species}/{species}.hrr"
    os.system(f"python3 main/v3.3HRR.py -i {tpm_path} -o {hrr_path}")
    end = dt.datetime.now()
    time = end - start
    print(f"Time taken: {time}")
    print("Started:", start.strftime("%Y-%m-%d %H:%M:%S"))
    print("Finished:", end.strftime("%Y-%m-%d %H:%M:%S"))
    print()


print("Done for all! :)")
