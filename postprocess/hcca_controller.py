import os
import datetime as dt

# species_job_list = ["Cpa", "Cre", "Osa", "Pab", "Sly", "Smo", "Vvi", "Zma", "Ath"]
species_job_list = ["Osa", "Pab", "Sly", "Smo", "Vvi", "Zma", "Ath"]

for species in species_job_list:
    start = dt.datetime.now()
    print('######################################################################')
    print(f"Running HCCA script for {species}...")
    print('######################################################################')
    # mapman_path = f"../plants/{species}/{species}.mercator.txt"
    hcca_path = f"../plants/{species}/{species}.hcca"
    hrr_path = f"../plants/{species}/{species}.hrr"
    os.system(f"python3 main/HCCA.py -i {hrr_path} -o {hcca_path} -s 3 -c 100")
    end = dt.datetime.now()
    time = end - start
    print(f"Time taken: {time}")
    print("Started:", start.strftime("%Y-%m-%d %H:%M:%S"))
    print("Finished:", end.strftime("%Y-%m-%d %H:%M:%S"))
    print()

print("Done for all! :)")
