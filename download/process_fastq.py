# Module Header ################################################################
import sys
from os.path import realpath, dirname
import re

if __name__ == '__main__':
    abspath = realpath(dirname(__file__))
    parent_module = re.search('^(.*plants-pipeline)', abspath).group()
    sys.path.insert(0, parent_module)
################################################################################

# Module imports
import concurrent.futures
import datetime as dt
import os
import time
import pdb
# Relative imports of CONSTANTS in config/local_paths.py
from config.constants import *

def initiate_logfile():
    log_path = f"{DATA_PATH}/download/logs/time/{get_timestamp()}-time-log.log"
    with open(log_path, 'w') as f:
        f.write('\t'.join(['timestamp', 'runid', 'ascp_time', 'kallisto_time', 'library_layout']) + '\n')
    return log_path

def get_timestamp():
    return dt.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')

def get_fastq_routes(runid):
    """Returns tuple of trailing path of fastq file in vol1/fastq/ server's directory,
    for paired and unpaired libraries"""
    dir2 = ""
    if 9 < len(runid) <= 12:
        dir2 = "0" * (12 - len(runid)) + runid[-(len(runid) - 9):] + "/"
    dirs = runid[:6] + "/" + dir2 + runid
    return dirs + "/" + runid + "_1.fastq.gz", dirs + "/" + runid + ".fastq.gz"

def run_sh_cmd(cmd):
    """Returns sh_function which runs the shell command cmd string formatted with kawrgs."""
    def sh_function(**kwargs):
        """Returns runtime and os.system exit code from running shell command cmd,
        formatted with kawrgs."""
        start = time.time()
        kwargs['ASPERA_SSH_KEY'] = ASPERA_SSH_KEY
        exit_code = os.system(cmd.format(**kwargs))
        runtime = round(time.time() - start, 2)
        return runtime, exit_code
    return sh_function

ascp_transfer = run_sh_cmd("ascp -QTd -l 300m -P33001 -@ 0:1000000000 -i '{ASPERA_SSH_KEY}' era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/{route} '{out_path}'")
kallisto_quant = run_sh_cmd("kallisto quant -i '{idx_path}' -t 2 -o '{out_dir}' --single -l 200 -s 20 '{fastq_path}'")
kallisto_index = run_sh_cmd("kallisto index -i {idx_path} {cds_path}")

def dl_fastq(runid):
    """Attempts downloading runid fastq as paired file if possible, unpaired otherwise.
    Returns tuple of runtime of ascp download, string of library layout
    and out_path where fastq file is stored."""
    paired, unpaired = get_fastq_routes(runid)
    out_path = f"{DATA_PATH}/download/fastq/{'/'.join(unpaired.split('/')[:-1])}"
    runtime, _ = ascp_transfer(route=paired, out_path=out_path)
    if os.path.exists(f"{DATA_PATH}/download/fastq/{paired}"):
        return runtime, 'paired', paired
    runtime, _ = ascp_transfer(route=unpaired, out_path=out_path)
    if os.path.exists(f"{DATA_PATH}/download/fastq/{unpaired}"):
        return runtime, 'unpaired', unpaired
    return runtime, 'failed', ''

def run_a_job(runid, idx_path):
    """Downloads fastq file, quantify with kallisto and
    returns ascp_runtime, kal_runtime, run_layout"""
    ascp_runtime, layout, route = dl_fastq(runid)
    if layout == 'failed':
        return runid, ascp_runtime, 0, layout
    else:
        out_dir = f"{DATA_PATH}/download/kallisto_out/{'/'.join(route.split('/')[:-1])}"
        os.makedirs(out_dir, exist_ok=True)
        kal_runtime, _ = kallisto_quant(idx_path=idx_path,
            out_dir=out_dir, fastq_path=f"{DATA_PATH}/download/fastq/{route}")
        return runid, ascp_runtime, kal_runtime, layout

def write_log(to_write, log_path):
    with open(log_path, 'a') as f:
        f.write(to_write)

def parallelize(job_fn, runids, idx_path):
    """Executes job_fn on each elements in runids in parallel.
    Returns a list of the return values of job_fn."""
    with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
        futures = [executor.submit(job_fn, runid, idx_path) for runid in runids]
    results = []
    for f in concurrent.futures.as_completed(futures):
        results.append(f.result())
    return results

def process_batch(runids, idx_path):
    log_path = initiate_logfile()
    start = time.time()
    results = parallelize(run_a_job, runids, idx_path)
    for runid, ascp_runtime, kal_runtime, layout in results:
        to_write = f"{get_timestamp()}\t{runid}\t{ascp_runtime}\t{kal_runtime}\t{layout}\n"
        write_log(to_write, log_path)
    total_runtime = round(time.time() - start, 2)
    write_log(f"Total runtime\t{total_runtime}\n", log_path)

__all__ = ['process_batch', 'kallisto_index']
