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
# Relative imports of CONSTANTS in config/constants.py
from config.constants import *

def initiate_logfile(log_type, headers, spe=''):
    log_path = f"{DATA_PATH}/download/logs/{log_type}/{get_timestamp()}-{spe}{log_type}.log"
    with open(log_path, 'w') as f:
        f.write('\t'.join(headers) + '\n')
    return log_path

def write_log(to_write, log_path):
    with open(log_path, 'a') as f:
        f.write(to_write)

def get_timestamp():
    return dt.datetime.now().strftime('%Y%m%d-%H%M%S')

def get_fastq_routes(runid):
    """Returns tuple of trailing path of fastq file in vol1/fastq/ server's directory, for paired and unpaired libraries,
    and also file names for paired and unpaired libraries"""
    p_file, up_file = f"{runid}_1.fastq.gz", f"{runid}.fastq.gz"
    dir2 = ""
    if 9 < len(runid) <= 12:
        dir2 = "0" * (12 - len(runid)) + runid[-(len(runid) - 9):] + "/"
    dirs = f"{runid[:6]}/{dir2}{runid}/"
    return f"{dirs}{p_file}", f"{dirs}{up_file}", p_file, up_file

def wrap_sh_command(cmd):
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

ascp_transfer = wrap_sh_command("ascp -QTd -l 300m -P33001 -@ 0:1000000000 -i '{ASPERA_SSH_KEY}' era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/{route} '{out_path}'")
kallisto_quant = wrap_sh_command("kallisto quant -i '{idx_path}' -t 2 -o '{out_dir}' --single -l 200 -s 20 '{fastq_path}'")
kallisto_index = wrap_sh_command("kallisto index -i {idx_path} {cds_path}")

def dl_fastq(runid):
    """Attempts downloading runid fastq as paired file if possible, unpaired otherwise.
    Returns tuple of runtime of ascp download, string of library layout,
    and filename of the downloaded fastq."""
    p_route, up_route, p_file, up_file = get_fastq_routes(runid)
    f"{runid}_1.fastq.gz", f"{runid}.fastq.gz"
    out_path = f"{DATA_PATH}/download/fastq_tmp/"
    runtime, _ = ascp_transfer(route=p_route, out_path=out_path)
    if os.path.exists(f"{DATA_PATH}/download/fastq-tmp/{p_file}"):
        return runtime, 'paired', p_file
    runtime, _ = ascp_transfer(route=up_route, out_path=out_path)
    if os.path.exists(f"{DATA_PATH}/download/fastq-tmp/{up_file}"):
        return runtime, 'unpaired', up_file
    return runtime, 'failed', ''

def run_a_job(runid, idx_path, init_log_path, runtime_log_path):
    """Downloads fastq file, quantify with kallisto and
    returns runid, ascp_runtime, kal_runtime, run_layout.
    Writes initiation time for this runid in init_log_path and
    write runtimes and layout in runtime_log_path"""
    write_log(f"{get_timestamp()}\t{runid}\n", init_log_path)
    ascp_runtime, layout, filename = dl_fastq(runid)
    if layout == 'failed':
        kal_runtime = 0
    else:
        kal_runtime, _ = kallisto_quant(idx_path=idx_path,
            out_dir=f"{DATA_PATH}/download/kallisto-tmp/{runid}/",
            fastq_path=f"{DATA_PATH}/download/fastq-tmp/{filename}")
    fields = [get_timestamp(), runid, str(ascp_runtime), str(kal_runtime), layout]
    write_log('\t'.join(fields) + '\n', runtime_log_path)
    return runid, ascp_runtime, kal_runtime, layout

def parallelize(job_fn, runids, idx_path, init_log_path):
    """Executes job_fn on each elements in runids in parallel.
    Returns a list of the return values of job_fn."""
    with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
        futures = [executor.submit(job_fn, runid, idx_path) for runid in runids]
    results = []
    for f in concurrent.futures.as_completed(futures):
        results.append(f.result())
    return results

def process_batch(runids, idx_path, spe):
    init_log_path = initiate_logfile('initiation', ['timestamp', 'runid'], spe=f"{spe}-")
    runtime_log_path = initiate_logfile('runtime', ['timestamp', 'runid', 'ascp_time', 'kallisto_time', 'library_layout'], spe=f"{spe}-")
    batch_start = time.time()
    results = parallelize(run_a_job, runids, idx_path, init_log_path)
    batch_runtime = round(time.time() - batch_start, 2)
    write_log(f"Total runtime\t{batch_runtime}\n", runtime_log_path)

__all__ = ['process_batch', 'kallisto_index']

