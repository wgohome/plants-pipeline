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
from config.constants import DATA_PATH, ASPERA_SSH_KEY
import helpers

def get_ftp_paths(runid):
    p_route, up_route, _, _ = helpers.get_fastq_routes(runid)
    return f"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{p_route}", f"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{up_route}"

def run_bash_command(runid, cmd):
    """Creates a .sh file for cmd, runs it with bash and deletes the .sh file"""
    bash_script_path = f"{DATA_PATH}/download/bash-tmp/{runid}.sh"
    with open(bash_script_path, 'w') as f:
        f.write(cmd)
    os.system(f"bash {bash_script_path}")
    os.remove(bash_script_path)

def wrap_sh_command(cmd, bash=False):
    """Returns sh_function which runs the shell command cmd string formatted with kawrgs."""
    def sh_function(**kwargs):
        """Returns runtime and os.system exit code from running shell command cmd,
        formatted with kawrgs."""
        start = time.time()
        kwargs['ASPERA_SSH_KEY'] = ASPERA_SSH_KEY
        kwargs['DATA_PATH'] = DATA_PATH
        if bash:
            exit_code = run_bash_command(kwargs['runid'], cmd.format(**kwargs))
        else:
            exit_code = os.system(cmd.format(**kwargs))
        runtime = round(time.time() - start, 2)
        return runtime, exit_code, cmd.format(**kwargs)
    return sh_function

ascp_transfer = wrap_sh_command("ascp -QT -l 300m -P33001 -@ 0:1000000000 -i '{ASPERA_SSH_KEY}' era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/{route} '{out_path}'")
kallisto_quant = wrap_sh_command("kallisto quant -i '{idx_path}' -t 2 -o '{out_dir}' --single -l 200 -s 20 '{fastq_path}'")
kallisto_index = wrap_sh_command("kallisto index -i {idx_path} {cds_path}")
kallisto_quant_curl = wrap_sh_command("kallisto quant -i {idx_path} -o {out_dir} --single -l 200 -s 20 -t 2 <(curl -L -r 0-1000000000 -m 600 --speed-limit 1000000 --speed-time 30 '{ftp_path}' 2> '{DATA_PATH}/download/kallisto-tmp/{runid}.log')", bash=True)

def dl_fastq(runid, layout):
    """Attempts downloading runid fastq as paired file if possible, unpaired otherwise.
    Returns tuple of runtime of ascp download, string of library layout,
    and filename of the downloaded fastq."""
    p_route, up_route, p_file, up_file = helpers.get_fastq_routes(runid)
    out_path = f"{DATA_PATH}/download/fastq-tmp/"
    if layout.upper() == 'PAIRED':
        runtime, _, cmd = ascp_transfer(route=p_route, out_path=out_path)
        if os.path.exists(f"{DATA_PATH}/download/fastq-tmp/{p_file}"):
            return runtime, 'paired', p_file, cmd
    else:
        runtime, _, cmd = ascp_transfer(route=up_route, out_path=out_path)
        if os.path.exists(f"{DATA_PATH}/download/fastq-tmp/{up_file}"):
            return runtime, 'unpaired', up_file, cmd
    return runtime, 'failed', '', cmd

def ascp_job(runid, layout, idx_path, init_log_path, runtime_log_path):
    """One of two download modes to choose from.
    Downloads fastq file, quantify with kallisto and
    returns runid, ascp_runtime, kal_runtime, run_layout.
    Writes initiation time for this runid in init_log_path and
    write runtimes and layout in runtime_log_path"""
    helpers.write_log(f"{helpers.get_timestamp()}\t{runid}\n", init_log_path)
    ascp_runtime, layout, filename, ascp_cmd = dl_fastq(runid, layout)
    if layout == 'failed':
        kal_runtime = 0
    else:
        out_dir = f"{DATA_PATH}/download/kallisto-tmp/{runid}/"
        os.makedirs(out_dir, exist_ok=True)
        kal_runtime, _, kal_cmd= kallisto_quant(idx_path=idx_path,
            out_dir=out_dir,
            fastq_path=f"{DATA_PATH}/download/fastq-tmp/{filename}")
    fields = [helpers.get_timestamp(), runid, str(ascp_runtime), str(kal_runtime), layout, ascp_cmd, kal_cmd]
    helpers.write_log('\t'.join(fields) + '\n', runtime_log_path)
    return [runid, ascp_runtime, kal_runtime, layout, ascp_cmd, kal_cmd]

def check_error_type(runid):
    """curl_job helper"""
    log_path = f"{DATA_PATH}/download/kallisto-tmp/{runid}.log"
    log = open(log_path,'r').read()
    if 'curl: (28)' in log:
        return 'operation_timeout'
    elif 'curl: (78)' in log:
        return 'remote_file_not_found'
    elif 'curl: (9)' in log:
        return 'remote_access_denied'
    else:
        return 'no_download_error'

def check_zero_processed(runid):
    """curl_job helper"""
    run_info_path = f"{DATA_PATH}/download/kallisto-tmp/{runid}/run_info.json"
    return '"n_processed": 0' in open(run_info_path,'r').read()

def curl_job(runid, layout, idx_path, init_log_path, runtime_log_path):
    """One of two download modes to choose from.
    Streams fastq with curl and directly pass to kallisto quant.
    Bypasses downloading of fastq."""
    helpers.write_log(f"{helpers.get_timestamp()}\t{runid}\n", init_log_path)
    p_ftp_path, up_ftp_path = get_ftp_paths(runid)
    out_dir = f"{DATA_PATH}/download/kallisto-tmp/{runid}/"
    os.makedirs(out_dir, exist_ok=True)
    if layout.upper() == 'PAIRED':
        runtime, _, _ = kallisto_quant_curl(runid=runid, idx_path=idx_path, out_dir=out_dir, ftp_path=p_ftp_path)
        # if not check_zero_processed(runid):
        error_type = check_error_type(runid)
        if error_type == 'no_download_error':
            fields = [helpers.get_timestamp(), runid, str(runtime), 'paired', p_ftp_path]
        else:
            fields = [helpers.get_timestamp(), runid, str(runtime), error_type, up_ftp_path]
        os.remove(f"{DATA_PATH}/download/kallisto-tmp/{runid}.log")
    else:
        runtime, _, _ = kallisto_quant_curl(runid=runid, idx_path=idx_path, out_dir=out_dir, ftp_path=up_ftp_path)
        error_type = check_error_type(runid)
        if error_type == 'no_download_error':
            fields = [helpers.get_timestamp(), runid, str(runtime), 'unpaired']
        else:
            fields = [helpers.get_timestamp(), runid, str(runtime), error_type]
        os.remove(f"{DATA_PATH}/download/kallisto-tmp/{runid}.log")
    helpers.write_log('\t'.join(fields) + '\n', runtime_log_path)
    return fields

def parallel_loop(job_fn, runids, layouts, idx_path, init_log_path, runtime_log_path, workers=8):
    """One of two loop functions to choose from.
    Executes job_fn on each elements in runids in parallel.
    Returns a list of the return values of job_fn."""
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(job_fn, runid, layout, idx_path, init_log_path, runtime_log_path) for runid, layout in zip(runids, layouts)]
    results = []
    for f in concurrent.futures.as_completed(futures):
        results.append(f.result())
    return results

def linear_loop(job_fn, runids, layouts, idx_path, init_log_path, runtime_log_path, workers=None): #workers is a dummy var here
    """One of two loop functions to choose from.
    Executes job_fn on each elements in runids linearly.
    Returns a list of the return values of job_fn."""
    results = []
    for runid, layout in zip(runids, layouts):
        results.append(job_fn(runid, layout, idx_path, init_log_path, runtime_log_path))
    return results

def process_batch(runids, layouts, idx_path, spe, curl=False, linear=False, workers=8):
    init_log_path = helpers.initiate_logfile('initiation', ['timestamp', 'runid'], spe=f"{spe}-")
    runtime_log_path = helpers.initiate_logfile('runtime', ['timestamp', 'runid', 'ascp_time', 'kallisto_time', 'library_layout'], spe=f"{spe}-")
    batch_start = time.time()
    # Define download modes
    loop_fn = linear_loop if linear else parallel_loop
    job_mode = curl_job if curl else ascp_job
    # Run the loop_fn in job_mode and pass other required parameters
    results = loop_fn(job_mode, runids, layouts, idx_path, init_log_path, runtime_log_path, workers)
    batch_runtime = round(time.time() - batch_start, 2)
    helpers.write_log(f"Total runtime\t{batch_runtime}\n", runtime_log_path)

__all__ = ['process_batch', 'kallisto_index']

