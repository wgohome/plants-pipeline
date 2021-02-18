# F. Postprocessing

## F1. Extracting the TPM matrices

Running this script will extract the TPM values for each species from `pipeline-data/download/kallisto-out` and compile them into a TPM matrix per species, to be stored in `pipeline-data/postprocess/tpm-matrices`. The file will be timestamped and labelled as `YYYYMMDD-HHMMSS-taxidXXXX_tpm_matrix.tsv`. In the case when there are multiple TPM matrices for a species, the latest timestamped file will be used.

Example:
```
python main/update_all_tpm_matrices.py -t 12 -n
```

Arguments:

- `-t` specified the number of threads to used to perform the tpm matrix construction for multiple species in parallel. If `-t 0` is set, then the job will be run linearly.
- `-n` is an optional argument to specify that only species for which its tpm matrices has not been constructed should be processed.
- `-m` is an optional method argument. The default method is method 1, the other methods are experimental methods to try to optimise the process but to no avail. Feel free to submit pull request if you have a more efficient way of doing it.

> TODO: head of a sample TPM matrix

## F2. Getting quality control (QC) statistics

This script is to be run twice - first to obtain the recommended cutoffs and visualisation of the resulting filter, then the user can review and edit the recommended cutoffs, and run the script again to confirm the cutoffs.

Arguments for the `get_qc_stats.py`:

- `-t` specifies the taxid(s) that QC is to be performed for. The flag accepts one or more arguments, separated by a single whitespace. `-u` should not be specified together here.
- `-u` is used only for the second step of the QC when the user has checked and edited the QC thresholds and wants to confirm it. When this is specified, `-t` does not need to be specified.

<ins>**Step 1**: Get automated QC recommendation</ins>

Get recommended QC cutoffs, its corresponsing number of samples accepted by the filter, and a jointplot visualisation for the given cutoff.

```
python postprocess/get_qc_stats.py -t XXXX YYYY ZZZZ
```

> TODO: ss of QC plot and summary table

<ins>**Step 2**: Confirm QC cutoffs</ins>

Confirm QC cutoffs for each species, to generate a json for the accepted run IDs for each species. The json file will be stored in `pipeline-data/postprocess/qc-out`

```
python postprocess/get_qc_stats.py -u
```

Sample of a qc-out json file:
```
{
  "3197": ["SRR4450256", "SRR4450257", "SRR4450265", "SRR4450266", ....],
  "3702": ["SRR5166024", "SRR5166028", "SRR5166023", "SRR5166027", "SRR5166035", "SRR5166034", ...]
}
```

## F3. Calculating F1 stats

> This part is done for the purpose of benchmarking the coexpression output for the publication. However, scripts can be modified using the available functions here and in other script files to cater to your unique needs for downstream analyses.

```
python postprocess/coexpression1.py -t XXXX YYYY ZZZZ
```

- `-t` species the taxid(s) to be processed. It accepts a list of one or more taxids, separated by a single whitespace.

> Sample output of the percentage matrix

## F4.

```
python postprocess/calc_f1_scores.py -t XXXX YYYY ZZZZ
```

- `-t` species the taxid(s) to be processed. It accepts a list of one or more taxids, separated by a single whitespace.

> Sample output of the f1 scores
