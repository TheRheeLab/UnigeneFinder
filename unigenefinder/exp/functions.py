import os
import subprocess
import time
import pandas as pd
import numpy as np
from sys import exit, stderr
import functools
from typing import List, Dict

# Create a partial print function with sys.stderr as the default file and flush=True by default
print = functools.partial(print, file=stderr, flush=True)

def run_rsem_prepare_reference(
    ref_fasta: str, 
    map_file: str, 
    out_prefix: str, 
    threads: int, 
    log_file: str
) -> None:
    """
    Run the rsem-prepare-reference command with the specified arguments.

    Parameters:
    - ref_fasta (str): Path to the reference FASTA file.
    - map_file (str): Path to the transcript-to-gene map file.
    - out_prefix (str): Output prefix for the RSEM reference files.
    - threads (int): Number of threads to use for the process.
    - log_file (str): Path to the log file where output will be written.
    """
    cmd = [
        'rsem-prepare-reference', 
        '--bowtie2', 
        '--num-threads', str(threads), 
        '--transcript-to-gene-map', map_file, 
        ref_fasta, 
        out_prefix
    ]
    
    print(f'Running command: {" ".join(cmd)}')
    start_time = time.time()
    
    with open(log_file, 'w') as log:
        result = subprocess.run(cmd, stdout=log, stderr=log)
    
    end_time = time.time()
    runtime = end_time - start_time
    print(f'Command completed in {runtime:.2f} seconds.\n')
    
    if result.returncode != 0:
        print(
            f'Error: RSEM reference preparation failed; check output in log file {log_file}.'
        )
        exit(1)

def run_rsem_scan_paired_end_reads(
    threads: int, 
    input_bam: str, 
    output_bam: str, 
    log_file: str
) -> None:
    """
    Run the rsem-scan-for-paired-end-reads command with the specified arguments.

    Parameters:
    - threads (int): Number of threads to use for the process.
    - input_bam (str): Path to the input BAM file.
    - output_bam (str): Path to the output BAM file.
    - log_file (str): Path to the log file where output will be written.
    """
    cmd = [
        'rsem-scan-for-paired-end-reads', 
        str(threads), 
        input_bam, 
        output_bam
    ]
    
    print(f'Running command: {" ".join(cmd)}')
    start_time = time.time()
    
    with open(log_file, 'w') as log:
        result = subprocess.run(cmd, stdout=log, stderr=log)
    
    end_time = time.time()
    runtime = end_time - start_time
    print(f'Command completed in {runtime:.2f} seconds.\n')
    
    if result.returncode != 0:
        print(
            f'Error: RSEM paired-end reads scanning failed; check output in log file {log_file}.'
        )
        exit(1)

def run_rsem_calculate_expression(
    threads: int, 
    bam_file: str, 
    ref_name: str, 
    out_prefix: str, 
    paired: bool, 
    log_file: str
) -> None:
    """
    Run the rsem-calculate-expression command with the specified arguments.

    Parameters:
    - threads (int): Number of threads to use for the process.
    - bam_file (str): Path to the BAM file to be used.
    - ref_name (str): Name of the reference used for RSEM.
    - out_prefix (str): Output prefix for the RSEM expression files.
    - paired (bool): Indicates if the input data is paired-end.
    - log_file (str): Path to the log file where output will be written.
    """
    cmd = [
        'rsem-calculate-expression', 
        '--num-threads', str(threads), 
        '--bam'
    ]
    
    if paired:
        cmd.append('--paired-end')
    
    cmd.extend([bam_file, ref_name, out_prefix])
    
    print(f'Running command: {" ".join(cmd)}')
    start_time = time.time()
    
    with open(log_file, 'w') as log:
        result = subprocess.run(cmd, stdout=log, stderr=log)
    
    end_time = time.time()
    runtime = end_time - start_time
    print(f'Command completed in {runtime:.2f} seconds.\n')
    
    if result.returncode != 0:
        print(
            f'Error: RSEM expression calculation failed; check output in log file {log_file}.'
        )
        exit(1)

def compile_expression_data(
    results_folder: str, 
    tpm_output_csv_path: str, 
    counts_output_csv_path: str, 
    mode: str = 'gene'
) -> None:
    """
    Compile data frames from all .genes.results or .isoforms.results files in the specified results directory,
    and write them to separate CSV files for TPM and counts.

    Parameters:
    - results_folder (str): Path to the directory containing RSEM result files.
    - tpm_output_csv_path (str): Path to save the compiled TPM data as a CSV file.
    - counts_output_csv_path (str): Path to save the compiled counts data as a CSV file.
    - mode (str, optional): Specifies whether to compile 'gene' or 'isoform' data. Defaults to 'gene'.

    Raises:
    - ValueError: If the mode is not 'gene' or 'isoform'.
    """
    if mode == 'gene':
        file_extension = '.genes.results'
        id_column = 'gene_id'
        output_id_column = 'GeneID'
    elif mode == 'isoform':
        file_extension = '.isoforms.results'
        id_column = 'transcript_id'
        output_id_column = 'TranscriptID'
    else:
        raise ValueError('Invalid mode. Use "gene" or "isoform".')

    files_in_directory = os.listdir(results_folder)
    
    results_files = [
        os.path.join(results_folder, file)
        for file in files_in_directory
        if file.endswith(file_extension)
    ]
    
    tpm_data: Dict[str, List] = {
        output_id_column: [],
        'Length': [],
        'EffectiveLengthMean': [],
        'EffectiveLengthStDev': []
    }
    
    counts_data: Dict[str, List] = {
        output_id_column: [],
        'Length': [],
        'EffectiveLengthMean': [],
        'EffectiveLengthStDev': []
    }
    
    effective_length_data: Dict[str, List[float]] = {}
    
    for file in sorted(results_files):
        df = pd.read_csv(file, sep='\t')
        
        sample_name = os.path.basename(file).replace(file_extension, '')
        tpm_data[sample_name] = df['TPM'].values
        counts_data[sample_name] = df['expected_count'].round().astype(int).values
        
        if len(tpm_data[output_id_column]) == 0:
            tpm_data[output_id_column] = df[id_column].values
            tpm_data['Length'] = df['length'].values
            counts_data[output_id_column] = df[id_column].values
            counts_data['Length'] = df['length'].values
        
        for transcript_id, effective_length in zip(
            df[id_column], df['effective_length']
        ):
            if transcript_id not in effective_length_data:
                effective_length_data[transcript_id] = []
            effective_length_data[transcript_id].append(effective_length)
    
    for transcript_id in tpm_data[output_id_column]:
        lengths = effective_length_data[transcript_id]
        mean_length = np.mean(lengths)
        
        if len(lengths) >= 3:
            std_length = np.std(lengths, ddof=1)
            tpm_data['EffectiveLengthMean'].append(mean_length)
            tpm_data['EffectiveLengthStDev'].append(std_length)
            counts_data['EffectiveLengthMean'].append(mean_length)
            counts_data['EffectiveLengthStDev'].append(std_length)
        else:
            tpm_data['EffectiveLengthMean'].append(mean_length)
            tpm_data['EffectiveLengthStDev'].append(np.nan)
            counts_data['EffectiveLengthMean'].append(mean_length)
            counts_data['EffectiveLengthStDev'].append(np.nan)
        
    tpm_final_df = pd.DataFrame(tpm_data)
    counts_final_df = pd.DataFrame(counts_data)
    
    tpm_final_df.to_csv(tpm_output_csv_path, index=False)
    counts_final_df.to_csv(counts_output_csv_path, index=False)
    print(f'Compiled TPM data written to {tpm_output_csv_path}')
    print(f'Compiled counts data written to {counts_output_csv_path}')
