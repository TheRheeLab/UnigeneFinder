import os
import functools
import subprocess
import time
import pandas as pd
import numpy as np
from sys import exit, stderr
from typing import List, Dict
from concurrent.futures import ProcessPoolExecutor

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
    bam_file: str, 
    ref_name: str, 
    out_prefix: str, 
    paired: bool, 
    log_file: str
) -> None:
    """
    Run the rsem-calculate-expression command with the specified arguments,
    without specifying an internal thread count.

    Parameters:
    - bam_file (str): Path to the BAM file to be used.
    - ref_name (str): Name of the reference used for RSEM.
    - out_prefix (str): Output prefix for the RSEM expression files.
    - paired (bool): Indicates if the input data is paired-end.
    - log_file (str): Path to the log file where output will be written.
    """
    cmd = [
        'rsem-calculate-expression',
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
        print(f'Error: RSEM expression calculation failed; check output in log file {log_file}.')
        exit(1)

def process_bam_file(
    bam_file: str,
    bam_dir: str,
    reference_out_prefix: str,
    results_dir: str,
    logs_dir: str,
    paired: bool
) -> None:
    """
    Process a single BAM file by constructing the input and output paths,
    then calling run_rsem_calculate_expression to calculate expression.

    Parameters:
      - bam_file (str): The name of the BAM file.
      - bam_dir (str): Directory containing the BAM files.
      - reference_out_prefix (str): Path prefix for the prepared RSEM reference.
      - results_dir (str): Directory where result files will be written.
      - logs_dir (str): Directory where log files will be stored.
      - paired (bool): True if the reads are paired-end; False if single-end.
    """
    # Construct the full path to the BAM file
    bam_path = os.path.join(bam_dir, bam_file)
    # Remove file extension to obtain the base name for outputs
    bam_name = os.path.splitext(bam_file)[0]
    # Create the output prefix for this BAM file
    out_prefix = os.path.join(results_dir, bam_name)
    # Construct the log file path
    log_file = os.path.join(logs_dir, f'{bam_name}_expression_log.txt')

    # Call the RSEM expression calculation function
    run_rsem_calculate_expression(
        bam_file=bam_path,
        ref_name=reference_out_prefix,
        out_prefix=out_prefix,
        paired=paired,
        log_file=log_file
    )

def process_all_bam_files(
    bam_dir: str,
    reference_out_prefix: str,
    results_dir: str,
    logs_dir: str,
    paired: bool,
    threads: int
) -> None:
    """
    Process all BAM files in the specified directory in parallel using ProcessPoolExecutor.

    For each BAM file in bam_dir (files ending with .bam), this function constructs
    the full input path, output prefix, and log file path, and then calls
    run_rsem_calculate_expression (via process_bam_file) for that file.

    Parameters:
      - bam_dir (str): Directory containing the BAM files.
      - reference_out_prefix (str): Path prefix for the prepared RSEM reference.
      - results_dir (str): Directory where result files will be written.
      - logs_dir (str): Directory where log files will be stored.
      - paired (bool): True if the reads are paired-end; False for single-end.
      - threads (int): Maximum number of parallel processes.
    """
    # List all BAM files (case-insensitive)
    bam_files = [f for f in os.listdir(bam_dir) if f.lower().endswith('.bam')]

    # Bind constant parameters to process_bam_file using functools.partial.
    bound_process = functools.partial(
        process_bam_file,
        bam_dir=bam_dir,
        reference_out_prefix=reference_out_prefix,
        results_dir=results_dir,
        logs_dir=logs_dir,
        paired=paired
    )

    # Execute in parallel using ProcessPoolExecutor.
    with ProcessPoolExecutor(max_workers=threads) as executor:
        executor.map(bound_process, bam_files)

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
