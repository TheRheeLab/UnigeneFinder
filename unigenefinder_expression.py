#!/usr/bin/env python

import os
import argparse
import functools
import multiprocessing
from sys import stderr
from unigenefinder import check_programs_in_path
from unigenefinder.exp import (
    run_rsem_prepare_reference,
    run_rsem_scan_paired_end_reads,
    run_rsem_calculate_expression,
    compile_expression_data
)

def parse_args():
    parser = argparse.ArgumentParser(description='Run RSEM commands on a folder of BAM files representing mapped reads.')
    
    # Required arguments
    parser.add_argument('--ref', required=True, help='Path to FASTA file for reference sequences')
    parser.add_argument('--map', required=True, help='Path to a transcript-to-gene map file')
    parser.add_argument('--bamdir', required=True, help='Path to a folder of BAM files representing mapped reads')
    parser.add_argument('--type', required=True, choices=['single', 'paired'], help='Type of reads: "single" or "paired"')
    
    # Optional arguments
    parser.add_argument('--threads', type=int, default=multiprocessing.cpu_count(), help='Number of threads to use (default: number of available cores)')
    parser.add_argument('--out', type=str, default='unigenefinder_expression', help='Output folder (default: "unigenefinder_expression")')
    
    return parser.parse_args()

def main():

    args = parse_args()

    # UNIGENEFINDER BANNER
    print('\n#######################################################')
    print('############## UNIGENE FINDER EXPRESSION ##############')
    print('#######################################################')

    # INITIALIZING
    print('\nINITIALIZING:\n')

    # Check that necessary programs are available in the system's PATH
    print('Checking for required programs in system path...')
    check_programs_in_path(['bowtie2-build',
                            'rsem-prepare-reference',
                            'rsem-scan-for-paired-end-reads',
                            'rsem-calculate-expression'])

    # Check if the output directory exists; if it does, exit with an error
    print('Creating output directories... ', end='')
    if os.path.exists(args.out):
        print(f'Error: Output directory "{args.out}" already exists. Please specify a different directory or remove the existing one.')
        exit(1)
    else:
        # Create the output directory and necessary subdirectories
        os.makedirs(os.path.join(args.out, 'reference'))
        os.makedirs(os.path.join(args.out, 'logs'))
        os.makedirs(os.path.join(args.out, 'results'))
        if args.type == 'paired':
                os.makedirs(os.path.join(args.out, 'converted_bam'))
    print('done.')

    # PREPARING RSEM REFERENCE
    print('\nPREPARING RSEM REFERENCE:\n')

    # Run the rsem-prepare-reference command
    reference_log_file = os.path.join(args.out, 'logs', '00_reference_preparation_log.txt')
    reference_out_prefix = os.path.join(args.out, 'reference', 'rsem_reference')
    run_rsem_prepare_reference(args.ref, args.map, reference_out_prefix, args.threads, reference_log_file)

    # CONVERTING PAIRED-END BAM FILES
    print('CONVERTING PAIRED-END BAM FILES:\n')

    # If the read type is 'paired', process each BAM file in the supplied directory
    if args.type == 'paired':
        bam_files = [f for f in os.listdir(args.bamdir) if f.lower().endswith('.bam')]
        
        for bam_file in bam_files:
            input_bam = os.path.join(args.bamdir, bam_file)
            output_bam = os.path.join(args.out, 'converted_bam', bam_file)
            log_file = os.path.join(args.out, 'logs', f'{bam_file}_scan_log.txt')
            
            # Run the rsem-scan-for-paired-end-reads command
            run_rsem_scan_paired_end_reads(args.threads, input_bam, output_bam, log_file)

    # CALCULATING EXPRESSION
    print('CALCULATING EXPRESSION:\n')

    # Set bam_dir based on whether the reads are paired-end or single-end
    bam_dir = os.path.join(args.out, 'converted_bam') if args.type == 'paired' else args.bamdir

    # Loop over all BAM files in the determined bam_dir, calculating expression
    bam_files = [f for f in os.listdir(bam_dir) if f.lower().endswith('.bam')]

    for bam_file in bam_files:
        input_bam = os.path.join(bam_dir, bam_file)
        bam_name = os.path.splitext(bam_file)[0]  # Get the BAM file name without the extension
        out_prefix = os.path.join(args.out, 'results', bam_name)
        log_file = os.path.join(args.out, 'logs', f'{bam_name}_expression_log.txt')
        
        # Run the rsem-calculate-expression command
        run_rsem_calculate_expression(args.threads, input_bam, reference_out_prefix, out_prefix, args.type == 'paired', log_file)

    # COMPILING EXPRESSION CSV TABLES
    print('COMPILING EXPRESSION CSV TABLES:\n')

    # Compile gene expression data into CSV files
    print('Per-Gene Expression Data:')
    compile_expression_data(
        results_folder=os.path.join(args.out, 'results'),
        tpm_output_csv_path=os.path.join(args.out, 'compiled_tpm_genes.csv'),
        counts_output_csv_path=os.path.join(args.out, 'compiled_counts_genes.csv'),
        mode='gene'
    )

    # Compile isoform expression data into CSV files
    print('\nPer-Isoform Expression Data:')
    compile_expression_data(
        results_folder=os.path.join(args.out, 'results'),
        tpm_output_csv_path=os.path.join(args.out, 'compiled_tpm_isoforms.csv'),
        counts_output_csv_path=os.path.join(args.out, 'compiled_counts_isoforms.csv'),
        mode='isoform'
    )

    print('\nAll actions completed.\n')

if __name__ == '__main__':

    # Create a partial print function with sys.stderr as the default file and flush=True by default
    print = functools.partial(print, file=stderr, flush=True)

    main()
