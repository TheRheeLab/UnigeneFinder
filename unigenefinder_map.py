#!/usr/bin/env python

import os
import argparse
import functools
import multiprocessing
from sys import exit, stderr
from unigenefinder import check_programs_in_path
from unigenefinder.map import (
    build_bowtie2_index,
    map_single_end_reads_bowtie2,
    map_paired_end_reads_bowtie2,
)

def main():
    parser = argparse.ArgumentParser(
        description='Build Bowtie2 index and run Bowtie2 on FASTQ files in a directory.'
    )
    
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '--ref', required=True, help='Input reference FASTA file'
    )
    required.add_argument(
        '--reads', required=True, help='Directory containing FASTQ files'
    )
    required.add_argument(
        '--type', required=True, choices=['single', 'paired'],
        help='Type of reads: "single" or "paired"'
    )

    parser.add_argument(
        '--threads', type=int, default=multiprocessing.cpu_count(),
        help='Number of threads (default: number of CPUs detected)'
    )
    parser.add_argument(
        '--outdir', default=os.path.join(os.getcwd(), 'unigenefinder_mapping'),
        help='Master output directory (default: ./unigenefinder_mapping)'
    )

    args = parser.parse_args()

    fasta = args.ref
    fastq_dir = args.reads
    read_type = args.type
    threads = args.threads
    fasta_base = os.path.basename(fasta).split('.')[0]
    outdir = args.outdir

    # UNIGENEFINDER BANNER
    print('\n####################################################')
    print('############## UNIGENE FINDER MAPPING ##############')
    print('####################################################')

    # INITIALIZING
    print('\nINITIALIZING:\n')

    # Check that necessary programs are available in the system's PATH
    print('Checking for required programs in system path...')
    check_programs_in_path(['bowtie2-build',
                            'bowtie2',
                            'samtools'
                            ])

    # Set up directories under the master output directory
    print(f'Creating output directory "{outdir}"')

    if os.path.exists(outdir):
        print(f'Error: Output directory "{outdir}" already exists. Exiting.')
        exit(1)
    os.makedirs(outdir)
    
    index_dir = os.path.join(outdir, f'{fasta_base}_index')
    bam_dir = os.path.join(outdir, f'{fasta_base}_BAM')
    log_dir = os.path.join(outdir, f'{fasta_base}_logs')

    os.makedirs(index_dir)
    os.makedirs(bam_dir)
    os.makedirs(log_dir)

    # BUILDING INDEX
    print('\nBUILDING INDEX:\n')
    build_bowtie2_index(fasta, index_dir, threads)

    # MAPPING READS
    print('\nMAPPING READS:\n')
    if read_type == 'single':
        map_single_end_reads_bowtie2(
            fastq_dir, index_dir, bam_dir, log_dir, fasta_base, threads
        )
    elif read_type == 'paired':
        map_paired_end_reads_bowtie2(
            fastq_dir, index_dir, bam_dir, log_dir, fasta_base, threads
        )
    
    print('\nAll actions completed.\n')

if __name__ == '__main__':

    # Create a partial print function with sys.stderr as the default file and flush=True by default
    print = functools.partial(print, file=stderr, flush=True)

    main()
