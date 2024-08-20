#!/usr/bin/env python

import os
import argparse
import multiprocessing
from unigenefinder.map import check_programs, build_hisat2_index, map_single_end_reads, map_paired_end_reads

def main():
    """Main function to build HISAT2 index and run HISAT2 on FASTQ files in a directory."""
    
    # Argument parser setup
    parser = argparse.ArgumentParser(description="Build HISAT2 index and run HISAT2 on FASTQ files in a directory.")
    
    # Required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument("--ref", required=True, help="Input reference FASTA file")
    required.add_argument("--reads", required=True, help="Directory containing FASTQ files")
    required.add_argument("--type", required=True, choices=['single', 'paired'], help="Type of reads: 'single' or 'paired'")

    # Optional arguments
    parser.add_argument("--threads", type=int, default=multiprocessing.cpu_count(), help="Number of threads (default: number of CPUs detected)")
    parser.add_argument("--index_dir", help="Index directory name (default: X_index where X is the input FASTA file name without extension)")
    parser.add_argument("--output_dir", help="Output directory (default: X_BAM where X is the input FASTA file name without extension)")
    parser.add_argument("--log_dir", help="Log output directory (default: X_logs where X is the input FASTA file name without extension)")

    # Parse the arguments
    args = parser.parse_args()

    # Extract arguments
    fasta = args.ref
    fastq_dir = args.reads
    read_type = args.type
    threads = args.threads
    fasta_base = os.path.basename(fasta).split('.')[0]
    index_dir = args.index_dir or f"{fasta_base}_index"
    output_dir = args.output_dir or f"{fasta_base}_BAM"
    log_dir = args.log_dir or f"{fasta_base}_logs"

    # Create necessary directories
    os.makedirs(index_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    # Check for required programs
    check_programs()

    # Build HISAT2 index
    build_hisat2_index(fasta, index_dir, threads)

    # Map reads
    if read_type == 'single':
        map_single_end_reads(fastq_dir, index_dir, output_dir, log_dir, fasta_base, threads)
    elif read_type == 'paired':
        map_paired_end_reads(fastq_dir, index_dir, output_dir, log_dir, fasta_base, threads)

if __name__ == "__main__":
    main()