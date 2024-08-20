#!/usr/bin/env python

import sys
from sys import stderr, argv
import os
from unigenefinder import (
    check_programs_in_path,
    get_system_resources,
    run_transdecoder_longorfs,
    process_orfs_and_transcripts,
    rename_fasta_headers,
    run_cd_hit_est,
    process_cd_hit_clusters,
    run_compacta,
    process_compacta_clusters,
    fasta_header_to_length,
    create_graph_from_clusterings,
    find_longest_transcripts
)

def main():
    """Main function to run the UniGeneFinder pipeline."""
    
    # Check for required programs in system path
    print('Checking for required programs in system path...', file=stderr)
    check_programs_in_path(['TransDecoder.LongOrfs', 'cd-hit-est', 'Compacta'])

    # Get system resource information
    print('Getting system resource information...', file=stderr)
    threads_avail, mem_avail = get_system_resources()

    # Load configuration from the provided config file
    print('Loading configuration...', file=stderr)
    config_file = sys.argv[1]
    config = {}

    with open(config_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line and not line.startswith('#'):
                key, value = line.split('=', 1)
                config[key.strip()] = value.strip()

    trinity_assemblies = config.get('trinity_assemblies')
    trinity_map = config.get('trinity_map')
    threads = int(config.get('threads', 0))
    mem = int(config.get('mem', 0))
    c = float(config.get('c', 0.9))
    bamdir = config.get('bamdir')
    workdir = config.get('workdir')

    # Use available resources if not specified in the config
    if threads == 0:
        threads = threads_avail
    if mem == 0:
        mem = int(round(mem_avail * 0.5))

    # Create the working directory
    print(f'Creating working directory {workdir}...')
    os.mkdir(workdir)

    # Define log file paths for each command line tool
    transdecoder_log = os.path.join(workdir, 'transdecoder.log')
    cdhit_log = os.path.join(workdir, 'cdhit.log')
    compacta_log = os.path.join(workdir, 'compacta.log')

    # Run TransDecoder to find long open reading frames (ORFs)
    print('Running TransDecoder...', file=stderr)
    transdecoder_dir = os.path.join(workdir, 'TransDecoder_Output')
    run_transdecoder_longorfs(trinity_assemblies, trinity_map, transdecoder_dir, 50, transdecoder_log)

    # Locate the TransDecoder output file
    transcripts_file_name = os.path.basename(trinity_assemblies)
    transdecoder_output_subdir = f"{transcripts_file_name}.transdecoder_dir"
    transdecoder_output = os.path.join(transdecoder_dir, transdecoder_output_subdir, 'longest_orfs.cds')

    # Process ORFs and transcripts
    print('Selecting best ORF/transcript per isoform...', file=stderr)
    best_per_isoform_path = os.path.join(workdir, 'BestPerIsoform')
    process_orfs_and_transcripts(transdecoder_output, trinity_assemblies, best_per_isoform_path)

    # Shorten sequence names for CD-HIT
    print('Shortening sequence names for CD-HIT...', file=stderr)
    short_to_full = rename_fasta_headers(os.path.join(workdir, 'BestPerIsoform_orfs.fasta'),
                                         os.path.join(workdir, 'Shortened.fasta'))

    # Run CD-HIT to cluster sequences
    print('Running CD-HIT...', file=stderr)
    run_cd_hit_est(threads, mem, c, os.path.join(workdir, 'Shortened.fasta'),
                   os.path.join(workdir, 'CDHIT_Output'), cdhit_log)

    # Process CD-HIT clusters
    print('Formatting CD-HIT clusters...', file=stderr)
    process_cd_hit_clusters(os.path.join(workdir, 'CDHIT_Output.clstr'), 
                            short_to_full, os.path.join(workdir, 'Clusters_CDHIT.txt'))

    # Run Compacta to cluster sequences based on coverage
    print('Running Compacta...', file=stderr)
    run_compacta(threads, bamdir, os.path.join(workdir, 'Compacta_Output'), compacta_log)

    # Process Compacta clusters
    print('Formatting Compacta clusters...', file=stderr)
    process_compacta_clusters(os.path.join(workdir, 'Compacta_Output_clusters.txt'), 
                              os.path.join(workdir, 'Clusters_Compacta.txt'))

    # Get sequence lengths for ORFs and transcripts
    print('Getting sequence lengths...', file=stderr)
    transcripts_to_orf_lens = fasta_header_to_length(os.path.join(workdir, 'BestPerIsoform_orfs.fasta'))
    transcripts_to_transcript_lens = fasta_header_to_length(os.path.join(workdir, 'BestPerIsoform_transcripts.fasta'))

    # Generate a combined cluster network
    print('Generating combined cluster network...', file=stderr)
    G = create_graph_from_clusterings(trinity_map, os.path.join(workdir, 'Clusters_CDHIT.txt'),
                                      os.path.join(workdir, 'Clusters_Compacta.txt'))

    # Select longest transcripts per cluster
    print('Choosing best sequences per cluster...', file=stderr)
    clusters, best = find_longest_transcripts(G, transcripts_to_orf_lens, transcripts_to_transcript_lens)

    # Write output files
    print('Writing output files...', file=stderr)
    with open('UniGeneFinder_Clusters.txt', 'w') as outfile:
        for cluster, seqIDs in clusters.items():
            print(cluster, ','.join(seqIDs), sep='\t', file=outfile)
        
    with open('UniGeneFinder_Best.txt', 'w') as outfile:
        for cluster, seqID in best.items():
            print(cluster, seqID, sep='\t', file=outfile)
        
    print('All actions complete.', file=stderr)

if __name__ == "__main__":
    main()