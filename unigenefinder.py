#!/usr/bin/env python

import functools
import os
import time
import psutil
import resource
from sys import stderr, argv
from unigenefinder import (
    check_programs_in_path,
    get_system_resources,
    run_transdecoder_longorfs,
    process_orfs_and_transcripts,
    run_cd_hit_est,
    run_compacta,
    process_cd_hit_clusters,
    process_compacta_clusters,
    rename_fasta_headers,
    get_fasta_lengths_and_seqs,
    create_graph_from_clusterings,
    find_longest_transcripts,
    translate_cds
)

def main():

    # Start monitoring
    process = psutil.Process(os.getpid())
    total_start_time = time.time()

    # UNIGENEFINDER BANNER
    print('\n############################################')
    print('############## UNIGENE FINDER ##############')
    print('############################################')

    # INITIALIZING
    print('\nINITIALIZING:')
    print('Checking for required programs in system path...')
    check_programs_in_path(['TransDecoder.LongOrfs', 'cd-hit-est', 'Compacta'])

    print('Getting system resource information.')
    threads_avail, mem_avail = get_system_resources()

    print('Loading configuration.')
    config_file = argv[1]
    config = {}

    with open(config_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line and not line.startswith('#'):
                key, value = line.split('=', 1)
                config[key.strip()] = value.strip()

    # Initialize variables from config dict
    trinity_assemblies = config.get('trinity_assemblies')
    trinity_map = config.get('trinity_map')
    threads = int(config.get('threads', 0))
    mem = int(config.get('mem', 0))
    c = float(config.get('c', 0.9))
    bamdir = config.get('bamdir')
    workdir = config.get('workdir')

    if threads == 0:
        threads = threads_avail
    if mem == 0:
        mem = int(round(mem_avail * 0.5))

    print(f'Creating working directory:')
    if os.path.exists(workdir):
        print(f"Error: The directory '{workdir}' already exists.")
        print('Please either remove this directory or modify your config file.')
        exit(1)
    else:
        os.mkdir(workdir)
        print(f'{workdir}')

    log_dir = os.path.join(workdir, '00_commandline_tool_logs')
    os.mkdir(log_dir)

    # TRANSDECODER
    print('\nTRANSDECODER:')
    transdecoder_start_time = time.time()
    transdecoder_dir = os.path.join(workdir, '01_transdecoder_output')
    transdecoder_log_file = os.path.join(log_dir, 'transdecoder_longorfs.log')
    print('Running TransDecoder... ', end='')
    run_transdecoder_longorfs(trinity_assemblies, trinity_map, transdecoder_dir, 50, transdecoder_log_file)
    transdecoder_time = time.time() - transdecoder_start_time
    print(f'done ({transdecoder_time:.2f} seconds).')

    transcripts_file_name = os.path.basename(trinity_assemblies)
    transdecoder_output_subdir = f'{transcripts_file_name}.transdecoder_dir'
    transdecoder_output = os.path.join(transdecoder_dir, transdecoder_output_subdir, 'longest_orfs.cds')

    print('Selecting best ORF/transcript per isoform... ', end='')
    best_per_isoform_path = os.path.join(workdir, 'BestPerIsoform')
    process_orfs_and_transcripts(transdecoder_output, trinity_assemblies, best_per_isoform_path)
    print('done.')

    # CD-HIT
    print('\nCD-HIT:')
    print('Shortening sequence names for CD-HIT... ', end='')
    cd_hit_dir = os.path.join(workdir, '02_cd-hit_output')
    os.mkdir(cd_hit_dir)
    short_to_full = rename_fasta_headers(os.path.join(trinity_assemblies),
                                         os.path.join(cd_hit_dir, 'Shortened.fasta'))
    print('done.')

    cd_hit_start_time = time.time()
    print('Running CD-HIT... ', end='')
    cd_hit_log_file = os.path.join(log_dir, 'cd_hit_est.log')
    run_cd_hit_est(threads, mem, c, os.path.join(cd_hit_dir, 'Shortened.fasta'),
                   os.path.join(cd_hit_dir, 'CDHIT_Output'), cd_hit_log_file)
    cd_hit_time = time.time() - cd_hit_start_time
    print(f'done ({cd_hit_time:.2f} seconds).')

    print('Formatting CD-HIT clusters... ', end='')
    process_cd_hit_clusters(os.path.join(cd_hit_dir, 'CDHIT_Output.clstr'), 
                            short_to_full, os.path.join(workdir, 'Clusters_CDHIT.txt'))
    print('done.')

    # COMPACTA
    print('\nCOMPACTA:')
    compacta_dir = os.path.join(workdir, '03_compacta_output')
    os.mkdir(compacta_dir)
    compacta_start_time = time.time()
    print('Running Compacta... ', end='')
    compacta_log_file = os.path.join(log_dir, 'compacta.log')
    run_compacta(threads, bamdir, os.path.join(compacta_dir, 'Compacta_Output'), compacta_log_file)
    compacta_time = time.time() - compacta_start_time
    print(f'done ({compacta_time:.2f} seconds).')

    print('Formatting Compacta clusters... ', end='')
    process_compacta_clusters(os.path.join(compacta_dir, 'Compacta_Output_clusters.txt'), 
                              os.path.join(workdir, 'Clusters_Compacta.txt'))
    print('done.')

    print('Getting sequences and  lengths... ', end='')
    transcripts_to_orf_lens, sequences_cds = get_fasta_lengths_and_seqs(
        os.path.join(workdir, 'BestPerIsoform_orfs.fasta')
    )
    transcripts_to_transcript_lens, sequences_transcripts = get_fasta_lengths_and_seqs(
        os.path.join(workdir, 'BestPerIsoform_transcripts.fasta')
    )
    sequences_protein = translate_cds(sequences_cds)
    print('done.')

    # COMBINING CLUSTERS
    print('\nCOMBINING CLUSTERS:')
    print('Generating combined cluster network...')
    G = create_graph_from_clusterings(trinity_map, os.path.join(workdir, 'Clusters_CDHIT.txt'),
                                      os.path.join(workdir, 'Clusters_Compacta.txt'))

    # FINAL STEPS
    print('\nFINAL STEPS:')
    print('Choosing best sequences per cluster... ', end='')
    clusters, best = find_longest_transcripts(G, transcripts_to_orf_lens, transcripts_to_transcript_lens)
    print('done.')

    results_dir = os.path.join(workdir, '04_results')
    os.mkdir(results_dir)

    print('Writing output files... ', end='')
    with open(os.path.join(results_dir, 'UniGeneFinder_Clusters.txt'), 'w') as outfile:
        for cluster, seqIDs in clusters.items():
            for seqID in seqIDs:
                print(cluster, seqID, sep='\t', file=outfile)

    with open(os.path.join(results_dir, 'UniGeneFinder_Best_Transcripts.fasta'), 'w') as outfile:
        for cluster, seqID in best.items():
            if seqID in sequences_transcripts:
                outfile.write(f'>{cluster}\n{sequences_transcripts[seqID]}\n')

    with open(os.path.join(results_dir, 'UniGeneFinder_Best_CDS.fasta'), 'w') as outfile:
        for cluster, seqID in best.items():
            if seqID in sequences_cds:
                outfile.write(f'>{cluster}\n{sequences_cds[seqID]}\n')

    with open(os.path.join(results_dir, 'UniGeneFinder_Best_Proteins.fasta'), 'w') as outfile:
        for cluster, seqID in best.items():
            if seqID in sequences_protein:
                outfile.write(f'>{cluster}\n{sequences_protein[seqID]}\n')

    print('done.')

    # ALL ACTIONS COMPLETE
    print('\nALL ACTIONS COMPLETE.')

    peak_memory = process.memory_info().rss
    print(f'Peak memory usage: {peak_memory / 1024 / 1024:.2f} MiB')

    total_time = time.time() - total_start_time
    print(f'Total runtime: {total_time:.2f} seconds')

if __name__ == '__main__':

    # Create a partial print function with sys.stderr as the default file and flush=True by default
    print = functools.partial(print, file=stderr, flush=True)

    main()
