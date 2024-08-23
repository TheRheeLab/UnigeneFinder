from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
import os
import re
from collections import defaultdict
import sys
import networkx as nx
from shutil import which
import psutil

###### FUNCTION TO GET SYSTEM RESOURCE INFORMATION

def get_system_resources() -> tuple:
    """Returns the number of threads available and the amount of memory in megabytes."""
    # Get the number of CPU threads
    threads = os.cpu_count()

    # Get the total memory in bytes and convert to megabytes
    memory_bytes = psutil.virtual_memory().total
    memory_megabytes = round(memory_bytes / (1024 * 1024))

    return threads, memory_megabytes

###### FUNCTION TO CHECK IF ALL REQUIRED PROGRAMS ARE IN SYSTEM PATH ######

def check_programs_in_path(programs: list) -> None:
    """Checks if specified programs are in the system path and exits if any are not found."""
    for program in programs:
        sys.stderr.write(f'Checking if {program} is in the system path... ')
        if which(program) is None:
            sys.stderr.write(f'Error: {program} not found in the system path.\n')
            sys.exit(1)
        else:
            sys.stderr.write(f'{program} found.\n')

###### FUNCTION TO LOAD CONFIGURATION OPTIONS ######

def load_config(config_file: str) -> None:
    """Loads configuration variables from a file into global variables."""
    global trinity_assemblies, trinity_map, threads, mem, workdir
    config = {}
    
    with open(config_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line and not line.startswith('#'):
                key, value = line.split('=', 1)
                config[key.strip()] = value.strip()
    
    print(config)
    
    trinity_assemblies = config.get('trinity_assemblies')
    trinity_map = config.get('trinity_map')
    threads = int(config.get('threads', 0))
    mem = int(config.get('mem', 0))
    workdir = config.get('workdir')

###### FUNCTION(S) FOR RUNNING TRANSDECODER ######

def run_transdecoder_longorfs(
    fasta_file, map_file, out_dir, min_len, log_file
):
    """
    Run TransDecoder.LongOrfs with the given parameters and log output to a file.
    """
    
    # Construct the command
    cmd = [
        'TransDecoder.LongOrfs',
        '-t', fasta_file,
        '-m', str(min_len),
        '--gene_trans_map', map_file,
        '--output_dir', out_dir
    ]
    
    # Open the log file for writing
    with open(log_file, 'w') as log:
        # Execute the command, routing stdout and stderr to the log file
        subprocess.run(cmd, stdout=log, stderr=log, check=True)

###### FUNCTION(S) FOR GETTING LONGEST ORF/TRANSCRIPT FROM TRANSDECODER OUTPUT ######        

def group_fasta_by_pattern(fasta_path):
    """
    Group sequences by pattern, removing everything following the first whitespace in headers.
    """
    pattern = re.compile(r'TRINITY_DN[0-9]+_c[0-9]+_g[0-9]+_i[0-9]+')
    grouped_sequences = defaultdict(dict)
    
    for record in SeqIO.parse(fasta_path, 'fasta'):
        header = record.id.split()[0]
        match = pattern.search(header)
        if match:
            key = match.group()
            grouped_sequences[key][header] = str(record.seq)
    
    return grouped_sequences

def find_longest_sequence(sequences):
    """
    Helper function to find the longest sequences from a dictionary of sequences.
    """
    sorted_sequences = sorted(
        sequences.items(), key=lambda item: len(item[1]), reverse=True
    )
    longest_length = len(sorted_sequences[0][1])
    longest_sequences = [
        item for item in sorted_sequences if len(item[1]) == longest_length
    ]
    return longest_sequences

def find_associated_transcript(header, transcript_groups, transcript_id):
    """
    Helper function to find the associated transcript for a given ORF header.
    """
    return next(iter(transcript_groups[transcript_id].items()))

def process_orf_group(transcript_id, orf_dict, transcript_groups):
    """
    Helper function to process a group of ORFs for a single transcript ID.
    """
    new_orfs_dict = {}
    new_transcripts_dict = {}
    
    longest_orfs = find_longest_sequence(orf_dict)
    longest_orf_header, longest_orf_seq = longest_orfs[0]
    new_orfs_dict[transcript_id] = longest_orf_seq

    if len(longest_orfs) == 1:
        transcript_header, transcript_seq = find_associated_transcript(
            longest_orf_header, transcript_groups, transcript_id
        )
        new_transcripts_dict[transcript_id] = transcript_seq
    else:
        associated_transcripts = []
        for header, seq in longest_orfs:
            transcript_header, transcript_seq = find_associated_transcript(
                header, transcript_groups, transcript_id
            )
            associated_transcripts.append((transcript_header, transcript_seq))
        longest_transcript_header, longest_transcript_seq = max(
            associated_transcripts, key=lambda item: len(item[1])
        )
        new_transcripts_dict[transcript_id] = longest_transcript_seq

    return new_orfs_dict, new_transcripts_dict

def process_orfs_and_transcripts(
    orf_fasta, transcript_fasta, output_prefix
):
    """
    Process ORF and transcript FASTA files and write output files with selected sequences.
    """
    # Group sequences by pattern in both ORF and transcript FASTA files
    orf_groups = group_fasta_by_pattern(orf_fasta)
    transcript_groups = group_fasta_by_pattern(transcript_fasta)

    new_orfs_dict = {}
    new_transcripts_dict = {}

    for transcript_id, orf_dict in orf_groups.items():
        orfs, transcripts = process_orf_group(
            transcript_id, orf_dict, transcript_groups
        )
        new_orfs_dict.update(orfs)
        new_transcripts_dict.update(transcripts)

    # Write the new ORFs and transcripts to output files
    with open(f'{output_prefix}_orfs.fasta', 'w') as orf_out:
        for header, seq in new_orfs_dict.items():
            orf_out.write(f'>{header}\n{seq}\n')

    with open(f'{output_prefix}_transcripts.fasta', 'w') as transcript_out:
        for header, seq in new_transcripts_dict.items():
            transcript_out.write(f'>{header}\n{seq}\n')

###### FUNCTION TO SHORTEN FASTA NAMES FOR CD-HIT ######

def rename_fasta_headers(input_fasta, output_fasta):
    """
    Rename FASTA headers to 1, 2, 3, etc. and save to a new file.
    """
    header_mapping = {}
    with open(output_fasta, 'w') as out_f:
        for idx, record in enumerate(SeqIO.parse(input_fasta, 'fasta'), start=1):
            original_header = record.id.split()[0]
            new_header = str(idx)
            header_mapping[new_header] = original_header
            record.id = new_header
            record.description = ''
            SeqIO.write(record, out_f, 'fasta')
    
    return header_mapping

###### FUNCTION TO RUN CD-HIT_EST ######

def run_cd_hit_est(
    threads: int, mem: int, c: float, fasta: str, output: str, log_file: str
) -> None:
    """
    Run cd-hit-est with the given parameters and log output to a file.
    """
    cmd = [
        'cd-hit-est',
        '-T', str(threads),
        '-M', str(mem),
        '-i', fasta,
        '-o', output,
        '-c', str(c),
        '-n', '8',
        '-G', '0',
        '-aL', '0.5',
        '-aS', '0.75'
    ]
    
    # Open the log file for writing
    with open(log_file, 'w') as log:
        # Execute the command, routing stdout and stderr to the log file
        subprocess.run(cmd, stdout=log, stderr=log, check=True)

###### FUNCTION TO PROCESS CD-HIT OUTPUT INTO A COMMON CLUSTERING FORMAT ######    
    
def process_cd_hit_clusters(
    cluster_info_path, name_mapping, output_filename
):
    """
    Process CD-HIT cluster info and write the mapping of cluster IDs to full transcript IDs.
    """
    cluster_dict = defaultdict(list)
    cluster_id = None
    
    with open(cluster_info_path, 'r') as cluster_file:
        for line in cluster_file:
            if line.startswith('>Cluster'):
                cluster_id = line.strip().replace('>', '').replace(' ', '')
            else:
                match = re.search(r'>\d+', line)
                if match:
                    short_id = match.group().replace('>', '')
                    if short_id in name_mapping:
                        full_id = name_mapping[short_id]
                        cluster_dict[cluster_id].append(full_id)
    
    with open(output_filename, 'w') as out_file:
        for cluster_id, transcript_ids in cluster_dict.items():
            for transcript_id in transcript_ids:
                out_file.write(f'{cluster_id}\t{transcript_id}\n')
                
###### FUNCTION TO RUN COMPACTA ######                

def run_compacta(
    threads, bam_dir, output, log_file
):
    """
    Run Compacta with the specified parameters and log output to a file.
    """
    # Get list of .bam files in the directory
    bam_files_list = [
        os.path.join(bam_dir, f)
        for f in os.listdir(bam_dir) if f.endswith('.bam')
    ]
    n = len(bam_files_list)
    bam_files = ','.join(bam_files_list)
    
    # Construct the command
    cmd = [
        'Compacta',
        '-t', str(threads),
        '-n', str(n),
        '-b', bam_files,
        '-o', str(output)
    ]
    
    # Open the log file for writing
    with open(log_file, 'w') as log:
        # Execute the command, routing stdout and stderr to the log file
        subprocess.run(cmd, stdout=log, stderr=log, check=True)

###### FUNCTION TO PROCESS COMPACTA OUTPUT INTO A COMMON CLUSTERING FORMAT ######        
        
def process_compacta_clusters(
    compacta_cluster_path, output_file_path
):
    """
    Process Compacta's output cluster text file and generate a general clustering file.
    """
    cluster_dict = defaultdict(list)
    
    with open(compacta_cluster_path, 'r') as cluster_file:
        for line in cluster_file:
            parts = line.strip().split('\t')
            if len(parts) > 1:
                transcript_id = parts[0]
                cluster_id = parts[1]
                cluster_dict[cluster_id].append(transcript_id)
    
    with open(output_file_path, 'w') as out_file:
        for cluster_id, transcript_ids in cluster_dict.items():
            for transcript_id in transcript_ids:
                out_file.write(f'{cluster_id}\t{transcript_id}\n')

###### FUNCTION TO GET MAPPING OF FASTA HEADERS TO SEQUENCE LENGTHS ######

def get_fasta_lengths_and_seqs(
    fasta_path
):
    """
    Read a FASTA file and return two dictionaries:
    1. A dictionary mapping headers to sequence lengths.
    2. A dictionary mapping headers to sequences.
    """
    header_to_length = {}
    header_to_sequence = {}
    
    for record in SeqIO.parse(fasta_path, 'fasta'):
        header = record.id
        sequence = str(record.seq)
        length = len(sequence)
        
        header_to_length[header] = length
        header_to_sequence[header] = sequence
    
    return header_to_length, header_to_sequence

###### FUNCTION TO TRANSLATE A DICT OF HEADERS TO CDS TO HEADERS TO PROT ######

def translate_cds(
    cds_dict
):
    """
    Translate nucleotide CDS sequences to protein sequences.
    """
    proteins_cds = {}

    for header, sequence in cds_dict.items():
        # Translate the nucleotide sequence to a protein sequence
        protein_seq = Seq(sequence).translate(to_stop=True)
        proteins_cds[header] = str(protein_seq)

    return proteins_cds

###### FUNCTION(S) TO COMBINE ALL CLUSTERING INFO INTO A SINGLE NETWORKX GRAPH ######

def read_clustering_info(
    file_path
):
    """
    Read a clustering info file and return a dictionary mapping cluster IDs to lists of transcript IDs.
    """
    cluster_dict = defaultdict(list)
    
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                cluster_id, transcript_id = parts
                cluster_dict[cluster_id].append(transcript_id)
    
    return cluster_dict

def create_graph_from_clusterings(
    trinity_path, cdhit_path, compacta_path
):
    """
    Create a graph from clustering info files and return the graph.
    """
    # Read clustering info files
    trinity_clusters = read_clustering_info(trinity_path)
    cdhit_clusters = read_clustering_info(cdhit_path)
    compacta_clusters = read_clustering_info(compacta_path)
    
    # Create the graph
    G = nx.Graph()
    
    # Function to add edges for clusters
    def add_edges_from_clusters(clusters):
        for transcript_ids in clusters.values():
            G.add_nodes_from(transcript_ids)
            for i in range(len(transcript_ids)):
                for j in range(i + 1, len(transcript_ids)):
                    G.add_edge(transcript_ids[i], transcript_ids[j])
    
    # Add nodes and edges for each clustering info file
    for clustering_name, clusters in [
        ('Trinity', trinity_clusters),
        ('CD-HIT', cdhit_clusters),
        ('Compacta', compacta_clusters)
    ]:
        add_edges_from_clusters(clusters)
        # Print the status to stderr
        sys.stderr.write(f'\tProcessing {clustering_name} clustering: ')
        sys.stderr.write(
            f'{G.number_of_nodes()} transcripts in '
            f'{nx.number_connected_components(G)} clusters.\n'
        )

    return G

###### FUNCTION TO SELECT LONGEST TRANSCRIPTS (FIRST BY ORF LENGTH, THEN TRANSCRIPT LENGTH) PER CLUSTER ######

def find_longest_transcripts(
    graph, orf_lengths, transcript_lengths
):
    """
    Find longest transcript in each connected component of the graph.
    """
    clusters_to_transcripts = {}
    clusters_to_longest_transcript = {}
    
    for i, component in enumerate(nx.connected_components(graph)):
        cluster_id = f'Unigene_{i+1}'
        transcript_list = list(component)
        clusters_to_transcripts[cluster_id] = transcript_list
        
        # Find the longest transcript by ORF length
        longest_orf_transcripts = []
        max_orf_length = -1
        for transcript in transcript_list:
            orf_length = orf_lengths.get(transcript, 0)
            if orf_length > max_orf_length:
                longest_orf_transcripts = [transcript]
                max_orf_length = orf_length
            elif orf_length == max_orf_length:
                longest_orf_transcripts.append(transcript)
        
        # If there's a tie, find the longest transcript by transcript length
        if len(longest_orf_transcripts) > 1:
            longest_transcript = max(
                longest_orf_transcripts, key=lambda t: transcript_lengths.get(t, 0)
            )
        else:
            longest_transcript = longest_orf_transcripts[0]
        
        clusters_to_longest_transcript[cluster_id] = longest_transcript
    
    return clusters_to_transcripts, clusters_to_longest_transcript
