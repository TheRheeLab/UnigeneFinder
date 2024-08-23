import os
import sys
import subprocess
import shutil

###### FUNCTION TO CHECK IF REQUIRED PROGRAMS ARE IN SYSTEM PATH ######

def check_programs() -> None:
    """
    Checks if specified programs are in the system path and exits if any are not found.
    """
    for program in ['bowtie2-build', 'bowtie2', 'samtools']:
        if not shutil.which(program):
            print(
                f'Error: {program} is not in the system path',
                file=sys.stderr
            )
            sys.exit(1)

###### FUNCTION TO BUILD BOWTIE2 INDEX ######

def build_bowtie2_index(fasta: str, index_dir: str, threads: int) -> None:
    """
    Builds a Bowtie2 index for the given reference FASTA file.
    """
    print('Building Bowtie2 index...', file=sys.stderr)
    cmd = [
        'bowtie2-build', '--threads', str(threads), fasta,
        os.path.join(index_dir, os.path.basename(fasta).split('.')[0])
    ]
    print(f'Running command: {" ".join(cmd)}', file=sys.stderr)
    result = subprocess.run(
        cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )
    if result.returncode != 0:
        print('Error building Bowtie2 index', file=sys.stderr)
        sys.exit(1)
    print('Index building completed', file=sys.stderr)

###### FUNCTION TO MAP SINGLE-END READS USING BOWTIE2 ######

def map_single_end_reads_bowtie2(
    fastq_dir: str, index_dir: str, output_dir: str,
    log_dir: str, fasta_base: str, threads: int
) -> None:
    """
    Maps single-end reads using Bowtie2 and sorts the output using samtools.
    """
    fastq_files = [
        f for f in os.listdir(fastq_dir)
        if f.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz'))
    ]
    print(
        f'Found {len(fastq_files)} FASTQ files in the specified folder.',
        file=sys.stderr
    )

    for fastq_file in fastq_files:
        fastq_path = os.path.join(fastq_dir, fastq_file)
        bam_output = os.path.join(
            output_dir, f'{os.path.splitext(fastq_file)[0]}.bam'
        )
        log_output = os.path.join(
            log_dir, f'{os.path.splitext(fastq_file)[0]}.log'
        )
        
        print(f'Mapping {fastq_file}...', file=sys.stderr)
        cmd = (
            f'bowtie2 -q --phred33 --sensitive --dpad 0 --gbar 99999999 '
            f'--mp 1,1 --np 1 --score-min L,0,-0.1 -I 1 -X 1000 '
            f'--no-mixed --no-discordant --no-unal --threads {threads} '
            f'-k 200 -x {index_dir}/{fasta_base} -U {fastq_path} | '
            f'samtools sort -n -o {bam_output} -'
        )
        print(f'Running command: {cmd}', file=sys.stderr)
        with open(log_output, 'w') as log_file:
            result = subprocess.run(
                cmd, shell=True, stdout=subprocess.DEVNULL, stderr=log_file
            )
        if result.returncode != 0:
            print(f'Error mapping {fastq_file}', file=sys.stderr)
            sys.exit(1)
        print(' Done', file=sys.stderr)

    print('All actions have completed', file=sys.stderr)

###### FUNCTION TO MAP PAIRED-END READS USING BOWTIE2 ######

def map_paired_end_reads_bowtie2(
    fastq_dir: str, index_dir: str, output_dir: str,
    log_dir: str, fasta_base: str, threads: int
) -> None:
    """
    Maps paired-end reads using Bowtie2 and sorts the output using samtools.
    """
    fastq_files = [
        f for f in os.listdir(fastq_dir)
        if f.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz'))
        and '_R1' in f
    ]
    print(
        f'Found {len(fastq_files)} FASTQ pairs in the specified folder.',
        file=sys.stderr
    )

    for fastq_file in fastq_files:
        fastq_file_R1 = fastq_file
        fastq_file_R2 = fastq_file.replace('_R1', '_R2')
        if not os.path.exists(os.path.join(fastq_dir, fastq_file_R2)):
            print(
                f'Error: Paired file {fastq_file_R2} not found for {fastq_file_R1}',
                file=sys.stderr
            )
            sys.exit(1)

        fastq_path_R1 = os.path.join(fastq_dir, fastq_file_R1)
        fastq_path_R2 = os.path.join(fastq_dir, fastq_file_R2)
        base_name = os.path.splitext(fastq_file_R1)[0].replace('_R1', '')
        bam_output = os.path.join(output_dir, f'{base_name}.bam')
        log_output = os.path.join(log_dir, f'{base_name}.log')
        
        print(
            f'Mapping {fastq_file_R1} and {fastq_file_R2}...',
            file=sys.stderr
        )
        cmd = (
            f'bowtie2 -q --phred33 --sensitive --dpad 0 --gbar 99999999 '
            f'--mp 1,1 --np 1 --score-min L,0,-0.1 -I 1 -X 1000 '
            f'--no-mixed --no-discordant --no-unal --threads {threads} '
            f'-k 200 -x {index_dir}/{fasta_base} -1 {fastq_path_R1} '
            f'-2 {fastq_path_R2} | samtools sort -n -o {bam_output} -'
        )
        print(f'Running command: {cmd}', file=sys.stderr)
        with open(log_output, 'w') as log_file:
            result = subprocess.run(
                cmd, shell=True, stdout=subprocess.DEVNULL, stderr=log_file
            )
        if result.returncode != 0:
            print(
                f'Error mapping {fastq_file_R1} and {fastq_file_R2}',
                file=sys.stderr
            )
            sys.exit(1)
        print(' Done', file=sys.stderr)

