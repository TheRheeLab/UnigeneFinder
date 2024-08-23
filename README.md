# UnigeneFinder

## Overview

**UnigeneFinder** is an automated pipeline designed to simplify and improve gene prediction from de novo transcriptome assemblies, especially for non-model organisms that lack a reference genome. By integrating multiple clustering methods, UnigeneFinder significantly reduces the redundancy that often inflates the number of transcripts in raw assemblies. The pipeline outputs a set of primary transcripts, coding sequences, and protein sequences analogous to those typically available for high-quality reference genomes. This enables downstream analyses such as differential gene expression, ortholog identification, and evolutionary studies to be conducted with greater accuracy. It also calculates gene expression for putative unigenes as well as all transcript isoforms.

The pipeline is implemented in Python and is fully automated, making it accessible to a wide range of users. It is designed to run efficiently on both high-performance computing (HPC) systems and personal computers. For ease of use, **UnigeneFinder** is distributed with all necessary dependencies encapsulated in a Singularity container, ensuring smooth installation and execution on both Linux and Windows systems.

## Installation and Usage

### Requirements

- Linux or Windows (via WSL)
- Singularity (for containerized execution)

### Installation

Detailed instructions for installation on both Linux and Windows systems can be found in the operating system-specific installation guides. The included Singularity container in the ZIP archive under the Releases section contains all dependencies, simplifying the setup process.

### Running UnigeneFinder

Once installed, running UnigeneFinder is straightforward. Example commands and configurations are provided in the Usage Guide, along with a sample dataset available for download from the Releases section. This dataset includes example data and the Singularity container for quick testing and demonstration of the pipeline's capabilities.

**Note:** The ZIP archive under releases includes both the Singularity container and example data. These are not included when cloning the repository or downloading the source code archives.


## Testing the Pipeline on Example Data

To test the pipeline using the provided example data, follow these steps:

First, extract the example data from the provided archive:

```bash
tar -xvf ExampleData.tar.gz
```

To test the core UnigeneFinder functionality on the example data, which represents a subset of Arabidopsis transcripts and a small number of paired-end reads from two samples used in the benchmarking, run the following command:

```bash
singularity run unigenefinder.sif unigenefinder.py ExampleData/config.txt
```

This will create a folder called `UnigeneFinder_Working_Directory` under `ExampleData` where all output is stored. The main outputs of the program are the cluster information `.txt` file and the three `.fasta` files.

To test the mapping functionality, run:

```bash
singularity run unigenefinder_singularity.sif unigenefinder_map.py --ref ExampleData/Arabidopsis.fasta --reads ExampleData/FASTQ --type paired
```

This will create a folder called `unigenefinder_mapping` which will include the BAM file output. These files are the same as the BAM files in the `ExampleData` folder.

To test the expression calculation functionality, run:

```bash
singularity run unigenefinder.sif unigenefinder_expression.py --ref ExampleData/Arabidopsis_thaliana_Subset_500.fasta --map ExampleData/Arabidopsis_thaliana_Subset_500.map --bamdir ExampleData/BAM --type paired
```

This will create a folder called `unigenefinder_expression` which will include four CSV files: counts and TPM values for both unigenes and all isoforms.
