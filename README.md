# UnigeneFinder

## Overview

**UnigeneFinder** is an automated pipeline designed to simplify and improve gene prediction from de novo transcriptome assemblies, especially for non-model organisms that lack a reference genome. By integrating multiple clustering methods, UnigeneFinder significantly reduces the redundancy that often inflates the number of transcripts in raw assemblies. The pipeline outputs a set of primary transcript coding sequences and protein sequences analogous to those typically available for high-quality reference genomes. This enables downstream analyses such as differential gene expression, ortholog identification, and evolutionary studies to be conducted with greater accuracy.

The pipeline is implemented in Python and is fully automated, making it accessible to a wide range of users. It is designed to run efficiently on both high-performance computing (HPC) systems and personal computers. For ease of use, **UnigeneFinder** is distributed with all necessary dependencies encapsulated in a Singularity container, ensuring smooth installation and execution on both Linux and Windows systems.

## Installation and Usage

### Requirements

- Linux or Windows (via WSL)
- Singularity (for containerized execution)

### Installation

Detailed instructions for installation on both Linux and Windows systems can be found in the operating system-specific installation guides. The included Singularity container in the ZIP archive under the Releases section contains all dependencies, simplifying the setup process.

### Running UnigeneFinder

Once installed, running UnigeneFinder is straightforward. Example commands and configurations are provided in the Usage Guide, along with a sample dataset available for download from the Releases section. This dataset includes example data and the Singularity container for quick testing and demonstration of the pipeline's capabilities.

**Note:** The ZIP archive under releases includes both the Singularity container and example data, making it easy to get started with UnigeneFinder right away.
