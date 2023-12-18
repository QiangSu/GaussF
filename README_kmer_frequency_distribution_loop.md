# Kmer Frequency Distribution Loop

## Overview
The `kmer_frequency_distribution_loop.py` script is designed to analyze transcript sequences at the isoform level to determine the frequency of each k-mer (where k is a fixed length). The script calculates both the local frequency of k-mers within individual transcripts (isoform-specific) and the global frequency across all analyzed transcripts, highlighting the complexity of isoform overlapping regions.

## Purpose
The local frequency reflects the occurrence of each k-mer within a specific isoform, while the global frequency provides insight into the commonality of each k-mer across all isoforms. The output of the script can be used to study transcript diversity, uncover patterns of sequence conservation, and examine potential isoform overlaps.

## Features
- Processes multiple FASTA files containing transcript sequences.
- Calculates the frequency of k-mers of a fixed length (default is 50) for each isoform.
- Determines global k-mer frequency across all isoforms from all processed FASTA files.
- Outputs the results for each isoform to individual CSV files.
- Maintains the original order of k-mers as they appear in the transcript sequences.

## Usage
To use the script, ensure that your FASTA files are named consistently, following the pattern `mart_export1.txt` through `mart_exportN.txt`, and are located in the specified directory. Update the script to point to this directory if necessary.

The `output_directory` variable within the script defines the path to the directory where the output CSV files will be saved. Change this to an appropriate location on your system.

Run the script with Python 3:

```bash
python kmer_frequency_distribution_loop.py


Output
After running the script, the specified output directory will contain CSV files for each transcript from each processed FASTA file. Each CSV file is named using the isoform name from the FASTA header and contains the following columns:

kmer: The specific k-mer sequence.
Local Frequency: The count of how many times the k-mer appears in the specific isoform.
Global Frequency: The count of how many times the k-mer appears across all analyzed isoforms.
Dependencies
The script requires Python 3 and does not rely on any third-party libraries or modules beyond the Python Standard Library.