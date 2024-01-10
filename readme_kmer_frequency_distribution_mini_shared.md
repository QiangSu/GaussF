This script processes a FASTA formatted file containing nucleotide sequences and exports kmer occurrence data for each sequence into individual CSV files.

Overview
A kmer is a substring of length k and kmer analysis is a common technique in bioinformatics for various applications including genome assembly, sequence alignment, and error correction. This script reads the nucleotide sequences from a FASTA file, counts the kmers within each sequence, and outputs CSV files containing the kmer analysis for each individual sequence.

Features
FASTA file parsing to extract headers and sequences
Calculation of kmer occurrences within each sequence with a set kmer length (default: 50)
Evaluation of kmer global frequency (across all sequences) and local frequency (within a single sequence)
Generation of CSV files for each sequence containing information about the kmers found within the sequence and their respective frequency statistics
Sanitization of headers for valid filenames during CSV file creation
Option to include only kmers with a minimum global occurrence (useful for filtering common kmers)
Command-line interface for specifying input and output locations
Directory creation for output files

Usage
The script is executed from the command line and requires two arguments: the path to the input FASTA file and the path to the output directory where the CSV files will be saved.

python kmer_analysis.py --input ./path/to/input.fasta --output ./path/to/output_directory

Output
The output consists of a series of CSV files saved in the specified output directory. Each file is named after the corresponding sequence header found in the FASTA file, sanitized to ensure valid filename format. The CSV files contain four columns:

kmer: the specific kmer sequence
Local_Frequency: the number of times the kmer appears in that specific sequence
Global_Frequency: the number of times the kmer appears across all sequences in the input file
Present_in_Transcripts: a list of sequence headers where the kmer is found (if the global frequency is greater than 1)
The script filters kmers to only include those with the minimum global frequency for each isoform while writing the CSV files, thus focusing on unique or characteristic kmers for each sequence.
