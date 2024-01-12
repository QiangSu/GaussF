# GaussF
GaussF Transcript Quantification Pipeline

Overview:

The GaussF pipeline is designed to accurately quantify transcript abundance at the isoform level using RNA-seq data. The algorithm leverages unique k-mer signatures to overcome common RNA-seq biases, employing a parametric Gaussian model for sophisticated bias correction. This methodology enables precise RPKM (Reads Per Kilobase of transcript, per Million mapped reads) estimates, facilitating in-depth transcriptomic analyses.

Step 1: transcriptome K-mer Analysis

kmer_frequency_distribution_mini_shared.py (isoform_unique_sequence-loop.py )

This tool processes a FASTA file containing transcript sequences and outputs a set of CSV files that summarize the k-mer content for each transcript. Each CSV file contains a list of k-mers of specified length that are present in the transcript, along with their local and global frequency, and information on which transcripts each k-mer is present in.

Features
K-mer Counting: For a given transcriptome FASTA file, count all k-mers of a specified length (default is set to 50).
Minimal Shared K-mers Output: For each transcript, output the k-mers that have the minimum global frequency—the smallest number of transcripts in which the k-mer appears.
CSV Output Content: Generate a CSV file for each isoform with the following columns:

'kmer': The k-mer sequence.
'Local_Frequency': The number of times the k-mer appears in the specific isoform.
'Global_Frequency': The number of transcripts that contain the k-mer across the entire transcriptome.
'Present_in_Transcripts': A list of transcript identifiers that share the k-mer if its global frequency is more than 1. For unique k-mers, the identifier of the single transcript is given.

Usage
To use this tool, you need to have Python installed on your system. The script requires a FASTA file with the transcript sequences as input and a directory path where the CSV files will be saved as output.

Execute the script with the necessary arguments from the command line. For example:

python kmer_frequency_distribution_mini_shared.py --input path/to/your/ACTB_reference/mart_export_ACTB.txt --output path/to/output/directory/

Command-Line Arguments
--input: Path to the input FASTA file containing transcript sequences.
--output: Path to the output directory where CSV files for each transcript will be saved.

Output File Details
For each transcript in the input FASTA file, the script will create a corresponding CSV file in the output directory with a name derived from the transcript header, sanitized to be filesystem-friendly.

In the output CSV files for each transcript, only k-mers that have the smallest global frequency for that transcript are included. If multiple k-mers share the same smallest global frequency, then all such k-mers are included in the CSV file. The 'Present_in_Transcripts' field in the CSV may include multiple transcript names, indicating that those transcripts share the k-mer.

If the global frequency of a k-mer is 1, indicating that it is unique to a single transcript, then the 'Present_in_Transcripts' field will only contain the identifier of that specific transcript.

Step 2.1: k-mer Counting and Normalization

kmer_counting_loop.py

Introduction
This tool is designed for bioinformatics analysis to count k-mer frequencies in sequencing data stored in FASTQ format. It is particularly useful when dealing with large datasets as it leverages Python's multiprocessing capabilities for parallel processing, thus enhancing performance and reducing computation time.
Features
Count specified k-mer sizes in FASTQ files (compressed with gzip).
Use input CSV files containing k-mer sequences to filter and count only relevant k-mers.
Handle large datasets efficiently with chunk-based parallel processing.
Utilize multiple CPU cores for faster computation.
Generate output CSV files containing the count of each k-mer.

Usage
To use this tool, you need to provide several command-line arguments. Here is the syntax for running the script:

python kmer_counter.py --k <kmer_size> --chunk_size <chunk_size> --fastq <fastq_file_path> --kmer_dir <kmer_directory> --output <output_directory> [--threads <number_of_threads>]

Command-Line Arguments
--k: Size of the k-mer you wish to count (required).
--chunk_size: Number of records from the FASTQ file to be processed in each parallel chunk (required).
--fastq: Path to the compressed FASTQ file that contains the sequencing data (required).
--kmer_dir: Directory path containing input CSV files with k-mer sequences; each CSV must have k-mers listed in the first column (required).
--output: Path to the output directory where the k-mer count CSV files will be saved (required). The output directory should be same as the input directory (kmer_dir)!!!
--threads: Number of threads to use for processing; defaults to the number of CPU cores available on the system.
Example usage:
python kmer_counting_loop.py --k 50 --threads 30 --chunk_size 10000000 --fastq /path/to/data.fastq.gz --kmer_dir /path/to/directory --output /path/to/directory

Output
The script will output CSV files in the specified output directory, with each file named according to the original k-mer CSV file but appended with _kmer_counts. For example, if there is an input file named sample_kmers.csv, the output file will be sample_kmers_kmer_counts.csv. Each output file will contain two columns: K-mer and Count, where Count is the frequency of that k-mer in the FASTQ file.
Performance
This tool is designed to handle large FASTQ files efficiently. By using parallel processing, the script splits the FASTQ file into chunks and processes each chunk in a separate CPU core, speeding up the counting operation significantly. The time taken will be printed at the end of the execution for each k-mer count task.

Step 2.2: K-mer Counts Merging and Normalization

merge_mormalize_isoform_count_v1.py

This script is designed to further process the output of a previous k-mer counting script. Its purpose is to merge the k-mer count data into the original k-mer CSV files and to normalize these counts to account for differences in the total number of k-mers and read counts. This is a necessary step in many bioinformatics workflows, particularly those involving comparative genomics or quantitative assessment of sequence representation.

Features
Merges k-mer count data with the original k-mer list CSV files.
Normalizes k-mer frequencies using the total k-mer counts and read lengths.
Supports input from gzipped FASTQ files for read count determination.
Efficiently calculates normalization factors and processes large datasets.

Example usage:
This script accepts command-line arguments to specify the input and output directories, the FASTQ file path, the read length, and the k-mer size. Here's how to run the script:

python merge_normalize_isoform_count_v1.py --directory <input_directory as the output directory of last kmer_counter.py script> --output_directory <output_directory new directory> --fastq <path_to_FASTQ.GZ> --read_length 150 --k 50

Command-Line Arguments
--directory: The directory containing the *_kmers.csv and corresponding *_kmer_counts.csv files (required). This directory is same as the output directory from the last script (kmer_counting_loop.py).
--output_directory: The directory where the merged and normalized CSV files will be saved (required). The output directory should be to a new directory for further GaussF workflow.
--fastq: The path to the gzipped FASTQ file for which k-mer counts were computed (required).
--read_length: The length of the reads in the FASTQ sequences, necessary for normalization (default is 150).
--k: The length of the k-mers used during the counting process (default is 50).
Output
For each *_kmers.csv file in the input directory, the script will save a corresponding *_merged_normalized.csv file in the output directory. This file will contain the original k-mer data, the raw count, and an additional column with normalized k-mer counts.

Step 3: Gaussian CDF Fitting for GC Content and Abundance Estimation

pipeline_abundance_GaussF_esti_loop.py

Introduction
This Python script is designed to analyze GC content distribution in sequence data and estimate the sequence abundance by fitting a cumulative distribution function (CDF) of a Gaussian to the GC content profile. It serves as a post-processing tool following k-mer counting, allowing researchers to derive meaningful biological insights based on the GC composition and k-mer abundance patterns.

Features
Analyzes the GC content of sequences represented by k-mers.
Performs fitting of a Gaussian CDF to the sum of normalized k-mer counts grouped by GC content percentage.
Extracts gene and transcript information from the input CSV filenames.
Produces structured output for quick assessment of fit success and estimated parameters.
Offers flexibility through user-defined minimum thresholds for k-mer counts appropriate for fitting.

Example usage:

python pipeline_abundance_GaussF_esti_loop.py --threshold 5 --input /path/to/merge_data --output / path/to/merge_data/results_file.txt

Command-Line Arguments
--input: The path to the input folder containing the k-mer CSV files where each file should have a filename format including gene and transcript IDs (e.g., GENE_ENST00001234567_kmers.csv) (required).
--output: The full path and name of the output CSV file where the results will be saved (required).
--threshold: The minimum number of k-mers required for performing fitting; the default value is 10 if not specified.
Output
The script will output a CSV file containing the following columns:

Gene: The gene identifier extracted from the input filename.
Transcript ID: The transcript identifier extracted from the input filename.
Global Frequency: The global frequency of the k-mer.
Present in Transcripts: The transcripts in which the k-mer is present, aggregated with a dash (-) separator.
Sum or Abundance RPKM: The cumulative sum of normalized k-mer counts if fitting failed, or the abundance represented by the amplitude parameter of the Gaussian CDF if fitting succeeded.
Mean: The mean of the Gaussian distribution, representing the center of the GC content distribution (only included if fitting succeeded).
SD: The standard deviation of the Gaussian distribution, representing the spread of the GC content distribution (only included if fitting succeeded).
The results file will also have entries with messages indicating whether fitting failed or whether data did not meet criteria for fitting (e.g., below threshold, not enough data points).


The interpretation of other files:

The 'human_pc_gene.tsv' file contains an annotated list of human protein-coding genes, each associated with a unique sequence region at the isoform level. This dataset comprises a total of 20,818 distinct genes derived from Homo sapiens, systematically cataloged to facilitate research on genetic variation and isoform expression.

The 'hm_mrna_isoform_summary.txt' details a comprehensive inventory of isoforms, including 
Transcript ID: Extracted from the filename by removing the '_kmers.csv' part, representing the gene name and transcript ID.
Number of Kmers: The total count of unique kmers found in each corresponding CSV file.
Global_Frequency: The frequency of the kmers across all sequences as found in the first data row of each corresponding CSV file.
Present_in_Transcripts: A list of transcripts in which the kmers are present as found in the first data row of each corresponding CSV file.

The summary encompasses 175,572 isoforms subjected to unique region screening. 

We welcome contributions to this project and encourage users to submit issues or pull requests on our GitHub repository. Additionally, if you have any questions or comments about this script, please feel free to contact us via email.
