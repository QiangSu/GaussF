Parallel K-mer Counting Tool
This Python script, titled parallel_kmer_counting.py, is designed to efficiently count the occurrences of specified k-mers within a gzipped FASTQ file, utilizing parallel processing to speed up the analysis.

Overview
K-mer analysis is a fundamental process in genomics and bioinformatics, which has applications in sequence assembly, alignment, and more. This tool takes a list of k-mers provided in CSV files and computes their counts within a given FASTQ dataset, with multiprocessing used to handle large datasets more effectively. It partitions the FASTQ file into chunks, allowing multiple processor cores to work on different sections of the data simultaneously.

Features
Parallel processing of FASTQ data to quickly count k-mers.
Accommodation of gzipped FASTQ files to save disk space.
Customizable chunk size to optimize memory usage during processing.
Utilization of multiprocessing based on the number of available CPU cores or a user-defined number of threads.
Output of k-mer count data into CSV format for further analysis.

Usage
The script is invoked from the command line and requires specific arguments pertaining to the FASTQ file and k-mer details:

python parallel_kmer_counting.py --k <k-mer size> --chunk_size <chunk size> --fastq <path to gzipped fastq file> --kmer_dir <directory with input CSVs> --output <output directory> --threads <number of threads>

--k: The size of the k-mer (required).
--chunk_size: The number of FASTQ records to process at once (required).
--fastq: The path to the gzipped FASTQ file for analysis (required).
--kmer_dir: The directory containing input CSV files, one k-mer per line, with the k-mer in the first column (required).
--output: The directory where the output CSV files will be saved (required).
--threads: The number of threads to use for processing (optional, defaults to the number of CPUs).

Output
The output is a collection of CSV files stored in the specified output directory, named according to the original k-mer CSV files with "_kmer_counts.csv" appended. Each output CSV file contains two columns: K-mer and Count. Each row represents a unique k-mer from the input file and its corresponding count within the FASTQ data.

Additional Notes
The script handles gzipped FASTQ files and expects .fastq.gz as input.
The performance gain from parallel processing depends on the computational resources available and the size of the dataset.
The chunk size should be chosen based on the available system memory, as larger chunk sizes consume more RAM.
Ensure that all input CSV files are correctly formatted with k-mers in the first column and no header.
