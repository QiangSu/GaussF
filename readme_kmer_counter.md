# K-mer Counter

## Overview

This repository contains a C++ program designed for bioinformatics analysis. The program counts the occurrences of specific substrings known as k-mers within DNA sequences. These k-mers are provided in a comma-separated values (CSV) file, and the DNA sequences are read from FASTQ files. The program is designed to be efficient and can process large datasets by utilizing multiple threads.

## Features

- Count occurrences of k-mers from a list specified in a CSV file.
- Process sequencing data in FASTQ format, supporting gzip compression.
- Parallel processing using threads for improved performance on large datasets.
- Preservation of k-mer order between input and output files.

## How It Works

The program reads the k-mers from the first column of a CSV file, skipping the header line. It then processes a FASTQ file where the DNA sequences are read and counted for matches against the list of k-mers. The results are output to a new CSV file which includes the k-mer and its count, keeping the original order from the input file.

## Usage
To compile the program, ensure you have `g++`, `zlib`, and `pthreads` installed on your system, and then run:

```bash
g++ -std=c++17 -o kmer_counter kmer_counter.cpp -lz -lpthread

Once compiled, you can run the program with the following command:

./kmer_counter --kmer_dir <kmer_csv_dir> --fastq_file <fastq_gz_file> --output_dir <output_csv_dir> --threads <number_of_threads>


