import argparse
import os
import tarfile
from collections import Counter
import gzip
import multiprocessing
import pandas as pd
from Bio import SeqIO
from io import StringIO
import glob
import sys


# Define a function to count k-mers in a sequence
def count_kmers(sequence, k):
    return [sequence[i:i + k] for i in range(len(sequence) - k + 1)]


# Define a function to process a chunk of FASTQ data and return k-mers
def process_chunk(start, end, file_path, k):
    chunk_kmers = []
    with gzip.open(file_path, 'rt') as fastq_file:
        records = SeqIO.parse(fastq_file, "fastq")
        for i, record in enumerate(records):
            if i >= start and i < end:
                sequence = str(record.seq)
                chunk_kmers.extend(count_kmers(sequence, k))
    return chunk_kmers


# Update the normalization function to include read length and k-mer size in the multiplier calculation
def normalize_frequency(kmer_counts, unique_kmer_count, total_read_count, read_length, k):
    multiplier = 1000 * 1000000 / (unique_kmer_count * total_read_count * (read_length - k + 1))
    return {kmer: count * multiplier for kmer, count in kmer_counts.items()}


# Process the CSV string, read and count FASTQ, normalize frequencies, and save to a file
def process_csv_string(csv_string, fastq_file_path, k, chunk_size, num_cores, read_length, output_directory,
                       output_filename):
    if not csv_string.strip():
        print(f"No data in CSV string for file {output_filename}. Skipping.")
        return

    try:
        unique_kmers_df = pd.read_csv(StringIO(csv_string), header=None)
        unique_kmers = unique_kmers_df[0].tolist()
    except pd.errors.EmptyDataError:
        print(f"Empty CSV data encountered for file {output_filename}. Skipping.")
        return

    num_records, all_kmers = read_and_count_fastq(fastq_file_path, chunk_size, num_cores, k)
    num_unique_kmers = len(unique_kmers)

    kmer_counts = Counter(all_kmers)
    normalized_frequencies = normalize_frequency(kmer_counts, num_unique_kmers, num_records, read_length, k)

    output_file_path = os.path.join(output_directory, output_filename)
    with open(output_file_path, "w") as output_file:
        output_file.write("K-mer,Normalized Frequency\n")
        for kmer in unique_kmers:
            frequency = normalized_frequencies.get(kmer, 0)
            output_file.write(f"{kmer},{frequency}\n")

    print(f"Normalized frequency data saved to {output_file_path}")


def read_and_count_fastq(fastq_file_path, chunk_size, num_cores, k):
    with gzip.open(fastq_file_path, 'rt') as fastq_file:
        num_records = sum(1 for line in fastq_file) // 4

    pool = multiprocessing.Pool(processes=num_cores)
    num_chunks = (num_records // chunk_size) + (1 if num_records % chunk_size else 0)
    chunks = [(i * chunk_size, min((i + 1) * chunk_size, num_records)) for i in range(num_chunks)]

    results = pool.starmap(process_chunk, [(start, end, fastq_file_path, k) for start, end in chunks])
    pool.close()
    pool.join()

    all_kmers = [kmer for result in results for kmer in result]
    return num_records, all_kmers


def process_csv_directory(csv_dir, fastq_file_path, k, chunk_size, num_cores, read_length, output_directory):
    csv_files = glob.glob(os.path.join(csv_dir, '*.csv'))
    for csv_file_path in csv_files:
        with open(csv_file_path, 'r') as file:
            csv_string = file.read()
            output_filename = os.path.splitext(os.path.basename(csv_file_path))[
                                  0] + "_kmer_normalized_frequency_counts.csv"
            process_csv_string(csv_string, fastq_file_path, k, chunk_size, num_cores, read_length, output_directory,
                               output_filename)


def main():
    parser = argparse.ArgumentParser(
        description='K-mer normalization from CSV files in a directory or a tar.gz archive.')
    parser.add_argument('--k', type=int, default=50, help='Size of the k-mer (default: 50).')
    parser.add_argument('--chunk_size', type=int, default=10000,
                        help='Number of records processed per chunk (default: 10000).')
    parser.add_argument('--fastq', type=str, required=True, help='Path to the FASTQ file.')
    parser.add_argument('--tar_gz', type=str, help='Path to the tar.gz file containing CSV files.')
    parser.add_argument('--csv_dir', type=str, help='Path to the directory containing CSV files.')
    parser.add_argument('--output', type=str, default='./unique_frequence_csv/',
                        help='Path to the output directory (default: ./unique_frequence_csv/).')
    parser.add_argument('--read_length', type=int, required=True, help='The length of the sequenced reads.')
    parser.add_argument('--num_threads', type=int, default=multiprocessing.cpu_count(),
                        help='Number of threads for parallel processing (default: number of CPUs).')
    args = parser.parse_args()

    # Ensure either tar_gz or csv_dir is provided
    if not (args.tar_gz or args.csv_dir):
        print("Error: Either --tar_gz or --csv_dir must be specified.")
        parser.print_help()
        sys.exit(1)

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    if args.tar_gz:
        with tarfile.open(args.tar_gz, "r:gz") as tar:
            members = tar.getmembers()
            for member in members:
                if member.name.endswith(".csv") and member.isfile():
                    file = tar.extractfile(member)
                    if file:
                        csv_string = file.read().decode('utf-8')
                        output_filename = os.path.splitext(member.name)[0] + "_kmer_normalized_frequency_counts.csv"
                        process_csv_string(csv_string, args.fastq, args.k, args.chunk_size, args.num_threads,
                                           args.read_length, args.output, output_filename)
    elif args.csv_dir:
        process_csv_directory(args.csv_dir, args.fastq, args.k, args.chunk_size, args.num_threads, args.read_length,
                              args.output)


if __name__ == "__main__":
    main()
