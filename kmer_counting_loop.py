import argparse
import os
import tarfile
from collections import Counter
import gzip
import multiprocessing
import pandas as pd
from Bio import SeqIO
from io import StringIO

# Define a function to count k-mers in a sequence
def count_kmers(sequence, k):
    return [sequence[i:i + k] for i in range(len(sequence) - k + 1)]

# Define a function to process a chunk of FASTQ data and return k-mers
def process_chunk(start, end, file_path, k):
    chunk_kmers = []
    # Open the FASTQ file and extract sequences within the specified chunk
    with gzip.open(file_path, 'rt') as fastq_file:
        records = SeqIO.parse(fastq_file, "fastq")
        for i, record in enumerate(records):
            if i >= start and i < end:
                sequence = str(record.seq)
                chunk_kmers.extend(count_kmers(sequence, k))
    return chunk_kmers

# Update the normalization function to include read length and k-mer size in the multiplier calculation
def normalize_frequency(kmer_counts, unique_kmer_count, total_read_count, read_length, k):
    # Adjust multiplier computation to correct for read length and k-mer size
    multiplier = 1000 * 1000000 / (unique_kmer_count * total_read_count * (read_length - k))
    return {kmer: count * multiplier for kmer, count in kmer_counts.items()}

# Process the CSV string, read and count FASTQ, normalize frequencies, and save to a file
def process_csv_string(csv_string, fastq_file_path, k, chunk_size, num_cores, read_length, output_directory, output_filename):
    if not csv_string.strip():
        print(f"No data in CSV string for file {output_filename}. Skipping.")
        return

    # Attempt to parse the CSV data
    try:
        unique_kmers_df = pd.read_csv(StringIO(csv_string), header=None)
        unique_kmers = unique_kmers_df[0].tolist()
    except pd.errors.EmptyDataError:
        print(f"Empty CSV data encountered for file {output_filename}. Skipping.")
        return

    # Read the FASTQ file and count the k-mers
    num_records, all_kmers = read_and_count_fastq(fastq_file_path, chunk_size, num_cores, k)
    num_unique_kmers = len(unique_kmers)

    # Counter for k-mer frequencies and normalization
    transcript2_kmer_counts = Counter(all_kmers)
    kmer_frequencies_normalized = normalize_frequency(transcript2_kmer_counts, num_unique_kmers, num_records, read_length, k)

    # Write the normalized frequencies to the output file
    output_file_path = os.path.join(output_directory, output_filename)
    with open(output_file_path, "w") as output_file:
        output_file.write("K-mer,Normalized Frequency\n")
        for kmer in unique_kmers:
            frequency = kmer_frequencies_normalized.get(kmer, 0)
            output_file.write(f"{kmer},{frequency}\n")

    print(f"Normalized frequency data saved to {output_file_path}")

# Define a function to read and count k-mers from a FASTQ file using multiprocessing
def read_and_count_fastq(fastq_file_path, chunk_size, num_cores, k):
    # Count the total number of records in the FASTQ file
    with gzip.open(fastq_file_path, 'rt') as fastq_file:
        num_records = sum(1 for line in fastq_file) // 4

    # Create multiprocessing pool and process each chunk
    pool = multiprocessing.Pool(processes=num_cores)
    num_chunks = (num_records // chunk_size) + (1 if num_records % chunk_size else 0)
    chunks = [(i * chunk_size, min((i + 1) * chunk_size, num_records)) for i in range(num_chunks)]

    # Using starmap to process multiple chunks in parallel
    results = pool.starmap(process_chunk, [(start, end, fastq_file_path, k) for start, end in chunks])
    pool.close()
    pool.join()

    # Combine all k-mers from processed chunks
    all_kmers = [kmer for result in results for kmer in result]
    return num_records, all_kmers

# Main function to parse arguments and invoke processing
def main():
    parser = argparse.ArgumentParser(description='K-mer normalization from CSV files in a tar.gz archive.')
    parser.add_argument('--k', type=int, default=50, help='Size of the k-mer (default: 50).')
    parser.add_argument('--chunk_size', type=int, default=10000, help='Number of records processed per chunk (default: 10000).')
    parser.add_argument('--fastq', type=str, required=True, help='Path to the FASTQ file.')
    parser.add_argument('--tar_gz', type=str, required=True, help='Path to the tar.gz file containing CSV files.')
    parser.add_argument('--output', type=str, default='./unique_frequence_csv/', help='Path to the output directory (default: ./unique_frequence_csv/).')
    parser.add_argument('--read_length', type=int, required=True, help='The length of the sequenced reads.')
    args = parser.parse_args()

    # Create the output directory if it does not exist
    output_directory = args.output
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Open tar.gz archive and process each CSV to generate normalized k-mer frequencies
    with tarfile.open(args.tar_gz, "r:gz") as tar:
        members = tar.getmembers()
        for member in members:
            if member.name.endswith(".csv") and member.isfile():
                file = tar.extractfile(member)
                if file:
                    csv_string = file.read().decode('utf-8')
                    output_filename = os.path.splitext(member.name)[0] + "_kmer_normalized_frequency_counts.csv"
                    process_csv_string(csv_string, args.fastq, args.k, args.chunk_size, multiprocessing.cpu_count(), args.read_length, args.output, output_filename)

if __name__ == "__main__":
    main()
