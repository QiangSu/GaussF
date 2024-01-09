import argparse
import time
import multiprocessing
import os
import gzip
from collections import OrderedDict
from Bio import SeqIO
import glob

def count_kmers(sequence, k, filter_set=None):
    kmers = [sequence[i:i + k] for i in range(len(sequence) - k + 1)]
    return kmers if filter_set is None else [kmer for kmer in kmers if kmer in filter_set]

def process_chunk(start, end, file_path, k, filter_set):
    chunk_kmers = []
    with gzip.open(file_path, 'rt') as fastq_file:
        records = SeqIO.parse(fastq_file, "fastq")
        for i, record in enumerate(records):
            if start <= i < end:
                sequence = str(record.seq)
                chunk_kmers.extend(count_kmers(sequence, k, filter_set))
    return chunk_kmers

def main():
    parser = argparse.ArgumentParser(description='Count k-mer frequencies in a FASTQ file.')
    parser.add_argument('--k', type=int, help='Size of the k-mer.', required=True)
    parser.add_argument('--chunk_size', type=int, help='Number of records processed per chunk.', required=True)
    parser.add_argument('--fastq', type=str, help='Path to the FASTQ file.', required=True)
    parser.add_argument('--kmer_dir', type=str, help='Directory containing input CSV files with k-mer sequences.', required=True)
    parser.add_argument('--output', type=str, help='Output directory for storing CSV files.', required=True)
    parser.add_argument('--threads', type=int, help='Number of threads to use for processing.', default=multiprocessing.cpu_count())
    args = parser.parse_args()

    # Get arguments
    k = args.k
    chunk_size = args.chunk_size
    fastq_file_path = args.fastq
    kmer_dir = args.kmer_dir
    output_directory = args.output
    num_cores = args.threads
    read_length = 150

    # Create the output directory if it doesn't exist
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Calculate the number of records in the FASTQ file
    with gzip.open(fastq_file_path, 'rt') as fastq_file:
        num_records = sum(1 for _ in SeqIO.parse(fastq_file, "fastq"))

    # Get a list of all CSV files in the directory
    kmer_files = glob.glob(os.path.join(kmer_dir, "*_kmers.csv"))

    for kmer_csv_file_path in kmer_files:
        # Load k-mers from the first column of the CSV file into an OrderedDict to preserve order
        kmers_from_csv = OrderedDict()
        with open(kmer_csv_file_path, 'r') as csvfile:
            for line in csvfile:
                kmer = line.split(',')[0].strip()
                if not kmer.lower().startswith("kmer") and kmer:  # Skip header or empty lines
                    kmers_from_csv[kmer] = 0  # Initialize the count as zero

        # Extract the base name of the CSV file to use in the output file name
        csv_base_name = os.path.splitext(os.path.basename(kmer_csv_file_path))[0]
        output_file_path = os.path.join(output_directory, f"{csv_base_name}_kmer_counts.csv")

        # Create chunk indices for parallel processing
        num_chunks = (num_records // chunk_size) + (1 if num_records % chunk_size != 0 else 0)
        chunks = [(i * chunk_size, min((i + 1) * chunk_size, num_records)) for i in range(num_chunks)]

        start_time = time.time()

        # Process the FASTQ file in parallel using multiprocessing
        with multiprocessing.Pool(processes=num_cores) as pool:
            results = [pool.apply_async(process_chunk, args=(start, end, fastq_file_path, k, set(kmers_from_csv))) for start, end in chunks]

            # Flatten the list of k-mers from all chunks and update the counts in the ordered dictionary
            for result in results:
                chunk_kmers = result.get()
                for kmer in chunk_kmers:
                    if kmer in kmers_from_csv:
                        kmers_from_csv[kmer] += 1

        # Write the k-mer counts to a CSV file in the specified output directory
        with open(output_file_path, "w") as output_file:
            output_file.write("K-mer,Count\n")
            for kmer, count in kmers_from_csv.items():
                output_file.write(f"{kmer},{count}\n")

        end_time = time.time()
        execution_time = end_time - start_time
        print(f"k-mer count data for {csv_base_name} saved to {output_file_path}")
        print(f"Execution Time: {execution_time:.2f} seconds")

if __name__ == "__main__":
    main()
