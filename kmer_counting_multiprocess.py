import time
from collections import Counter
import gzip
import multiprocessing
import pandas as pd
from Bio import SeqIO

def count_kmers(sequence, k):
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

def process_chunk(start, end, file_path, k):
    chunk_kmers = []
    with gzip.open(file_path, 'rt') as fastq_file:
        records = SeqIO.parse(fastq_file, "fastq")
        for i, record in enumerate(records):
            if i >= start:
                if i >= end:
                    break
                sequence = str(record.seq)
                chunk_kmers.extend(count_kmers(sequence, k))
    return chunk_kmers

def normalize_frequency(kmer_counts, unique_kmer_count, total_read_count, read_length, k):
    multiplier = 1000 * 1000000 / (unique_kmer_count * total_read_count * (read_length - k))
    return {kmer: count * multiplier for kmer, count in kmer_counts.items()}

read_length = 151  # Add this line with your correct read length
k = 50
chunk_size = 10000
num_cores = multiprocessing.cpu_count()


# File paths
csv_file_path = "GAPDH_unique_kmers.csv"
fastq_file_path = "1032245-N_cutadapt_trim_2P_100000.fastq.gz"
output_file_path = "kmer_normalized_frequency_counts1.csv"

# Calculate the number of records in the fastq file
with gzip.open(fastq_file_path, 'rt') as fastq_file:
    num_records = sum(1 for line in fastq_file) // 4

# Read the unique k-mers from the CSV file
unique_kmers_df = pd.read_csv(csv_file_path, header=None)
unique_kmers = unique_kmers_df[0].tolist()
num_unique_kmers = len(unique_kmers)

# Create chunk indices for parallel processing
num_chunks = (num_records // chunk_size) + (1 if num_records % chunk_size else 0)
chunks = [(i * chunk_size, min((i + 1) * chunk_size, num_records)) for i in range(num_chunks)]

start_time = time.time()

# Process the fastq file in parallel
pool = multiprocessing.Pool(processes=num_cores)
results = [pool.apply_async(process_chunk, args=(start, end, fastq_file_path, k)) for start, end in chunks]

# Collect k-mers from all chunks and count their occurrences
all_kmers = []
for result in results:
    all_kmers.extend(result.get())
transcript2_kmer_counts = Counter(all_kmers)

# Normalize the frequency counts
kmer_frequencies_normalized = normalize_frequency(transcript2_kmer_counts, num_unique_kmers, num_records, read_length, k)

# Write the normalized frequencies to a CSV file
with open(output_file_path, "w") as output_file:
    output_file.write("K-mer,Normalized Frequency\n")
    for kmer in unique_kmers:
        frequency = kmer_frequencies_normalized.get(kmer, 0)
        output_file.write(f"{kmer},{frequency}\n")

end_time = time.time()
execution_time = end_time - start_time
print(f"Normalized frequency data saved to {output_file_path}")
print(f"Execution Time: {execution_time:.2f} seconds")

pool.close()
pool.join()
