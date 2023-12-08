from Bio import SeqIO
import csv
import os

def count_kmers(sequence, k):
    kmers = []
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmers.append(kmer)
    return kmers

def gc_content(kmer):
    gc_count = kmer.count('G') + kmer.count('C')
    return (gc_count / len(kmer)) * 100

# Replace with the path to your FASTA file
fasta_file_path = "ACTB_NM_001101_1812nt.fasta"
k = 50

# Dictionary to keep track of GC content frequencies
gc_content_freq = {}

# Parse the FASTA file and calculate GC content frequencies
for record in SeqIO.parse(fasta_file_path, "fasta"):
    kmers = count_kmers(str(record.seq), k)
    for kmer in kmers:
        gc = round(gc_content(kmer))  # round to the nearest whole number
        gc_content_freq[gc] = gc_content_freq.get(gc, 0) + 1

# Sort the frequencies by GC content
sorted_gc_content_freq = sorted(gc_content_freq.items())

# Construct the output CSV file name to include the k-mer length and the name of the FASTA file
fasta_filename = os.path.basename(fasta_file_path)
fasta_filename_without_extension = os.path.splitext(fasta_filename)[0]
csv_file_path = f"{fasta_filename_without_extension}_k{k}_gc_content_frequencies.csv"

# Output the frequencies to a CSV file
with open(csv_file_path, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(["GC_Content", "Frequency"])  # Write header
    for gc, freq in sorted_gc_content_freq:
        csvwriter.writerow([gc, freq])  # Write GC content (without '%') and frequency

print(f"GC content frequencies have been written to {csv_file_path}")
