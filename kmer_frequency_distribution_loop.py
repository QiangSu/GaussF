from collections import defaultdict
import csv
import os


def read_fasta(file_path):
    headers = []
    sequences = []
    sequence = ""
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence:
                    sequences.append(sequence)
                sequence = ""
                headers.append(line[1:])  # Remove the '>' character
            else:
                sequence += line
        if sequence:
            sequences.append(sequence)
    return headers, sequences


def sanitize_filename(header):
    return header.replace('|', '_')


# Create a single output directory
output_directory = './kmer_frequency_distribution'
os.makedirs(output_directory, exist_ok=True)

# Loop through multiple FASTA files
for i in range(1, 39):
    fasta_file = f'./mart_export{i}.txt'
    transcript_headers, transcripts = read_fasta(fasta_file)

    kmer_length = 50
    global_kmer_counts = defaultdict(int)

    # Process each FASTA file
    for isoform_index, sequence in enumerate(transcripts):
        # Output CSV file for each transcript
        sanitized_header = sanitize_filename(transcript_headers[isoform_index])
        output_csv_path = os.path.join(output_directory, f"mart_export{i}_{sanitized_header}_kmers.csv")

        with open(output_csv_path, mode='w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(['kmer', 'Local Frequency', 'Global Frequency'])

            # Initialize local kmer count dictionary
            local_kmer_counts = defaultdict(int)

            # Iterate over kmers in the original order they appear in the transcript
            for j in range(len(sequence) - kmer_length + 1):
                kmer = sequence[j:j + kmer_length]
                local_kmer_counts[kmer] += 1
                global_kmer_counts[kmer] += 1

            # Write kmers and their frequencies to the CSV file
            for j in range(len(sequence) - kmer_length + 1):
                kmer = sequence[j:j + kmer_length]
                local_freq = local_kmer_counts[kmer]
                global_freq = global_kmer_counts[kmer]
                csv_writer.writerow([kmer, local_freq, global_freq])

    print(f"Kmers for mart_export{i} have been saved as CSV files in the directory: {output_directory}")
