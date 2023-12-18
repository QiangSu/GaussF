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
                sequence = ""  # Reset the sequence
                headers.append(line[1:])  # Append the header without the '>'
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

# Create a dictionary to map kmers to the transcripts they are found in
kmer_transcript_mapping = defaultdict(set)

# Loop through multiple FASTA files
for i in range(1, 39):
    fasta_file = f'./mart_export{i}.txt'
    transcript_headers, transcripts = read_fasta(fasta_file)

    kmer_length = 50

    # Initialize a global kmer count dictionary
    global_kmer_counts = defaultdict(int)

    # Process each transcript in the FASTA file
    for isoform_index, sequence in enumerate(transcripts):
        # Initialize local kmer count dictionary
        local_kmer_counts = defaultdict(int)

        # Iterate over kmers in the transcript
        for j in range(len(sequence) - kmer_length + 1):
            kmer = sequence[j:j + kmer_length]
            local_kmer_counts[kmer] += 1
            global_kmer_counts[kmer] += 1
            kmer_transcript_mapping[kmer].add(transcript_headers[isoform_index])  # Add header to kmer mapping

    # After processing all transcripts, output to CSV for each transcript
    for isoform_index, header in enumerate(transcript_headers):
        output_csv_path = os.path.join(output_directory, f"mart_export{i}_{sanitize_filename(header)}_kmers.csv")

        with open(output_csv_path, mode='w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(['kmer', 'Local Frequency', 'Global Frequency', 'Found In Transcripts'])

            # Iterate through all kmers and write them to the CSV
            for kmer, local_freq in local_kmer_counts.items():
                global_freq = global_kmer_counts[kmer]
                found_in_transcripts = '; '.join(kmer_transcript_mapping[kmer] - {header})  # Exclude current transcript
                csv_writer.writerow([kmer, local_freq, global_freq, found_in_transcripts])

    # Notify the user for each file processed
    print(f"Kmers for mart_export{i} have been saved as CSV files in the directory: {output_directory}")
