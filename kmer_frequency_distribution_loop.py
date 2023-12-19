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

def process_fasta_write_csv(input_txt_path, output_directory):
    # Read the FASTA file
    transcript_headers, transcripts = read_fasta(input_txt_path)

    # Ensure output directory exists
    os.makedirs(output_directory, exist_ok=True)
    kmer_length = 50
    global_kmer_counts = defaultdict(int)
    kmer_transcript_sets = defaultdict(set)  # Stores which transcripts contain the kmer

    # Count kmers globally and per transcript
    for isoform_index, sequence in enumerate(transcripts):
        for i in range(len(sequence) - kmer_length + 1):
            kmer = sequence[i:i + kmer_length]
            global_kmer_counts[kmer] += 1
            kmer_transcript_sets[kmer].add(isoform_index)

    # Output kmer data to CSV files
    for isoform_index, header in enumerate(transcript_headers):
        output_csv_path = os.path.join(output_directory, sanitize_filename(header) + '_kmers.csv')

        with open(output_csv_path, mode='w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(['kmer', 'Local Frequency', 'Global Frequency', 'Present in Transcripts'])

            # Local kmer occurrences for this transcript
            local_kmer_counts = defaultdict(int)
            for i in range(len(transcripts[isoform_index]) - kmer_length + 1):
                kmer = transcripts[isoform_index][i:i + kmer_length]
                local_kmer_counts[kmer] += 1

            # Write the kmers, frequencies and sets of transcripts where the kmer is present
            for kmer, local_freq in local_kmer_counts.items():
                global_freq = global_kmer_counts[kmer]
                transcripts_containing_kmer = ', '.join(transcript_headers[i] for i in kmer_transcript_sets[kmer])
                csv_writer.writerow([kmer, local_freq, global_freq, transcripts_containing_kmer])

    print(f"Kmers for each transcript from {input_txt_path} have been saved as CSV files in the directory: {output_directory}")

# Main loop to process each FASTA file
input_directory = './'
output_directory = './kmer_frequency_distribution'  # The single output directory
os.makedirs(output_directory, exist_ok=True)  # Ensure the output directory exists

for i in range(1, 39):
    input_txt_path = os.path.join(input_directory, f"mart_export{i}.txt")
    process_fasta_write_csv(input_txt_path, output_directory)
