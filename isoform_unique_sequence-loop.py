from Bio import SeqIO
import csv
import os


# Functions from your original code:
def sliding_window(seq, kmer_size):
    return [seq[i:i + kmer_size] for i in range(len(seq) - kmer_size + 1)]


def unique_kmers(records, kmer_size):
    kmers_dict = {}
    for i, record in enumerate(records):
        for kmer in sliding_window(str(record.seq), kmer_size):
            if kmer not in kmers_dict:
                kmers_dict[kmer] = {i}
            else:
                kmers_dict[kmer].add(i)
    return kmers_dict


# K-mer size
kmer_size = 50

# Directory containing the input files
input_directory = '/home/boot/qiangsu/90-sample/kmer_isoform_gene'  # Adjust this to the location of your files

# Create output directory if it doesn't already exist
output_directory = 'kmer_isoform_gene'
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Loop through each file index to process files named 'mart_export.txt' to 'mart_export35.txt'
for file_index in range(1, 39):
    # Generate the file name based on the current index
    input_filename = f'mart_export{file_index}.txt'
    # Combine the input directory with the file name
    input_filepath = os.path.join(input_directory, input_filename)

    # Check if the file exists before processing
    if not os.path.exists(input_filepath):
        print(f"File not found: {input_filepath}")
        continue

    # Parse the FASTA file
    records = list(SeqIO.parse(input_filepath, 'fasta'))
    # Create a dictionary with all unique k-mers for all records in the file
    kmers_dict = unique_kmers(records, kmer_size)

    # Process each record within the current FASTA file
    for i, record in enumerate(records):
        # Extract list of unique k-mers for the current record
        unique_kmers_list = [kmer for kmer, record_numbers in kmers_dict.items() if record_numbers == {i}]

        # Replace "|" in the record.id with "_"
        record_id_for_filename = record.id.replace('|', '_')
        # Generate the output CSV file name
        output_filename = os.path.join(output_directory, record_id_for_filename + '_unique_kmers.csv')

        # Write unique k-mers to a CSV file
        with open(output_filename, 'w', newline='') as f:
            writer = csv.writer(f)
            for kmer in unique_kmers_list:
                writer.writerow([kmer])

    print(f"Processed {input_filename}")

# End of the script
