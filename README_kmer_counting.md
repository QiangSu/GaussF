# K-mer Counting Tool

## Description

This script processes sequences from a provided FASTQ file and counts the occurrences of k-mers. The counts are then normalized based on k-mer frequencies provided in CSV files within a `.tar.gz` archive. This is useful for genomic and bioinformatics analysis.

## Prerequisites

Before running this script, make sure to have the following installed:
- Python 3
- Pandas library
- BioPython library

## Installation

Clone the repository to your local machine using the following command:

```shell
git clone https://github.com/your-username/your-repository-name.git


Navigate to the cloned directory and install the required Python packages using:

shell
pip install -r requirements.txt



Usage
Run the script from the command line by specifying the mandatory arguments.

shell
python kmer_counting_loop.py --k <k-mer size> --chunk_size <chunk size> --fastq <path to FASTQ file> --tar_gz <path to tar.gz file> --output <output directory>

Replace the placeholders with your specific values:

<k-mer size> - The length of the k-mer (default is 50).
<chunk size> - The chunk size used for processing the FASTQ file (default is 10000).
<path to FASTQ file> - The path to your FASTQ file.
<path to tar.gz file> - The path to the tar.gz archive containing your CSV files with k-mer counts.
<output directory> - The directory where the output CSV files will be saved.


Example Command
shell
python kmer_counting_loop.py --k 50 --chunk_size 10000 --fastq "C:/path/to/1032245-T_cutadapt_trim_2P_100000.fastq.gz" --tar_gz "C:/path/to/combined_isoform_50mer_csv_files.tar.gz" --output "C:/path/to/output_directory"

Contributing
Contributions to this project are welcome! 

Support and Feedback
If you encounter any issues or have feedback, please file an issue on the GitHub repository issue tracker.

License
For commercial purpose should be discussed with author.

Acknowledgements
Special thanks to all contributors and to the bioinformatics community for their valuable insights.
















