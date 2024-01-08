# GaussF
GaussF Transcript Quantification Pipeline

Overview:

The GaussF pipeline is designed to accurately quantify transcript abundance at the isoform level using RNA-seq data. The algorithm leverages unique k-mer signatures to overcome common RNA-seq biases, employing a parametric Gaussian model for sophisticated bias correction. This methodology enables precise RPKM (Reads Per Kilobase of transcript, per Million mapped reads) estimates, facilitating in-depth transcriptomic analyses.

Step 1: transcriptome K-mer Analysis

kmer_frequency_distribution_mini_shared.py 

This tool processes a FASTA file containing transcript sequences and outputs a set of CSV files that summarize the k-mer content for each transcript. Each CSV file contains a list of k-mers of specified length that are present in the transcript, along with their local and global frequency, and information on which transcripts each k-mer is present in.

Features
K-mer Counting: For a given transcriptome FASTA file, count all k-mers of a specified length (default is set to 50).
Minimal Shared K-mers Output: For each transcript, output the k-mers that have the minimum global frequency—the smallest number of transcripts in which the k-mer appears.
CSV Output Content: Generate a CSV file for each isoform with the following columns:

'kmer': The k-mer sequence.
'Local_Frequency': The number of times the k-mer appears in the specific isoform.
'Global_Frequency': The number of transcripts that contain the k-mer across the entire transcriptome.
'Present_in_Transcripts': A list of transcript identifiers that share the k-mer if its global frequency is more than 1. For unique k-mers, the identifier of the single transcript is given.

Usage
To use this tool, you need to have Python installed on your system. The script requires a FASTA file with the transcript sequences as input and a directory path where the CSV files will be saved as output.

Command-Line Arguments
--input: Path to the input FASTA file containing transcript sequences.
--output: Path to the output directory where CSV files for each transcript will be saved.
Running the Script
Execute the script with the necessary arguments from the command line. For example:
python kmer_frequency_distribution_mini_shared.py --input path/to/your/ACTB_reference/mart_export_ACTB.txt --output path/to/output/directory/

Output File Details
For each transcript in the input FASTA file, the script will create a corresponding CSV file in the output directory with a name derived from the transcript header, sanitized to be filesystem-friendly.

In the output CSV files for each transcript, only k-mers that have the smallest global frequency for that transcript are included. If multiple k-mers share the same smallest global frequency, then all such k-mers are included in the CSV file. The 'Present_in_Transcripts' field in the CSV may include multiple transcript names, indicating that those transcripts share the k-mer.

If the global frequency of a k-mer is 1, indicating that it is unique to a single transcript, then the 'Present_in_Transcripts' field will only contain the identifier of that specific transcript.

Step 2: k-mer Counting and Normalization

This step involves counting the occurrence of transcript-specific k-mers within quality-controlled fastq datasets, implemented through kmer_counting_loop.py or kmer_counting_individual.py. For each isoform, k-mers and their normalized counts are documented in separate CSV files. These counts are normalized according to base-calling read length, sequencing depth, isoform length, and k-mer length.

Example usage:

bash
python kmer_counting_loop.py \
--k 50 \
--chunk_size 1000000 \
--num_threads 10 \
--read_length 151 \
--fastq ./N1_cutadapt_trim_2P.fastq.gz \
--csv_dir /path/to//combined_isoform_50mer_csv_files/ \
--output /path/to/unique_frequence_csv/ \


It is crucial to ensure that the reference .tar.gz file containing unique k-mer CSVs is located within the same repository as the kmer_counting_loop.py script.

Step 3: GaussF Model Analysis

Finally, the unique-region-derived k-mer counts are analyzed by the pipeline_abundance_GaussF_esti_loop.py script, featuring a self-benchmarking algorithm. The output comprises unbiased RPKM values for each transcript. Results are systematically collated within a text file, ready for downstream applications.

Example usage:
bash
python pipeline_abundance_GaussF_esti_loop.py \
--threshold 10 \
--input /path/to/unique_frequence_csv \
--output /path/to/results_file.txt \

#threshold is defined as the minimum number of unique k-mers required for the Gaussian fitting process. It serves as a cutoff for considering the quality and quantity of k-mer data that are sufficient to proceed with the calculation of GC content distribution. This ensures that the fitting is performed on datasets with a representative sample of k-mers from a given transcript isoform's unique region.

The input directory should ideally reside in the same local environment established during k-mer counting in Step 2.

Example Usage of GaussF Pipeline

In the gaussf > example_GaussF/ folder, we have uploaded 14 example CSV files representing transcript k-mers and their sequencing counts. These files are prepared for use with the GaussF pipeline and provide a comprehensive dataset to process with the script for an unbiased total count estimation.

This pipeline description is structured to guide users through a sequence of operations for estimating transcript isoform abundances within a dataset. The outlined workflow is not only a representation of the analytical procedure but also serves as instructive documentation for replicability and application to diverse transcriptomic datasets. For further insights and updates, users are invited to consult the pipeline's repository on GitHub.

The interpretation of other files:

The 'human_pc_gene.tsv' file contains an annotated list of human protein-coding genes, each associated with a unique sequence region at the isoform level. This dataset comprises a total of 20,818 distinct genes derived from Homo sapiens, systematically cataloged to facilitate research on genetic variation and isoform expression.

The 'isoform_unique_50mer_summary_list.txt' details a comprehensive inventory of isoforms, including their names and corresponding 50-mer sequences associated with unique genomic regions. The summary encompasses 175,572 isoforms subjected to unique region screening. Out of these, 121,144 isoforms are characterized by unique sequences, while 54,428 isoforms lack distinct sequences, sharing regions with alternate isoforms.

The documentation provided in 'README_GaussF.md' offers a detailed guide for utilizing the GaussF Python script. The README delineates the functional scope of the script, prerequisites for execution, and step-by-step instructions to ensure accurate computational procedures for users pursuing studies in Gaussian filtering.

In 'README_kmer_counting.md', users will find thorough instructions for operating the k-mer counting Python script. This README document includes the script's purpose, system requirements, installation guidelines, and a comprehensive command-line interface (CLI) usage tutorial, facilitating precise k-mer analysis in genomic sequences.

The 'README_reference_unique_sequence.md' serves as an instructional manual for the associated Python script that identifies unique sequences within reference genomes. It provides a clear exposition of the script's objectives, dependencies, configuration steps, and detailed execution commands, designed to aid researchers in pinpointing and studying unique genomic regions.

We welcome contributions to this project and encourage users to submit issues or pull requests on our GitHub repository. Additionally, if you have any questions or comments about this script, please feel free to contact us via email.
