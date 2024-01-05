# GaussF
GaussF Transcript Quantification Pipeline

Overview:

The GaussF pipeline is designed to accurately quantify transcript abundance at the isoform level using RNA-seq data. The algorithm leverages unique k-mer signatures to overcome common RNA-seq biases, employing a parametric Gaussian model for sophisticated bias correction. This methodology enables precise RPKM (Reads Per Kilobase of transcript, per Million mapped reads) estimates, facilitating in-depth transcriptomic analyses.

Step 1: Unique k-mer Identification

The script isoform_unique_sequence-loop.py enumerates distinct k-mers for each transcript. The targeted transcript set is defined by the research focus and can be customized based on specific gene lists. Approximately 175,000 transcripts at the isoform level can be processed to assign unique genomic regions. The targeted transcript set can be customized based on specific gene lists, making it suitable for a wide range of research applications. To use this script, simply run it with the desired input file containing transcript sequences and specify the k-mer length. The output will be a list of unique sequences for each transcript, along with their corresponding genomic coordinates. If you require the unique reference file, please send an email to su@chemie.uni-siegen.de and we will provide it to you promptly.
The kmer_frequency_distribution_loop.py script is designed to analyze transcript sequences at the isoform level to determine the frequency of each k-mer (where k is a fixed length). The script calculates both the local frequency of k-mers within individual transcripts (isoform-specific) and the global frequency across all analyzed transcripts, and track of the sets of transcripts that each kmer appears in, highlighting the complexity of isoform overlapping regions. The local frequency reflects the occurrence of each k-mer within a specific isoform, while the global frequency provides insight into the commonality of each k-mer across all isoforms. The output of the script can be used to study transcript diversity, uncover patterns of sequence conservation, and examine potential isoform overlaps.

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
