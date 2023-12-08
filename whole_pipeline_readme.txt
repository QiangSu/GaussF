GaussF Transcript Quantification Pipeline

Overview:

The GaussF pipeline is designed to accurately quantify transcript abundance at the isoform level using RNA-seq data. The algorithm leverages unique k-mer signatures to overcome common RNA-seq biases, employing a parametric Gaussian model for sophisticated bias correction. This methodology enables precise RPKM (Reads Per Kilobase of transcript, per Million mapped reads) estimates, facilitating in-depth transcriptomic analyses.

Step 1: Unique k-mer Identification

The script isoform_unique_sequence-loop.py enumerates distinct k-mers for each transcript. The targeted transcript set is defined by the research focus and can be customized based on specific gene lists. Approximately 175,000 transcripts at the isoform level can be processed to assign unique genomic regions.

Step 2: k-mer Counting and Normalization

This step involves counting the occurrence of transcript-specific k-mers within quality-controlled fastq datasets, implemented through kmer_counting_loop.py or kmer_counting_individual.py. For each isoform, k-mers and their normalized counts are documented in separate CSV files. These counts are normalized according to base-calling read length, sequencing depth, isoform length, and k-mer length.

Example usage:

bash
python kmer_counting_loop.py \
--k 50 \
--chunk_size 10000 \
--read_length 151 \
--fastq ./1032245-T_cutadapt_trim_2P_100000.fastq.gz \
--tar_gz /path/to/combined_isoform_50mer_csv_files.tar.gz \
--output /path/to/unique_frequence_csv/

It is crucial to ensure that the reference .tar.gz file containing unique k-mer CSVs is located within the same repository as the kmer_counting_loop.py script.

Step 3: GaussF Model Analysis

Finally, the unique-region-derived k-mer counts are analyzed by the pipeline_abundance_GaussF_esti_loop.py script, featuring a self-benchmarking algorithm. The output comprises unbiased RPKM values for each transcript. Results are systematically collated within a text file, ready for downstream applications.

Example usage:
bash
python pipeline_abundance_GaussF_esti_loop.py \
--input /path/to/unique_frequence_csv \
--output /path/to/results_file.txt

The input directory should ideally reside in the same local environment established during k-mer counting in Step 2.


This pipeline description is structured to guide users through a sequence of operations for estimating transcript isoform abundances within a dataset. The outlined workflow is not only a representation of the analytical procedure but also serves as instructive documentation for replicability and application to diverse transcriptomic datasets. For further insights and updates, users are invited to consult the pipeline's repository on GitHub.