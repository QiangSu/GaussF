GC Content Analysis and Gaussian CDF Fitting Tool
This Python script, titled gc_content_gaussian_fit.py, performs a GC content analysis on k-mer sequences extracted from CSV files. It uses statistical modeling to fit the cumulative distribution of GC content percentages with Gaussian functions, providing insights into the nucleotide composition and its distribution.

Overview
GC content, the percentage of nucleotides in a DNA or RNA sequence that are guanine (G) or cytosine (C), is an important metric in genomic studies. This script calculates the GC content of k-mers and performs a curve-fitting analysis to model the distribution of these values across all k-mers in a dataset.

The script supports two models for fitting: a single Gaussian cumulative distribution function (CDF) and a bi-Gaussian CDF, the latter used in case the single Gaussian fit fails or is inadequate due to the distribution's complexity.

Features
Calculation of GC content for a list of k-mers obtained from CSV files.
Curve fitting to a single Gaussian CDF for GC content distribution analysis.
Alternative fitting to a bi-Gaussian CDF if the initial fitting is unsuccessful or if negative parameters are derived.
Handling of input data from multiple CSV files within a given directory.
Configuration of the minimum k-mer count threshold required for curve fitting to proceed.
Output of analysis results to a consolidated CSV file.

Usage
Run the script from the command line by providing the necessary arguments: the path to the input directory containing k-mer CSV files, the path and filename for the output CSV file, and an optional k-mer threshold for fitting:

python gc_content_gaussian_fit.py --input /path/to/kmer_csv_folder --output /path/to/results.csv --threshold 10

--input: Path to the input directory containing the CSV files with k-mer sequences (required).
--output: Path and name of the output file where results will be saved (required).
--threshold: Minimum number of k-mers required for fitting to proceed. Defaults to 10 (optional).

Output
The script outputs a CSV file containing the results of the analysis, which includes:

Gene: The gene name extracted from the CSV filename.
Transcript ID: The transcript ID extracted from the CSV filename.
Global Frequency: The frequency of occurrence across all k-mers.
Present in Transcripts: Transcript information associated with each k-mer.
Sum or Abundance RPKM: The sum of normalized k-mer counts or the calculated abundance in RPKM.
Mean: The mean value of the fitted Gaussian CDFs.
SD: The standard deviation of the fitted Gaussian CDFs.

Additional Notes
The script is configured to work with CSV files named in a specific pattern: GeneName_TranscriptID_*.csv. Ensure that your input files match this pattern for correct parsing of gene names and transcript IDs.
Before running the script, validate that the input CSV files are properly formatted with k-mers and the necessary data columns.
The script expects to fit the distribution of summed GC content percentages and may report issues if the input data is not supportive of a Gaussian distribution.
