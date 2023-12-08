# pipeline_abundance_GaussF_esti_loop Tool

## Description

This Python script processes CSV files containing k-mer frequency data to calculate the GC content of each k-mer and then performs statistical analysis by fitting the accumulated frequency data to a Gaussian cumulative distribution function (CDF). It calculates the amplitude of the Gaussian CDF for each CSV file and outputs the results into a text file.

## Features

- Calculation of GC content for k-mers
- Grouping and summing frequencies by GC content
- Fitting data to a Gaussian CDF
- Outputting the amplitude of the Gaussian CDF fit for each CSV file

## Prerequisites

To run this script, you will need:
- Python 3
- Pandas (Python library)
- NumPy (Python library)
- SciPy (Python library)

## Installation

First, clone this repository to your local machine:

```shell
git clone https://github.com/your-username/pipeline_abundance_GaussF_esti_loop.git


shell
cd kmer-gc-content-analysis
pip install pandas numpy scipy


Usage
Place your CSV files inside a subfolder named unique_frequence_csv in the project directory. Run the script:

shell
python pipeline_abundance_GaussF_esti_loop.py --input ./unique_frequence_csv --output ./reslut1.txt
 

To run this Python script, users will execute a command like this:

bash
python pipeline_abundance_GaussF_esti_loop.py unique_frequence_csv

The script will create an output file named amplitude_fit_cdf_results.txt in the project directory containing the amplitude of the Gaussian CDF fit for each k-mer CSV file.

Output
The results are stored in a file named amplitude_fit_cdf_results.txt in the following format:

file1.csv: amplitude_value1
file2.csv: amplitude_value2
...