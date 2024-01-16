import argparse
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import norm
import os
import re
import warnings
from scipy.optimize import OptimizeWarning

# Suppress the OptimizeWarning from SciPy
warnings.simplefilter("ignore", OptimizeWarning)

warnings.simplefilter("ignore", RuntimeWarning)

# Function to calculate the GC content percentage of a sequence
def calculate_gc_content(seq):
    gc_count = seq.count('G') + seq.count('C')
    gc_content_percentage = (gc_count / len(seq)) * 100
    return gc_content_percentage

# Single Gaussian CDF function
def gaussian_cdf(x, A0, A, xc, w):
    return A0 + A * norm.cdf((x - xc) / w)

# Function to extract gene and transcript information from filename, including gene names with hyphens
def extract_gene_transcript(file_name):
    match = re.match(r'([A-Za-z0-9-]+)_(ENST\d+)_.*\.csv$', file_name)
    if match:
        return match.groups()
    else:
        return None, None

# Function to create the output string for unfit data
def format_unfit_output(gene_name, transcript_id, global_frequency, present_transcript, sum_frequency, message):
    return f"{gene_name},{transcript_id},{global_frequency},{present_transcript},{sum_frequency:.2f},,{message}\n"

# Main procedure
def main(input_folder, output_file_path, row_threshold):
    input_folder_path = os.path.abspath(input_folder)
    csv_files = [file for file in os.listdir(input_folder_path) if file.endswith('.csv')]

    with open(output_file_path, 'w') as f_out:
        f_out.write('Gene,Transcript ID,Global Frequency,Present in Transcripts,Sum or Abundance RPKM,Mean,SD\n')

    for file_name in csv_files:
        df = pd.read_csv(os.path.join(input_folder_path, file_name))

        global_frequency = df.at[0, 'Global_Frequency']
        present_transcripts = "-".join(df['Present_in_Transcripts'].unique()) if global_frequency > 1 else df.at[0, 'Present_in_Transcripts']
        df['Normalized_K-mer_Count'] = pd.to_numeric(df['Normalized_K-mer_Count'], errors='coerce')
        df.dropna(subset=['Normalized_K-mer_Count'], inplace=True)
        df['GC Content (%)'] = df['kmer'].apply(calculate_gc_content)

        gc_content_grouped = df.groupby('GC Content (%)')['Normalized_K-mer_Count'].agg(['sum']).reset_index()
        gc_content_grouped.sort_values('GC Content (%)', inplace=True)

        sum_frequency = df['Normalized_K-mer_Count'].sum()
        gene_name, transcript_id = extract_gene_transcript(file_name)
        if not gene_name or not transcript_id:
            print(f"Warning: Could not parse gene name and transcript ID from {file_name}")
            continue

        if len(df['kmer']) < row_threshold:
            result_string = format_unfit_output(gene_name, transcript_id, global_frequency, present_transcripts, sum_frequency, "Below k-mer threshold: no fitting performed")
            with open(output_file_path, "a") as f_out:
                f_out.write(result_string)
            continue

        if len(gc_content_grouped) < 4:  # There are four parameters in the gaussian_cdf function
            result_string = format_unfit_output(gene_name, transcript_id, global_frequency, present_transcripts, sum_frequency, "Not enough GC content data points for fitting")
            with open(output_file_path, "a") as f_out:
                f_out.write(result_string)
            continue

        try:
            y_data = gc_content_grouped['sum'].cumsum()
            initial_guess = [0, y_data.iloc[-1], gc_content_grouped['GC Content (%)'].mean(),
                             gc_content_grouped['GC Content (%)'].std()]
            popt, pcov = curve_fit(gaussian_cdf, gc_content_grouped['GC Content (%)'], y_data, p0=initial_guess,
                                   maxfev=9000)

            # Check for a negative mean value and mean greater than standard deviation
            if popt[2] < 0 or popt[2] <= popt[3]:
                raise ValueError("Negative mean or mean not greater than standard deviation")

            # Check if the fitted amplitude is more than double the sum of individual values
            if popt[1] > 2 * sum_frequency:
                raise ValueError("Fitted amplitude is more than double the sum of individual values")

            # Check if the standard deviation (popt[3]) is negative
            if popt[3] < 0:
                raise ValueError("Standard deviation is negative")

            # Fitting succeeded, so output normal results
            result_string = f"{gene_name},{transcript_id},{global_frequency},{present_transcripts},{popt[1]:.2f},{popt[2]:.2f},{popt[3]:.2f}\n"
        except (RuntimeError, ValueError) as e:
            # Fitting failed or one of the checks did not pass; output the sum of individual values
            error_message = str(e)
            result_string = format_unfit_output(gene_name, transcript_id, global_frequency, present_transcripts,
                                                sum_frequency, f"Fitting failed: {error_message}")

        with open(output_file_path, "a") as f_out:
            f_out.write(result_string)

    print(f"Results have been written to {output_file_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze GC content and fit Gaussian CDF.")
    parser.add_argument('--input', type=str, required=True, help="Path to the input folder containing the CSV files.")
    parser.add_argument('--output', type=str, required=True, help="Path and name of the output file to save the results.")
    parser.add_argument('--threshold', type=int, default=10, help="Minimum number of k-mers required for fitting. Default is 10.")
    args = parser.parse_args()

    # Ensure the output directory exists
    output_directory = os.path.dirname(args.output)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    main(args.input, args.output, args.threshold)
