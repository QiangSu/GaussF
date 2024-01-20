import os
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import norm
import argparse
import re

# Parse command line arguments
parser = argparse.ArgumentParser(description="Analyze GC content and fit Gaussian CDF.")
parser.add_argument('--input', type=str, required=True, help="Path to the input folder containing the CSV files.")
parser.add_argument('--output', type=str, required=True, help="Path and name of the output file to save the results.")
parser.add_argument('--threshold', type=int, default=10, help="Minimum number of k-mers required for fitting. Default is 10.")
args = parser.parse_args()

# GC content calculation function
def calculate_gc_content(kmer):
    gc_count = kmer.count('G') + kmer.count('C')
    total_bases = len(kmer)
    gc_content_percent = (gc_count / total_bases) * 100
    return round(gc_content_percent, 2)

# Gaussian CDF definition
def gaussian_cdf(x, A0, A, xc, w):
    return A0 + A * norm.cdf((x - xc) / w)

# Gaussian CDF with fixed parameters
def gaussian_cdf_fixed(x, A0, A, xc_fixed, w_fixed):
    return A0 + A * norm.cdf((x - xc_fixed) / w_fixed)

# Extract gene name and transcript ID from filename
def extract_gene_transcript_id(filename):
    match = re.search(r'(\w+)_(\w+)_kmers', filename)
    if match:
        return match.group(1), match.group(2)
    else:
        return "Unknown", "Unknown"

# List to store the results
results = []

# Loop through each file in the directory
for filename in os.listdir(args.input):
    if filename.endswith("merged_normalized.csv"):
        filepath = os.path.join(args.input, filename)

        # Extract gene name and transcript ID from the filename
        gene_name, transcript_id = extract_gene_transcript_id(filename)

        # Read the CSV file into a DataFrame
        df = pd.read_csv(filepath)

        # Ensure the 'kmer' column exists
        if 'kmer' not in df.columns:
            print(f"'kmer' column is missing in file {filename}. Skipping this file.")
            continue

        # Calculate and add GC content to the DataFrame
        df['GC_Content'] = df['kmer'].apply(calculate_gc_content)

        # Group by GC content and sum frequencies
        gc_content_data = df.groupby('GC_Content').agg({
            'Local_Frequency': 'sum',
            'Normalized_K-mer_Count': 'sum'
        }).reset_index()

        # Add the first row's Global_Frequency and Present_in_Transcripts to every result
        global_frequency = df.at[0, 'Global_Frequency']
        present_in_transcripts = df.at[0, 'Present_in_Transcripts']

        # Check if there are at least args.threshold distinct GC contents
        if len(gc_content_data['GC_Content'].unique()) < args.threshold:
            sum_normalized_kmer_count = gc_content_data['Normalized_K-mer_Count'].sum()
            results.append({
                'File': filename,
                'Gene_Name': gene_name,
                'Transcript_ID': transcript_id,
                'Global_Frequency': global_frequency,
                'Present_in_Transcripts': present_in_transcripts,
                'Sum or Fitted A (Abundance)': '{:.2f}'.format(sum_normalized_kmer_count),
                'Fixed Mean (xc)': 'N/A',
                'Fixed Standard Deviation (w)': 'N/A',
                'Report': 'Insufficient Data'
            })
            continue

        # Sort and calculate cumulative sums
        gc_content_data_sorted = gc_content_data.sort_values(by='GC_Content')
        gc_content_data_sorted['Cumulative_Local_Frequency'] = gc_content_data_sorted['Local_Frequency'].cumsum()
        gc_content_data_sorted['Cumulative_Normalized_Count'] = gc_content_data_sorted['Normalized_K-mer_Count'].cumsum()

        # Get the data for fitting
        x_data = gc_content_data_sorted['GC_Content']
        y_data_local = gc_content_data_sorted['Cumulative_Local_Frequency']
        y_data_normalized = gc_content_data_sorted['Cumulative_Normalized_Count']

        initial_guesses_local = [min(y_data_local), max(y_data_local) - min(y_data_local), x_data.mean(), x_data.std()]
        try:
            # Fit the Gaussian CDF to Cumulative Local Frequency
            popt_local, pcov_local = curve_fit(gaussian_cdf, x_data, y_data_local, p0=initial_guesses_local)
            A0_fitted_local, A_fitted_local, xc_fitted_local, w_fitted_local = popt_local

            # Check if the fixed mean (xc) is more than twice the fixed standard deviation (w)
            if xc_fitted_local > 1 * w_fitted_local:
                # Proceed with normalized count fitting using fixed xc and w
                initial_guesses_normalized = [min(y_data_normalized), max(y_data_normalized) - min(y_data_normalized)]
                popt_normalized, pcov_normalized = curve_fit(
                    lambda x, A0, A: gaussian_cdf_fixed(x, A0, A, xc_fitted_local, w_fitted_local),
                    x_data,
                    y_data_normalized,
                    p0=initial_guesses_normalized
                )
                A0_fitted_normalized, A_fitted_normalized = popt_normalized

                # Append successful fitting result
                results.append({
                    'File': filename,
                    'Gene_Name': gene_name,
                    'Transcript_ID': transcript_id,
                    'Global_Frequency': global_frequency,
                    'Present_in_Transcripts': present_in_transcripts,
                    'Sum or Fitted A (Abundance)': '{:.2f}'.format(A_fitted_normalized),
                    'Fixed Mean (xc)': '{:.2f}'.format(xc_fitted_local),
                    'Fixed Standard Deviation (w)': '{:.2f}'.format(w_fitted_local),
                    'Report': 'OK'
                })
            else:
                sum_normalized_kmer_count = y_data_normalized.sum()
                results.append({
                    'File': filename,
                    'Gene_Name': gene_name,
                    'Transcript_ID': transcript_id,
                    'Global_Frequency': global_frequency,
                    'Present_in_Transcripts': present_in_transcripts,
                    'Sum or Fitted A (Abundance)': '{:.2f}'.format(sum_normalized_kmer_count),
                    'Fixed Mean (xc)': 'N/A',
                    'Fixed Standard Deviation (w)': 'N/A',
                    'Report': 'Mean not more than SD'
                })
        except RuntimeError as e:
            error_message = str(e)
            if "Optimal parameters not found" in error_message:
                # Output the sum of individual values from the 'Normalized_K-mer_Count' column
                sum_normalized_kmer_count = gc_content_data['Normalized_K-mer_Count'].sum()
                results.append({
                    'File': filename,
                    'Gene_Name': gene_name,
                    'Transcript_ID': transcript_id,
                    'Global_Frequency': global_frequency,
                    'Present_in_Transcripts': present_in_transcripts,
                    'Sum or Fitted A (Abundance)': '{:.2f}'.format(sum_normalized_kmer_count),
                    'Fixed Mean (xc)': 'N/A',
                    'Fixed Standard Deviation (w)': 'N/A',
                    'Report': 'Fit Failed - Optimal parameters not found'
                })
            else:
                # For other types of RuntimeError (not related to optimization), you can decide how to handle these cases
                results.append({
                    'File': filename,
                    'Gene_Name': gene_name,
                    'Transcript_ID': transcript_id,
                    'Global_Frequency': global_frequency,
                    'Present_in_Transcripts': present_in_transcripts,
                    'Sum or Fitted A (Abundance)': 'N/A',
                    'Fixed Mean (xc)': 'N/A',
                    'Fixed Standard Deviation (w)': 'N/A',
                    'Report': f'Fit Failed - {error_message}'
                })
                

# Create results DataFrame
results_df = pd.DataFrame(results)

# Save to CSV file specified by the command-line argument
results_df.to_csv(args.output, index=False)

# Print the results
#print(results_df)
