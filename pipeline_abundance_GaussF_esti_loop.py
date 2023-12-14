import argparse
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm
import os


# Function to calculate the GC content percentage of a sequence
def calculate_gc_content(seq):
    gc_count = seq.count('G') + seq.count('C')
    gc_content_percentage = (gc_count / len(seq)) * 100
    return gc_content_percentage


# Adjusted Gaussian CDF function according to the specified equation
def gaussian_cdf(x, A0, A, xc, w):
    return A0 + A * norm.cdf((x - xc) / w)


# Function to use the fixed mean and standard deviation for curve fitting
def gaussian_cdf_fixed(x, A0, A, fixed_xc, fixed_w):
    return gaussian_cdf(x, A0, A, fixed_xc, fixed_w)


def main(input_folder, output_file_path):
    # Get the full path of the input folder
    input_folder_path = os.path.abspath(input_folder)

    # Find all CSV files in the input folder
    csv_files = [file for file in os.listdir(input_folder_path) if file.endswith('.csv')]

    # Process each CSV file
    for file_name in csv_files:
        # Read the CSV file without a header
        df = pd.read_csv(os.path.join(input_folder_path, file_name), header=None)

        # If the first row is non-numeric (i.e., names of the columns), skip it
        if any(isinstance(val, str) and val.strip().isalpha() for val in df.iloc[0]):
            df = df.iloc[1:].reset_index(drop=True)

        # Assign appropriate column names after reading the CSV file
        df.columns = ['K-mer', 'Frequency']

        # Convert 'Frequency' column to float
        df['Frequency'] = pd.to_numeric(df['Frequency'], errors='coerce')

        # Drop NA values that could not be converted to float
        df.dropna(subset=['Frequency'], inplace=True)

        # Add a new column to the DataFrame for GC content
        df['GC Content (%)'] = df['K-mer'].apply(calculate_gc_content)

        # Group by 'GC Content (%)' and aggregate 'Frequency' column as sum and count
        gc_content_grouped = df.groupby('GC Content (%)')['Frequency'].agg(['sum', 'size']).reset_index()

        # Sort the grouped data by 'GC Content (%)' in ascending order
        gc_content_grouped.sort_values('GC Content (%)', inplace=True)

        # Calculate the accumulated 'Sum Accumulated' and 'Size Accumulated' values
        gc_content_grouped['Sum Accumulated'] = gc_content_grouped['sum'].cumsum()
        gc_content_grouped['Size Accumulated'] = gc_content_grouped['size'].cumsum()

        # Fit the Gaussian CDF to the Size Accumulated data
        initial_guess = [0, gc_content_grouped['Size Accumulated'].iloc[-1],
                         gc_content_grouped['GC Content (%)'].mean(), gc_content_grouped['GC Content (%)'].std()]
        popt, pcov = curve_fit(gaussian_cdf,
                               gc_content_grouped['GC Content (%)'],
                               gc_content_grouped['Size Accumulated'],
                               p0=initial_guess)
        A0_fit, A_fit, xc_fit, w_fit = popt

        # Use the fitted xc and w as fixed values to fit the Sum Accumulated data
        initial_guess_sum = [0, gc_content_grouped['Sum Accumulated'].iloc[-1]]
        popt_sum, pcov_sum = curve_fit(
            lambda x, A0, A: gaussian_cdf_fixed(x, A0, A, xc_fit, w_fit),
            gc_content_grouped['GC Content (%)'],
            gc_content_grouped['Sum Accumulated'],
            p0=initial_guess_sum
        )
        A0_fit_sum, A_fit_sum = popt_sum

        # Format the output result
        result_string = f"{os.path.splitext(file_name)[0]}: xc = {xc_fit}, w = {w_fit}, A = {A_fit_sum}\n"

        # Write the result for this file to the output file
        with open(output_file_path, "a") as f:
            f.write(result_string)

    print(f"Results have been written to {output_file_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze GC content and fit Gaussian CDF.")
    parser.add_argument(
        '--input',
        type=str,
        default="unique_frequency_csv",
        help="Path to the input folder containing the CSV files. Default is 'unique_frequency_csv'."
    )
    parser.add_argument(
        '--output',
        type=str,
        default="results.txt",
        help="Path and name of the output file to save the results. Default is 'results.txt'."
    )
    args = parser.parse_args()

    # Ensure directory exists for the output file
    output_directory = os.path.dirname(os.path.abspath(args.output))
    if output_directory and not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Check if the input folder exists
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input folder '{args.input}' not found.")

    # Clear previous content before writing new results
    if os.path.exists(args.output):
        open(args.output, 'w').close()

    main(args.input, args.output)
