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

# Function to calculate the Gaussian CDF
def gaussian_cdf(x, A0, amplitude, mean, stddev):
    return A0 + amplitude * norm.cdf(x, mean, stddev)

# Function with fixed mean and stddev
def gaussian_cdf_fixed(x, A0, amplitude, fixed_mean, fixed_stddev):
    return gaussian_cdf(x, A0, amplitude, fixed_mean, fixed_stddev)

def main(input_folder, output_file_path):
    # Get the full path of the input folder
    input_folder_path = os.path.abspath(input_folder)

    # Find all CSV files in the input folder
    csv_files = [file for file in os.listdir(input_folder_path) if file.endswith('.csv')]

    # Initialize a dictionary to store amplitude results
    amplitude_results = {}

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

        # Add a new column to the dataframe for GC content
        df['GC Content (%)'] = df['K-mer'].apply(calculate_gc_content)

        # Group by 'GC Content (%)' and sum up the 'Frequency' column
        gc_content_grouped = df.groupby('GC Content (%)')['Frequency'].sum().reset_index()

        # Sort the grouped data by 'GC Content (%)' in ascending order
        gc_content_grouped.sort_values('GC Content (%)', inplace=True)

        # Calculate the cumulative frequency
        gc_content_grouped['Accumulated Frequency'] = gc_content_grouped['Frequency'].cumsum()

        # Perform the Gaussian CDF fit assuming fixed values for mean and stddev
        initial_mean = np.mean(gc_content_grouped['GC Content (%)'])
        initial_stddev = np.std(gc_content_grouped['GC Content (%)'])

        x_data_cdf = gc_content_grouped['GC Content (%)']
        y_data_cdf = gc_content_grouped['Accumulated Frequency']

        # Perform the Gaussian CDF fit using fixed values for mean and stddev
        initial_guess_cdf = [1, max(y_data_cdf)]  # Baseline offset and amplitude

        # Fitting the Gaussian CDF with fixed mean and stddev
        popt_cdf, _ = curve_fit(
            lambda x, A0, amplitude: gaussian_cdf_fixed(x, A0, amplitude, initial_mean, initial_stddev),
            x_data_cdf,
            y_data_cdf,
            p0=initial_guess_cdf,
            bounds=(0, [np.inf, np.inf])
        )

        # Collect the amplitude_fit_cdf value
        A0_fit_cdf, amplitude_fit_cdf = popt_cdf
        amplitude_results[file_name] = amplitude_fit_cdf

    # Write the amplitude_fit_cdf results for each file to the specified output file
    with open(output_file_path, "w") as f:
        for file_name, amplitude in amplitude_results.items():
            f.write(f"{file_name}: {amplitude}\n")

    print(f"Amplitude fit CDF results have been saved to {output_file_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze GC content and fit Gaussian CDF.")
    parser.add_argument(
        '--input',
        type=str,
        default="unique_frequence_csv",
        help="Path to the input folder containing the CSV files. Default is 'unique_frequence_csv'."
    )
    parser.add_argument(
        '--output',
        type=str,
        default="amplitude_fit_cdf_results.txt",
        help="Path and name of the output file to save the results. Default is 'amplitude_fit_cdf_results.txt'."
    )
    args = parser.parse_args()

    # Check if the directory for the output file exists; if not, create it
    output_directory = os.path.dirname(os.path.abspath(args.output))
    if output_directory and not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Check if the input folder exists
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input folder '{args.input}' not found.")

    main(args.input, args.output)
