import argparse
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import norm
import os
import re


# Function to calculate the GC content percentage of a sequence
def calculate_gc_content(seq):
    gc_count = seq.count('G') + seq.count('C')
    gc_content_percentage = (gc_count / len(seq)) * 100
    return gc_content_percentage


# Adjusted Gaussian CDF function according to the specified equation
def gaussian_cdf(x, A0, A, xc, w):
    return A0 + A * norm.cdf((x - xc) / w)


def main(input_folder, output_file_path, row_threshold):
    # Get the full path of the input folder
    input_folder_path = os.path.abspath(input_folder)

    # Find all CSV files in the input folder
    csv_files = [file for file in os.listdir(input_folder_path) if file.endswith('.csv')]

    # Open the output file and write the header
    with open(output_file_path, 'w') as f_out:
        f_out.write('Gene,Transcript ID,Abundance RPKM,Mean,SD\n')  # Write the header

    # Process each CSV file
    for file_name in csv_files:
        # Read the CSV file without a header
        df = pd.read_csv(os.path.join(input_folder_path, file_name), header=None)

        # Convert 'Frequency' column to float
        df[1] = pd.to_numeric(df[1], errors='coerce')

        # Drop NA values that could not be converted to float
        df.dropna(subset=[1], inplace=True)

        # Add a new column to the DataFrame for GC content
        df['GC Content (%)'] = df[0].apply(calculate_gc_content)

        # Group by 'GC Content (%)' and aggregate 'Frequency' column as sum
        gc_content_grouped = df.groupby('GC Content (%)')[1].agg(['sum']).reset_index()

        # Skip if the number of rows is less than the threshold
        if len(gc_content_grouped) < row_threshold:
            print(
                f"Skipping file {file_name}: number of data points ({len(gc_content_grouped)}) below threshold ({row_threshold}).")
            continue

        # Sort the grouped data by 'GC Content (%)' in ascending order
        gc_content_grouped.sort_values('GC Content (%)', inplace=True)

        try:
            # Fit the Gaussian CDF to the accumulated sum of frequencies
            y_data = gc_content_grouped['sum'].cumsum()
            initial_guess = [0, y_data.iloc[-1], gc_content_grouped['GC Content (%)'].mean(),
                             gc_content_grouped['GC Content (%)'].std()]
            popt, _ = curve_fit(gaussian_cdf, gc_content_grouped['GC Content (%)'], y_data, p0=initial_guess,
                                maxfev=9000)

            # Extract gene name and transcript ID from the filename using regular expressions
            match = re.match(r'(\w+)_(ENST\d+\.\d+)_.*', file_name)
            if match:
                gene_name, transcript_id = match.groups()
            else:
                print(f"Could not parse gene name and transcript ID from {file_name}")
                continue

            # Format output to include gene name and transcript ID, and fitting parameters rounded to two decimal places
            result_string = f"{gene_name},{transcript_id},{popt[1]:.2f},{popt[2]:.2f},{popt[3]:.2f}\n"

            with open(output_file_path, "a") as f_out:
                f_out.write(result_string)
        except RuntimeError as e:
            print(f"Warning: Fitting did not converge for {file_name}. Error: {e}")

    print(f"Results have been written to {output_file_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze GC content and fit Gaussian CDF.")
    parser.add_argument('--input', type=str, required=True, help="Path to the input folder containing the CSV files.")
    parser.add_argument('--output', type=str, required=True,
                        help="Path and name of the output file to save the results.")
    parser.add_argument('--threshold', type=int, default=10,
                        help="Minimum number of data points required for fitting. Default is 10.")
    args = parser.parse_args()

    # Ensure the output directory exists
    output_directory = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(output_directory, exist_ok=True)

    # Run the main function with additional row_threshold argument
    main(args.input, args.output, args.threshold)
