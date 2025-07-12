import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse


def plot_overlapping_reads(filename, start_line, num_reads):
  """Plots overlapping reads with lines connecting starts and ends.

  Args:
      filename: The name of the CSV file containing the data.
      num_reads: The number of reads to plot.
  """
  data = pd.read_csv(filename, sep="\t")

  # Select the first 10 rows
  data_subset = data.iloc[start_line:start_line+num_reads]  

  # Extract the columns as NumPy arrays 
  starts = data_subset["Start"].to_numpy()
  ends = data_subset["End"].to_numpy()

  startsMatch = data_subset["StartMatch"].to_numpy()
  endsMatch = data_subset["EndMatch"].to_numpy()

  # Generate row numbers for the Y-axis 
  row_numbers = np.arange(1, len(starts) + 1)

  # Follow the rest of the plotting code
  fig, ax = plt.subplots(figsize=(10, 6))

  # Plot a scatter plot
  ax.scatter(starts, row_numbers, color='blue', alpha=0.5)
  ax.scatter(ends, row_numbers, color='blue', alpha=0.5)

  ax.scatter(startsMatch, row_numbers, color='green', alpha=0.5)
  ax.scatter(endsMatch, row_numbers, color='green', alpha=0.5)


  # Draw lines between corresponding points
  for start, end, row in zip(starts, ends, row_numbers):
      ax.plot([start, end], [row, row], color='gray', linestyle='-', linewidth=0.5)


  # Draw lines between corresponding points
  for start, end, row in zip(startsMatch, endsMatch, row_numbers):
      ax.plot([start, end], [row, row], color='purple', linestyle='-', linewidth=0.5)


  # Set axis labels and title
  ax.set_xlabel('Position')
  ax.set_ylabel('Read No')
  ax.set_title('Overlapping reads visualization')

  # Show the plot
  plt.tight_layout()
  plt.show()


if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser(description='Plot overlapping reads with connecting lines.')
    parser.add_argument('start_line', type=int, help='Line number from where to start reading the file.')
    parser.add_argument('num_reads', type=int, help='Number of reads to plot from the start line.')
    parser.add_argument('--filename', default='overlaps.csv', help='Name of the CSV file.')
    args = parser.parse_args()

    # Call the plotting function with command-line arguments
    plot_overlapping_reads(args.filename, args.start_line, args.num_reads)
