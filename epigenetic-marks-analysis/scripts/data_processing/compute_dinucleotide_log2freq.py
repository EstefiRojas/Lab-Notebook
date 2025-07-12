import argparse
import numpy as np

def calculate_log2_ratios(sequence):
    """Calculates log2 ratios of dinucleotide frequencies to corresponding nucleotide frequencies.

    Args:
        sequence (str): The DNA sequence to analyze.

    Returns:
        dict: A dictionary containing log2 ratios for each dinucleotide.
    """

    dinucleotides = ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]
    observed_counts = {dinuc: 0 for dinuc in dinucleotides}
    base_counts = {"A": 0, "C": 0, "G": 0, "T": 0}

    # Count dinucleotides and individual bases
    for i in range(len(sequence) - 1):
        dinuc = sequence[i:i+2]
        observed_counts[dinuc] += 1
    for base in sequence:
        base_counts[base] += 1

    # Calculate log2 ratios
    log2_ratios = {}
    for dinuc in dinucleotides:
        nuc1, nuc2 = dinuc[0], dinuc[1]
        dinuc_freq = observed_counts[dinuc] / (len(sequence) - 1)
        nuc1_freq = base_counts[nuc1] / len(sequence)
        nuc2_freq = base_counts[nuc2] / len(sequence)

        # Handle potential zero frequencies using a small pseudocount
        pseudocount = 1e-4
        dinuc_freq += pseudocount
        nuc1_freq += np.sqrt(pseudocount)
        nuc2_freq += np.sqrt(pseudocount)

        ratio = dinuc_freq / (nuc1_freq * nuc2_freq)
        log2_ratios[dinuc] = np.log2(ratio)

    return log2_ratios


def main():
    """Handles command-line argument parsing and script execution."""

    parser = argparse.ArgumentParser(description="Calculate log2 ratios of dinucleotide to nucleotide frequencies.")
    parser.add_argument("sequence", help="The DNA sequence to analyze.")
    args = parser.parse_args()

    try:
        # Validate input (basic check)
        if not set(args.sequence.upper()).issubset({"A", "C", "G", "T"}):
            raise ValueError("Invalid sequence. Only A, C, G, and T are allowed.")

        log2_ratios = calculate_log2_ratios(args.sequence.upper())
        # Create a comma-separated string for output
        output_string = ""
        for dinuc, ratio in list(log2_ratios.items())[:-1]:  # Iterate up to the second-to-last item
            output_string += f"{ratio:.4f},"  
        
        # Add the last item without a trailing comma
        last_dinuc, last_ratio = list(log2_ratios.items())[-1]
        output_string += f"{last_ratio:.4f}"  
        
        print(output_string)  # Print the formatted string

    except ValueError as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()

