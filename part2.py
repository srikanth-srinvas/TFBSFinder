import argparse
from collections import Counter
import os
import logging
import matplotlib.pyplot as plt
import seaborn as sns

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def count_nucleotides(sequence):
    """
    Count the occurrences of each nucleotide in a sequence.
    """
    if not sequence:
        logging.warning("Empty sequence received.")
        return {}

    try:
        counts = {nucleotide: sequence.count(nucleotide) for nucleotide in "ACGT"}
        logging.info(f"Nucleotide counts for sequence: {counts}")
        return counts
    except Exception as e:
        logging.error(f"Error while counting nucleotides: {e}")
        raise


def analyze_tfbs_output(input_file):
    """Analyze TFBS output to find top sequences and count sequences without matches."""
    logging.info(f"Analyzing TFBS results from file: {input_file}")
    tfbs_counts = Counter()
    no_tfbs_count = 0

    try:
        with open(input_file, "r") as file:
            next(file)  # Skip header
            for line in file:
                fields = line.strip().split("\t")
                if len(fields) < 2:
                    logging.error(f"Invalid line format: {line}")
                    raise ValueError(f"Invalid line format: {line}")

                sequence_id, *tfbs_fields = fields  # Capture extra TFBS values
                tfbs = ",".join(
                    tfbs_fields
                )  # Combine multiple TFBS fields into one string

                if tfbs == "None":
                    no_tfbs_count += 1
                else:
                    tfbs_list = tfbs.split(",")
                    tfbs_counts.update(tfbs_list)
    except FileNotFoundError:
        logging.error(f"File {input_file} not found.")
        raise
    except Exception as e:
        logging.error(f"Error while analyzing TFBS output: {e}")
        raise

    # Log results
    logging.info(f"Total sequences without matches: {no_tfbs_count}")
    logging.info(f"Top 10 TFBS counts: {tfbs_counts.most_common(10)}")

    # Get top 10 most common TFBS sequences
    top_tfbs = tfbs_counts.most_common(10)
    return top_tfbs, no_tfbs_count


def plot_top_tfbs(top_tfbs, output_file):
    """Plot the top 10 most common TFBS sequences and save as a PNG."""
    if not top_tfbs:
        logging.warning("No TFBS data to plot.")
        return

    sequences, counts = zip(
        *top_tfbs
    )  # Unpack the top TFBS into sequences and their counts
    logging.info("Creating plot for the top TFBS sequences.")

    try:
        # Create a bar plot using Seaborn
        plt.figure(figsize=(10, 6))
        sns.barplot(
            x=list(sequences),
            y=list(counts),
            hue=list(sequences),
            palette="viridis",
            legend=False,
        )

        # Add labels and title
        plt.xlabel("TFBS Sequence", fontsize=12)
        plt.ylabel("Count", fontsize=12)
        plt.title("Top 10 Most Common TFBS Sequences", fontsize=14)
        plt.xticks(
            rotation=45, ha="right"
        )  # Rotate x-axis labels for better readability
        plt.tight_layout()  # Ensure everything fits nicely on the plot

        # Save the plot to a file
        plt.savefig(output_file)
        logging.info(f"Plot saved to: {output_file}")
        plt.close()
    except Exception as e:
        logging.error(f"Error while creating or saving the plot: {e}")
        raise


def save_summary_to_file(top_tfbs, no_tfbs_count, output_file):
    """Save the summary table to a tab-separated text file."""
    try:
        with open(output_file, "w") as file:
            file.write("Top 10 Most Common TFBS Sequences:\n")
            for tfbs, count in top_tfbs:
                file.write(f"{tfbs}\t{count}\n")
            file.write(f"\nNumber of sequences without any TFBS:\t{no_tfbs_count}\n")
        logging.info(f"Summary saved to: {output_file}")
    except Exception as e:
        logging.error(f"Error while saving summary to file: {e}")
        raise


def main():
    parser = argparse.ArgumentParser(description="Analyze TFBS results.")
    parser.add_argument("input_file", help="Input file from Part 1.")
    parser.add_argument(
        "--output_dir",
        default="TFBS_results",
        help="Directory to save output files (default: 'TFBS_results').",
    )
    args = parser.parse_args()

    # Ensure output directory exists
    try:
        os.makedirs(args.output_dir, exist_ok=True)
        logging.info(f"Output directory created or already exists: {args.output_dir}")
    except Exception as e:
        logging.error(f"Error while creating output directory: {e}")
        raise

    try:
        # Analyze the TFBS file
        top_tfbs, no_tfbs_count = analyze_tfbs_output(args.input_file)

        # Output file paths
        summary_file = os.path.join(args.output_dir, "TFBS_summary.txt")
        plot_file = os.path.join(args.output_dir, "TFBS_plot.png")

        # Save the summary table and plot
        save_summary_to_file(top_tfbs, no_tfbs_count, summary_file)
        plot_top_tfbs(top_tfbs, plot_file)

        # Print to console
        logging.info("Analysis complete.")
        print("Top 10 most common TFBS sequences:")
        for tfbs, count in top_tfbs:
            print(f"{tfbs}: {count}")

        print(f"\nNumber of sequences without any TFBS: {no_tfbs_count}")
        print(f"\nSummary table saved to: {summary_file}")
        print(f"Plot saved to: {plot_file}")
    except Exception as e:
        logging.error(f"Error during analysis: {e}")
        raise


if __name__ == "__main__":
    main()
