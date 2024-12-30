import unittest
from unittest.mock import patch, mock_open
import os
from collections import Counter
import tempfile
from part2 import analyze_tfbs_output, plot_top_tfbs, save_summary_to_file, count_nucleotides

class TestPart2(unittest.TestCase):

    def setUp(self):
        # Create a temporary directory for output files
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_dir = self.temp_dir.name

        # Mock input file content for TFBS analysis
        self.mock_tfbs_input = (
            "Sequence_ID\tTFBS\n"
            "seq1\tGGGAAATTCC\n"
            "seq2\tGGGGAATTCC\n"
            "seq3\tNone\n"
            "seq4\tGGGAAATTCC,GGGGAATTCC\n"
            "seq5\tNone\n"
        )
        self.mock_tfbs_file = os.path.join(self.output_dir, "mock_tfbs_input.txt")
        with open(self.mock_tfbs_file, "w") as f:
            f.write(self.mock_tfbs_input)

    def tearDown(self):
        # Clean up temporary directory
        self.temp_dir.cleanup()

    def test_analyze_tfbs_output(self):
        """Test the analyze_tfbs_output function."""
        top_tfbs, no_tfbs_count = analyze_tfbs_output(self.mock_tfbs_file)

        # Check top TFBS
        expected_tfbs = Counter({"GGGAAATTCC": 2, "GGGGAATTCC": 2}).most_common(10)
        self.assertEqual(top_tfbs, expected_tfbs)

        # Check no TFBS count
        self.assertEqual(no_tfbs_count, 2)

    def test_plot_top_tfbs(self):
        """Test the plot_top_tfbs function."""
        top_tfbs = [("GGGAAATTCC", 5), ("GGGGAATTCC", 3)]
        plot_file = os.path.join(self.output_dir, "TFBS_plot.png")
        
        plot_top_tfbs(top_tfbs, plot_file)

        # Check if the plot file was created
        self.assertTrue(os.path.exists(plot_file))

    def test_save_summary_to_file(self):
        """Test the save_summary_to_file function."""
        top_tfbs = [("GGGAAATTCC", 5), ("GGGGAATTCC", 3)]
        no_tfbs_count = 2
        summary_file = os.path.join(self.output_dir, "TFBS_summary.txt")
        
        save_summary_to_file(top_tfbs, no_tfbs_count, summary_file)

        # Check if the summary file has the expected content
        with open(summary_file, "r") as file:
            content = file.read()
        
        expected_content = (
            "Top 10 Most Common TFBS Sequences:\n"
            "GGGAAATTCC\t5\n"
            "GGGGAATTCC\t3\n"
            "\nNumber of sequences without any TFBS:\t2\n"
        )
        self.assertEqual(content, expected_content)

    def test_count_nucleotides_normal_sequence(self):
        """Test count_nucleotides with a normal sequence."""
        result = count_nucleotides('ATCGATCG')
        self.assertEqual(result, {'A': 2, 'C': 2, 'G': 2, 'T': 2})

    def test_count_nucleotides_empty_sequence(self):
        """Test count_nucleotides with an empty sequence."""
        with self.assertLogs(level='WARNING'):
            result = count_nucleotides('')
        self.assertEqual(result, {})

    def test_count_nucleotides_large_sequence(self):
        """Test count_nucleotides with a large sequence."""
        sequence = 'A' * 1000000 + 'C' * 1000000 + 'G' * 1000000 + 'T' * 1000000
        result = count_nucleotides(sequence)
        self.assertEqual(result, {'A': 1000000, 'C': 1000000, 'G': 1000000, 'T': 1000000})

    def test_count_nucleotides_invalid_characters(self):
        """Test count_nucleotides with invalid characters."""
        result = count_nucleotides('ATCGXYZ')
        self.assertEqual(result, {'A': 1, 'C': 1, 'G': 1, 'T': 1})  # Ignores non-ACGT characters


if __name__ == "__main__":
    unittest.main()