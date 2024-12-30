import os
import unittest
import logging
from part1 import (
    validate_consensus,
    build_regex_from_consensus,
    process_fasta_file,
    parse_fasta,
)

# Configure logging for tests
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


class TestPart1(unittest.TestCase):
    def setUp(self):
        self.test_fasta = "test_sequences.fasta"
        self.test_consensus = "GGGRNWYYCC"
        self.output_file = "test_results.txt"
        self.empty_fasta = "test_empty.fasta"
        self.malformed_fasta = "test_malformed.fasta"

        # Create a mock FASTA file
        with open(self.test_fasta, "w") as f:
            f.write(">seq1\nATGGGAAATTCCGGGAAATTCC\n")
            f.write(">seq2\nCCGGGGAATTCCGTTCC\n")
            f.write(">seq3\nATGGTTCC\n")

        # Create an empty FASTA file
        open(self.empty_fasta, "w").close()

        # Create a malformed FASTA file
        with open(self.malformed_fasta, "w") as f:
            f.write("This is not a valid FASTA file.\n")

    def tearDown(self):
        # Cleanup test files
        for file in [
            self.test_fasta,
            self.output_file,
            self.empty_fasta,
            self.malformed_fasta,
        ]:
            if os.path.exists(file):
                os.remove(file)

    def test_validate_consensus(self):
        logging.info("Testing validate_consensus...")
        # Valid consensus sequence
        try:
            validate_consensus(self.test_consensus)
        except ValueError as e:
            self.fail(f"Validation failed unexpectedly: {e}")

        # Invalid consensus sequence
        with self.assertRaises(ValueError):
            validate_consensus("GGGZXYYCC")

    def test_build_regex_from_consensus(self):
        logging.info("Testing build_regex_from_consensus...")
        regex = build_regex_from_consensus(self.test_consensus)
        self.assertTrue(regex.match("GGGAAATTCC"))
        self.assertFalse(regex.match("TTTCCGGGAA"))

    def test_process_fasta_file(self):
        logging.info("Testing process_fasta_file...")
        process_fasta_file(self.test_fasta, self.test_consensus, self.output_file, "\t")
        with open(self.output_file, "r") as f:
            lines = f.readlines()
            self.assertEqual(lines[0], "Sequence_ID\tTFBS\n")
            # Ensure TFBS matches are present and accurate
            self.assertIn("seq1\tGGGAAATTCC\tGGGAAATTCC\n", lines)
            self.assertIn("seq2\tGGGGAATTCC\n", lines)
            self.assertIn("seq3\tNone\n", lines)

    def test_parse_fasta_valid(self):
        logging.info("Testing parse_fasta with valid FASTA...")
        sequences = parse_fasta(self.test_fasta)
        self.assertEqual(len(sequences), 3)
        self.assertIn("seq1", sequences)
        self.assertEqual(sequences["seq1"], "ATGGGAAATTCCGGGAAATTCC")

    def test_parse_fasta_empty(self):
        logging.info("Testing parse_fasta with empty FASTA...")
        with self.assertLogs(level="WARNING"):
            sequences = parse_fasta(self.empty_fasta)
        self.assertEqual(len(sequences), 0)

    def test_parse_fasta_malformed(self):
        logging.info("Testing parse_fasta with malformed FASTA...")
        with self.assertRaises(ValueError):
            parse_fasta(self.malformed_fasta)

    def test_parse_fasta_missing_file(self):
        logging.info("Testing parse_fasta with missing file...")
        with self.assertRaises(FileNotFoundError):
            parse_fasta("nonexistent_file.fasta")


if __name__ == "__main__":
    unittest.main()
