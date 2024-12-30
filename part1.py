import argparse
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# ASCII Art for Tool Name
ASCII_ART = r"""
 _________  ________  ______     ______       ___  _                 __                
|  _   _  ||_   __  ||_   _ \  .' ____ \    .' ..](_)               |  ]               
|_/ | | \_|  | |_ \_|  | |_) | | (___ \_|  _| |_  __   _ .--.   .--.| | .---.  _ .--.  
    | |      |  _|     |  __'.  _.____`.  '-| |-'[  | [ `.-. |/ /'`\' |/ /__\\[ `/'`\] 
   _| |_    _| |_     _| |__) || \____) |   | |   | |  | | | || \__/  || \__., | |     
  |_____|  |_____|   |_______/  \______.'  [___] [___][___||__]'.__.;__]'.__.'[___]    
                                                                                       

Extract transcription factor binding sites (TFBS) from DNA sequences based on IUPAC consensus sequences.
"""

# Map IUPAC consensus codes to possible bases
IUPAC_CODES = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
    'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
}

def parse_fasta(fasta_file):
    """
    Parse a FASTA file and return a dictionary of sequences.
    """
    sequences = {}
    try:
        with open(fasta_file, 'r') as file:
            current_seq_name = None
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    current_seq_name = line[1:]
                    sequences[current_seq_name] = ""
                elif current_seq_name:
                    sequences[current_seq_name] += line
                else:
                    logging.error("Malformed FASTA: Sequence data without header.")
                    raise ValueError("Malformed FASTA file.")
        if not sequences:
            logging.warning("The input FASTA file is empty.")
    except FileNotFoundError:
        logging.error(f"File {fasta_file} not found.")
        raise
    except Exception as e:
        logging.error(f"Error while parsing the FASTA file: {e}")
        raise
    return sequences

def build_regex_from_consensus(consensus):
    """
    Convert an IUPAC consensus sequence into a regex pattern.

    Args:
        consensus (str): Consensus sequence in IUPAC notation.

    Returns:
        re.Pattern: Compiled regex pattern corresponding to the consensus sequence.
    """
    regex_pattern = ''.join(IUPAC_CODES[base] for base in consensus)
    return re.compile(regex_pattern)

def validate_consensus(consensus):
    """
    Validate the IUPAC consensus sequence.

    Args:
        consensus (str): Consensus sequence in IUPAC notation.

    Raises:
        ValueError: If the sequence contains invalid characters.
    """
    for base in consensus:
        if base not in IUPAC_CODES:
            logging.error(f"Invalid base '{base}' in consensus sequence.")
            raise ValueError(f"Invalid base '{base}' in consensus sequence.")

def find_overlapping_matches(sequence, pattern):
    """
    Find all matches of a regex pattern in a given DNA sequence, including overlapping matches.

    Args:
        sequence (str): DNA sequence to search in.
        pattern (re.Pattern): Compiled regex pattern.

    Returns:
        list: List of starting positions (0-based) of all matches.
    """
    matches = []
    start = 0
    while start < len(sequence):
        match = pattern.search(sequence, start)
        if not match:
            break
        matches.append(sequence[match.start():match.end()])
        start = match.start() + 1  # Allow overlap by moving start forward by one
    return matches

def process_fasta_file(input_fasta, consensus, output_file, delimiter):
    """
    Process a FASTA file to find transcription factor binding sites (TFBS).

    Args:
        input_fasta (str): Path to the input FASTA file.
        consensus (str): Consensus sequence in IUPAC notation.
        output_file (str): Path to the output file.
        delimiter (str): Delimiter for the output file (e.g., tab or comma).
    """
    logging.info(f"Processing FASTA file: {input_fasta}")
    regex_pattern = build_regex_from_consensus(consensus)
    try:
        sequences = parse_fasta(input_fasta)
    except Exception as e:
        logging.error(f"Failed to parse FASTA file: {e}")
        raise
    
    with open(output_file, 'w') as output:
        # Write header to the output file
        header = f"Sequence_ID{delimiter}TFBS\n"
        output.write(header)
        
        for seq_id, seq in sequences.items():
            rev_comp = str(Seq(seq).reverse_complement())  # Reverse complement
            matches = find_overlapping_matches(seq, regex_pattern)
            matches += find_overlapping_matches(rev_comp, regex_pattern)
            
            tfbs_results = delimiter.join(matches) if matches else "None"
            output.write(f"{seq_id}{delimiter}{tfbs_results}\n")
    
    logging.info(f"Results written to {output_file}")

def main():
    """
    Main function to parse arguments and process the input FASTA file.
    """
    print(ASCII_ART)
    parser = argparse.ArgumentParser(
        usage="python3 part1.py input_fasta consensus [--output_file OUTPUT_FILE] [--delimiter {tab,comma}]",
        epilog="Example: python script.py input.fasta GGGRNWYYCC --output_file results.txt --delimiter comma"
    )
    parser.add_argument("input_fasta", help="Input FASTA file with DNA sequences.")
    parser.add_argument("consensus", help="Consensus sequence in IUPAC notation.")
    parser.add_argument(
        "--output_file",
        default="tfbs_results.txt",
        help="Output file for TFBS matches (default: tfbs_results.txt)."
    )
    parser.add_argument(
        "--delimiter",
        choices=["tab", "comma"],
        default="tab",
        help="Delimiter for the output file (default: tab)."
    )
    
    args = parser.parse_args()

    try:
        validate_consensus(args.consensus)
        delimiter = "\t" if args.delimiter == "tab" else ","
        process_fasta_file(args.input_fasta, args.consensus, args.output_file, delimiter)
    except ValueError as e:
        logging.error(f"Error: {e}")

if __name__ == "__main__":
    main()