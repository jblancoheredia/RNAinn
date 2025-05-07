#!/usr/bin/env python3

import argparse
from Bio import SeqIO

def filter_contigs(input_file, output_file, min_length=100, min_coverage=2):
    """Filter FASTA contigs based on minimum length and coverage requirements."""
    def parse_header(header):
        parts = header.split('_')
        length = int(parts[3])
        coverage = float(parts[5])
        return length, coverage
    
    filtered_records = []
    with open(input_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            length, coverage = parse_header(record.description)
            if length >= min_length and coverage >= min_coverage:
                filtered_records.append(record)
    
    SeqIO.write(filtered_records, output_file, "fasta")
    print(f"Filtered {len(filtered_records)} contigs")

def main():
    parser = argparse.ArgumentParser(description='Filter SPADES contigs by length and coverage.')
    parser.add_argument('input_file', help='Input FASTA file')
    parser.add_argument('output_file', help='Output filtered FASTA file')
    parser.add_argument('-l', '--min_length', type=int, default=100,
                        help='Minimum contig length (default: 100)')
    parser.add_argument('-c', '--min_coverage', type=float, default=2.0,
                        help='Minimum coverage (default: 2.0)')
    
    args = parser.parse_args()
    filter_contigs(args.input_file, args.output_file, args.min_length, args.min_coverage)

if __name__ == "__main__":
    main()
