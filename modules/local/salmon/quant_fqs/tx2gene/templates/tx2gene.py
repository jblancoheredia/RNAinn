#!/usr/bin/env python3

# Written by Lorena Pantano with subsequent reworking by Jonathan Manning.
# Modified by Juan Blanco Heredia blancoj@mskcc.org for mskcc/RNAinn pipeline. 
# Released under the MIT license.

import os
import re
import glob
import logging
import argparse
import platform
from typing import Dict
from collections.abc import Set
from collections import Counter, defaultdict, OrderedDict

# Configure logging
logging.basicConfig(format="%(name)s - %(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.

    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.

    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str

def read_top_transcripts(quant_file: str) -> Set[str]:
    """
    Read the top 100 transcripts from the quantification file.

    Parameters:
    quant_file (str): Path to the quantification file.

    Returns:
    set: A set containing the top 100 transcripts.
    """
    try:
        with open(quant_file, "r") as file_handle:
            # Read the file and extract the top 100 transcripts
            return {line.split()[0] for i, line in enumerate(file_handle) if i > 0 and i <= 100}
    except FileNotFoundError:
        # Log an error and raise if the quant file does not exist
        logger.error(f"Quantification file not found: {quant_file}")
        raise FileNotFoundError(f"Quantification file not found: {quant_file}")

def discover_transcript_attribute(gtf_file: str, transcripts: Set[str]) -> str:
    """
    Discover the attribute in the GTF that corresponds to transcripts, prioritizing 'transcript_id'.

    Parameters:
    gtf_file (str): Path to the GTF file.
    transcripts (Set[str]): A set of transcripts to match in the GTF file.

    Returns:
    str: The attribute name that corresponds to transcripts in the GTF file.
    """

    votes = Counter()
    with open(gtf_file) as inh:
        # Read GTF file, skipping header lines
        for line in filter(lambda x: not x.startswith("#"), inh):
            cols = line.split("\\t")

            # Use regular expression to correctly split the attributes string
            attributes_str = cols[8]
            attributes = dict(re.findall(r'(\\S+) "(.*?)(?<!\\\\)";', attributes_str))

            votes.update(key for key, value in attributes.items() if value in transcripts)

    if not votes:
        # Error out if no matching attribute is found
        logger.error("No attribute in GTF matching transcripts")

    # Check if 'transcript_id' is among the attributes with the highest votes
    if "transcript_id" in votes and votes["transcript_id"] == max(votes.values()):
        logger.info("Attribute 'transcript_id' corresponds to transcripts.")
        return "transcript_id"

    # If 'transcript_id' isn't the highest, determine the most common attribute that matches the transcripts
    if votes:
        attribute, _ = votes.most_common(1)[0]
        logger.info(f"Attribute '{attribute}' corresponds to transcripts.")
        return attribute
    else:
        # Fallback to transcript_id if no matches found
        logger.warning("No matching attribute found, defaulting to 'transcript_id'")
        return "transcript_id"

def parse_attributes(attributes_text: str) -> Dict[str, str]:
    """
    Parse the attributes column of a GTF file.

    Parameters:
    attributes_text: The attributes column as a string.

    Returns:
    attr_dict: A dictionary of the attributes.
    """
    # Split the attributes string by semicolon and strip whitespace
    attributes = attributes_text.strip().split(";")
    attr_dict = OrderedDict()

    # Iterate over each attribute pair
    for attribute in attributes:
        # Skip empty attributes
        if not attribute.strip():
            continue
        # Split the attribute into key and value, ensuring there are two parts
        parts = attribute.strip().split(" ", 1)
        if len(parts) == 2:
            key, value = parts
            # Remove any double quotes from the value
            value = value.replace('"', "")
            attr_dict[key] = value

    return attr_dict

def map_transcripts_to_gene(
    quant_type: str,
    gtf_file: str,
    quant_file: str,
    gene_id: str,
    extra_id_field: str,
    output_file: str,
) -> bool:
    """
    Map transcripts to gene names and write the output to a file.

    Parameters:
    quant_type (str): The quantification method used (e.g., 'salmon').
    gtf_file (str): Path to the GTF file.
    quant_file (str): Path to the quantification file.
    gene_id (str): The gene ID attribute in the GTF file.
    extra_id_field (str): Additional ID field in the GTF file.
    output_file (str): The output file path.

    Returns:
    bool: True if the operation was successful, False otherwise.
    """
    # Read the top transcripts based on quantification type
    transcripts = read_top_transcripts(quant_file)
    # Discover the attribute that corresponds to transcripts in the GTF
    transcript_attribute = discover_transcript_attribute(gtf_file, transcripts)

    # Open GTF and output file to write the mappings
    # Initialize the set to track seen combinations
    seen = set()

    with open(gtf_file) as inh, open(output_file, "w") as output_handle:
        output_handle.write(f"{transcript_attribute}\\t{gene_id}\\t{extra_id_field}\\n")
        # Parse each line of the GTF, mapping transcripts to genes
        for line in filter(lambda x: not x.startswith("#"), inh):
            cols = line.split("\\t")
            if len(cols) < 9:
                continue

            attr_dict = parse_attributes(cols[8])
            if gene_id in attr_dict and transcript_attribute in attr_dict:
                # Create a unique identifier for the transcript-gene combination
                transcript_gene_pair = (
                    attr_dict[transcript_attribute],
                    attr_dict[gene_id],
                )

                # Check if the combination has already been seen
                if transcript_gene_pair not in seen:
                    # If it's a new combination, write it to the output and add to the seen set
                    extra_id = attr_dict.get(extra_id_field, attr_dict[gene_id])
                    output_handle.write(f"{attr_dict[transcript_attribute]}\\t{attr_dict[gene_id]}\\t{extra_id}\\n")
                    seen.add(transcript_gene_pair)

    return True

# Main execution
if __name__ == "__main__":
    prefix = "${meta.id}"
    # Find the appropriate quantification file
    if "${quant_type}" == "salmon":
        quant_file = "${prefix}_quant.sf"
    else:  # kallisto
        quant_file = "${prefix}_abundance.tsv"

    if not os.path.exists(quant_file):
        logger.error(f"Quantification file not found: {quant_file}")
        raise FileNotFoundError(f"Quantification file not found: {quant_file}")
    
    if not map_transcripts_to_gene("${quant_type}", "${gtf}", quant_file, "${id}", "${extra}", f"{prefix}_tx2gene.tsv"):
        logger.error("Failed to map transcripts to genes.")

    # Write the versions
    versions_this_module = {}
    versions_this_module["${task.process}"] = {"python": platform.python_version()}
    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions_this_module))
