#!/usr/bin/env python
# Check the samplesheet for valid inputs
# Checks for
# 1. Whether minimum required files exist
# 2. Whether provided sample IDs are unique
# 3. Whether appropriate parameter combinations were provided
#   - if coverage is to be determined from read alignments, reads must also be provided
#   - if a path to a coverage table is already provided in the samplesheet
#
# This script is based on the example at: https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv

import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat nf-core/autometa samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path: str) -> None:
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error: str, context: str = "Line", context_str: str = "") -> None:
    error_str = f"ERROR: Please check samplesheet -> {error}"
    if context != "" and context_str != "":
        error_str = f"ERROR: Please check samplesheet -> {error}\n{context.strip()}: '{context_str.strip()}'"
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:

    NOTE: 0 and 1 are converted to boolean values (false and true, respectively)
     *and* although for sample 2 reads were provided, the input is specifying to
    retrieve the coverage information from the assembly's contigs' headers.

    sample,assembly,fastq_1,fastq_2,coverage_tab,cov_from_assembly
    SAMPLE_1,ASSEMBLY_1,fwd_reads.fastq.gz,rev_reads.fastq.gz,,0
    SAMPLE_2,ASSEMBLY_2,fwd_reads.fastq.gz,rev_reads.fastq.gz,,1
    SAMPLE_3,ASSEMBLY_3,,,,1
    SAMPLE_4,ASSEMBLY_4,,,/path/to/coverage.tsv,0

    For an example see:
    https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
    """

    samples = {}
    min_col_names = ["sample", "assembly", "cov_from_assembly"]
    min_num_cols = len(min_col_names)
    req_header_cols = [
        "sample",
        "assembly",
        "fastq_1",
        "fastq_2",
        "coverage_tab",
        "cov_from_assembly",
    ]
    num_required_cols = len(req_header_cols)
    with open(file_in, "r") as fh:
        ## Check header
        header = fh.readline().strip()
        header_cols = [header_col.strip('"') for header_col in header.split(",")]
        if header_cols[:num_required_cols] != req_header_cols:
            print(
                f"ERROR: Please check samplesheet header -> {','.join(header_cols)} != {','.join(req_header_cols)}"
            )
            sys.exit(1)

        ## Check sample entries
        for line in fh:
            sample_cols = [
                sample_col.strip().strip('"') for sample_col in line.strip().split(",")
            ]

            # Check valid number of columns per row
            num_cols = len([sample_col for sample_col in sample_cols if sample_col])
            if num_cols < min_num_cols:
                print_error(
                    f"Invalid number of populated columns (minimum = {min_num_cols})!",
                    "Line",
                    line,
                )

            if len(sample_cols) != num_required_cols:
                print_error(
                    f"Invalid number of columns (required = {num_required_cols})! Check that you have the correct number of commas and retry",
                    "Line",
                    line,
                )
            ## Check sample name entries
            (
                sample,
                assembly,
                fastq_1,
                fastq_2,
                coverage_tab,
                cov_from_assembly,
            ) = sample_cols

            sample = sample.replace(" ", "_")
            if not sample:
                print_error("Sample entry has not been specified!", "Line", line)

            cov_from_assembly_choices = set(["0", "spades"])
            if not cov_from_assembly:
                print_error(
                    f"Must provide a value ({', '.join(cov_from_assembly_choices)}) for cov_from_assembly column!",
                    "Line",
                    line,
                )
            if cov_from_assembly not in cov_from_assembly_choices:
                print_error(
                    f"cov_from_assembly column must be one of {', '.join(cov_from_assembly_choices)}!",
                    "Line",
                    line,
                )

            ## Check FastQ file extension
            for fastq in [fastq_1, fastq_2]:
                if fastq:
                    if fastq.find(" ") != -1:
                        print_error("FastQ file contains spaces!", "Line", line)
                    if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                        print_error(
                            "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                            "Line",
                            line,
                        )

            # NOTE: An empty string will fail with the file(...) method for groovy.. So we pass in "0" here
            # to be checked later with file(...).exists()
            # Same goes for the fastq_* files
            coverage_tab = coverage_tab if coverage_tab else "0"
            ## Auto-detect paired-end/single-end
            ## [ assembly, single_end, fastq_1, fastq_2, coverage_tab, cov_from_assembly ]
            sample_info = []
            if fastq_1 and fastq_2:  ## Paired-end short reads
                sample_info = [
                    assembly,
                    "0",
                    fastq_1,
                    fastq_2,
                    coverage_tab,
                    cov_from_assembly,
                ]
            elif fastq_1 and not fastq_2:  ## Single-end short reads
                sample_info = [
                    assembly,
                    "1",
                    fastq_1,
                    "0",
                    coverage_tab,
                    cov_from_assembly,
                ]
            else:
                sample_info = [assembly, "0", "0", "0", coverage_tab, cov_from_assembly]

            ## Create sample mapping dictionary = { sample: [ assembly, single_end, fastq_1, fastq_2, cov_from_assembly ] }
            if sample in samples:
                print_error("Samplesheet contains duplicate samples!", "Line", line)
            samples[sample] = sample_info

    if not samples:
        print_error("No entries to process!", f"Samplesheet: {file_in}")

    ## Write validated samplesheet with appropriate columns
    sample_lines = ""
    for sample, sample_info in samples.items():
        sample_info_line = ",".join(sample_info)
        sample_line = f"{sample},{sample_info_line}\n"
        sample_lines += sample_line
    out_dir = os.path.dirname(file_out)
    make_dir(out_dir)
    header = (
        ",".join(
            [
                "sample",
                "assembly",
                "single_end",
                "fastq_1",
                "fastq_2",
                "coverage_tab",
                "cov_from_assembly",
            ]
        )
        + "\n"
    )
    with open(file_out, "w") as outfh:
        outfh.write(header)
        outfh.write(sample_lines)


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
