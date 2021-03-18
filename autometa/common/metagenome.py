#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
COPYRIGHT
Copyright 2020 Ian J. Miller, Evan R. Rees, Kyle Wolf, Siddharth Uppal,
Shaurya Chanana, Izaak Miller, Jason C. Kwan

This file is part of Autometa.

Autometa is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Autometa is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with Autometa. If not, see <http://www.gnu.org/licenses/>.
COPYRIGHT

Script containing Metagenome class for general handling of metagenome assembly
"""


import gzip
import logging
import numbers
import os

import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import SeqUtils
from functools import lru_cache

from autometa.common.utilities import timeit

logger = logging.getLogger(__name__)


class Metagenome:
    """Autometa Metagenome Class.

    Parameters
    ----------
    assembly : str
        </path/to/metagenome/assembly.fasta>

    Attributes
    ----------
    sequences : list
        [seq,...]
    seqrecords : list
        [SeqRecord,...]
    nseqs : int
        Number of sequences in assembly.
    length_weighted_gc : float
        Length weighted average GC% of assembly.
    size : int
        Total assembly size in bp.
    largest_seq : str
        id of longest sequence in assembly

    Methods
    ----------
    * self.fragmentation_metric()
    * self.describe()
    * self.length_filter()

    """

    def __init__(
        self,
        assembly,
    ):
        self.assembly = os.path.realpath(assembly)

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return self.assembly

    @property
    @lru_cache(maxsize=None)
    def sequences(self) -> list:
        """Retrieve the sequences from provided `assembly`.

        Returns
        -------
        list
            [seq, seq, ...]

        """
        fh = (
            gzip.open(self.assembly, "rt")
            if self.assembly.endswith(".gz")
            else open(self.assembly)
        )
        sequences = [seq for __, seq in SimpleFastaParser(fh)]
        fh.close()
        return sequences

    @property
    @lru_cache(maxsize=None)
    def seqrecords(self) -> list:
        """Retrieve SeqRecord objects from provided `assembly`.

        Returns
        -------
        list
            [SeqRecord, SeqRecord, ...]

        """
        fh = (
            gzip.open(self.assembly, "rt")
            if self.assembly.endswith(".gz")
            else open(self.assembly)
        )
        records = [seq for seq in SeqIO.parse(fh, "fasta")]
        fh.close()
        return records

    @property
    def nseqs(self) -> int:
        """Retrieve the number of sequences in provided `assembly`.

        Returns
        -------
        int
            Number of sequences parsed from `assembly`

        """
        return len(self.sequences)

    @property
    @lru_cache(maxsize=None)
    def length_weighted_gc(self) -> float:
        """Retrieve the length weighted average GC percentage of provided `assembly`.

        Returns
        -------
        float
            GC percentage weighted by contig length.

        """
        weights = []
        add_weight = weights.append
        gc_contents = []
        add_gc_contents = gc_contents.append
        for seq in self.sequences:
            weight = len(seq) / self.size
            add_weight(weight)
            gc_content = SeqUtils.GC(seq)
            add_gc_contents(gc_content)
        return np.average(a=gc_contents, weights=weights)

    @property
    def size(self) -> int:
        """Retrieve the summation of sizes for each contig in the provided `assembly`.

        Returns
        -------
        int
            Total summation of contig sizes in `assembly`

        """
        return sum(len(seq) for seq in self.sequences)

    @property
    def largest_seq(self) -> str:
        """Retrieve the name of the largest sequence in the provided `assembly`.

        Returns
        -------
        str
            record ID of the largest sequence in `assembly`.

        """
        max = float("-inf")
        largest = None
        for rec in self.seqrecords:
            if len(rec) > max:
                largest = rec
                max = len(rec)
        return largest.id

    def fragmentation_metric(self, quality_measure: float = 0.50) -> int:
        """Describes the quality of assembled genomes that are fragmented in
        contigs of different length.

        Note
        ----
        For more information see this metagenomics `wiki <http://www.metagenomics.wiki/pdf/definition/assembly/n50>`_ from Matthias Scholz

        Parameters
        ----------
        quality_measure : 0 < float < 1
            Description of parameter `quality_measure` (the default is .50).
            I.e. default measure is N50, but could use .1 for N10 or .9 for N90

        Returns
        -------
        int
            Minimum contig length to cover `quality_measure` of genome (i.e. length
            weighted median)

        """
        if not 0 < quality_measure < 1:
            raise ValueError(
                f"quality measure must be b/w 0 and 1! given: {quality_measure}"
            )
        target_size = self.size * quality_measure
        total_length = 0
        # Sorts lengths from longest to shortest
        for length in sorted([len(seq) for seq in self.sequences], reverse=True):
            total_length += length
            if total_length >= target_size:
                return length

    @timeit
    def describe(self) -> pd.DataFrame:
        """Return dataframe of details.

        Columns
        -------

            # assembly : Assembly input into Metagenome(...) [index column]
            # nseqs : Number of sequences in assembly
            # size : Size or total sum of all sequence lengths
            # N50 :
            # N10 :
            # N90 :
            # length_weighted_gc_content : Length weighted average GC content
            # largest_seq : Largest sequence in assembly

        Returns
        -------
        pd.DataFrame

        """
        return pd.DataFrame(
            [
                {
                    "assembly": self.assembly,
                    "nseqs": self.nseqs,
                    "size (bp)": self.size,
                    "N50 (bp)": self.fragmentation_metric(0.5),
                    "N10 (bp)": self.fragmentation_metric(0.1),
                    "N90 (bp)": self.fragmentation_metric(0.9),
                    "length_weighted_gc_content (%)": self.length_weighted_gc,
                    "largest_seq": self.largest_seq,
                }
            ]
        ).set_index("assembly")

    @timeit
    def length_filter(self, out: str, cutoff: int = 3000, force: bool = False):
        """Filters sequences by length with provided cutoff.

        Note
        ----
        A WARNING will be emitted if the length filter is applied *after*
        the ORFs provided for the Metagenome are already called prompting the
        user to perform orf calling again to correspond to length filtered
        contigs.

        Parameters
        ----------
        out : str
            Path to write length filtered output fasta file.
        cutoff : int, optional
            Lengths above or equal to `cutoff` that will be retained (the default is 3000).
        force : bool, optional
            Overwrite existing `out` file (the default is False).

        Returns
        -------
        Metagenome
            autometa Metagenome object with only assembly sequences above the cutoff threshold.

        Raises
        -------
        TypeError
            cutoff value must be a float or integer
        ValueError
            cutoff value must be a positive real number
        FileExistsError
            filepath consisting of sequences that passed filter already exists
        """

        if not isinstance(cutoff, numbers.Number) or isinstance(cutoff, bool):
            # https://stackoverflow.com/a/4187220/13118765
            raise TypeError(f"cutoff: {cutoff} must be a float or int")
        if cutoff <= 0:
            raise ValueError(f"cutoff: {cutoff} must be a positive real number")
        if os.path.exists(out) and not force:
            raise FileExistsError(out)
        logger.info(f"Getting contigs greater than or equal to {cutoff:,} bp")
        records = [record for record in self.seqrecords if len(record.seq) >= cutoff]
        n_written = SeqIO.write(records, out, "fasta")
        if not n_written:
            logger.warning(
                f"No contigs pass {cutoff:,}. Returning original Metagenome."
            )
            return self
        else:
            return Metagenome(assembly=out)

    @timeit
    def gc_content(self) -> pd.DataFrame:
        """Retrieves GC content from sequences in assembly

        Returns
        -------
        pd.DataFrame
            index="contig", columns=["gc_content","length"]
        """
        return pd.DataFrame(
            [
                {
                    "contig": record.id,
                    "gc_content": SeqUtils.GC(record.seq),
                    "length": len(record.seq),
                }
                for record in self.seqrecords
            ]
        ).set_index("contig")


def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )
    parser = argparse.ArgumentParser(
        description="This script handles filtering by length and can calculate various metagenome statistics.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--assembly",
        help="Path to metagenome assembly (nucleotide fasta).",
        metavar="filepath",
        required=True,
    )
    parser.add_argument(
        "--output-fasta",
        help="Path to output length-filtered assembly fasta file.",
        metavar="filepath",
        required=True,
    )
    parser.add_argument(
        "--output-stats",
        help="Path to output length-filtered assembly fasta file.",
        metavar="filepath",
        required=False,
    )
    parser.add_argument(
        "--output-gc-content",
        help="Path to output assembly contigs' GC content and length.",
        metavar="filepath",
        required=False,
    )
    parser.add_argument(
        "--cutoff",
        help="Cutoff to apply to length filter",
        default=3000,
        type=int,
        metavar="int",
    )
    parser.add_argument(
        "--force", help="Overwrite existing files", action="store_true", default=False
    )
    parser.add_argument(
        "--verbose",
        help="Log more information to terminal.",
        action="store_true",
        default=False,
    )

    args = parser.parse_args()
    raw_mg = Metagenome(assembly=args.assembly)

    filtered_mg = raw_mg.length_filter(
        out=args.output_fasta, cutoff=args.cutoff, force=args.force
    )
    if args.output_stats:
        stats_df = filtered_mg.describe()
        stats_df.to_csv(args.output_stats, sep="\t", index=True, header=True)
    if args.output_gc_content:
        gc_df = filtered_mg.gc_content()
        gc_df.to_csv(args.output_gc_content, sep="\t", index=True, header=True)


if __name__ == "__main__":
    main()
