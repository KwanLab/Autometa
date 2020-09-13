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


import decimal
import logging
import numbers
import os

import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import SeqUtils
from functools import lru_cache

from autometa.common.external import prodigal
from autometa.common.utilities import timeit
from autometa.common.utilities import gunzip
from autometa.common.exceptions import ExternalToolError

# TODO: Should place all imports of Database paths to a config module so they
# all exist in one place

logger = logging.getLogger(__name__)


class Metagenome:
    """Autometa Metagenome Class.

    Parameters
    ----------
    assembly : str
        </path/to/metagenome/assembly.fasta>
    outdir : str
        </path/to/output/directory> (the default is None)
    nucl_orfs_fpath : str
        </path/to/assembly.orfs.fna>
    prot_orfs_fpath : str
        </path/to/assembly.orfs.faa>
    fwd_reads : list, optional
        [</path/to/forward_reads.fastq>, ...]
    rev_reads : list, optional
        [</path/to/reverse_reads.fastq>, ...]
    se_reads : list, optional
        [</path/to/single_end_reads.fastq>, ...]

    Attributes
    ----------
    orfs_called : bool
        True if both `nucl_orfs_fpath` and `prot_orfs_fpath` exist else False
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
    nucls : list
        SeqRecords all correspond to retrieved nucleotide ORFs. [SeqRecord, ...].
    prots : list
        SeqRecords all correspond to retrieved amino-acid ORFs. [SeqRecord, ...].

    Methods
    ----------
    * self.fragmentation_metric()
    * self.describe()
    * self.length_filter()
    * self.call_orfs()
    * self.orfs()

    """

    def __init__(
        self,
        assembly,
        outdir,
        nucl_orfs_fpath,
        prot_orfs_fpath,
        fwd_reads=None,
        rev_reads=None,
        se_reads=None,
    ):
        self.assembly = os.path.realpath(assembly)
        self.fwd_reads = fwd_reads
        self.rev_reads = rev_reads
        self.se_reads = se_reads
        self.outdir = outdir
        self.nucl_orfs_fpath = nucl_orfs_fpath
        self.prot_orfs_fpath = prot_orfs_fpath

    def __repr__(self):
        return str(self)

    def __str__(self):
        return self.assembly

    @property
    @lru_cache(maxsize=None)
    def sequences(self):
        """Retrieve the sequences from provided `assembly`.

        Returns
        -------
        list
            [seq, seq, ...]

        """
        with open(self.assembly) as fh:
            return [seq for title, seq in SimpleFastaParser(fh)]

    @property
    @lru_cache(maxsize=None)
    def seqrecords(self):
        """Retrieve SeqRecord objects from provided `assembly`.

        Returns
        -------
        list
            [SeqRecord, SeqRecord, ...]

        """
        return [seq for seq in SeqIO.parse(self.assembly, "fasta")]

    @property
    def nseqs(self):
        """Retrieve the number of sequences in provided `assembly`.

        Returns
        -------
        int
            Number of sequences parsed from `assembly`

        """
        return len(self.sequences)

    @property
    @lru_cache(maxsize=None)
    def length_weighted_gc(self):
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
    def size(self):
        """Retrieve the summation of sizes for each contig in the provided `assembly`.

        Returns
        -------
        int
            Total summation of contig sizes in `assembly`

        """
        return sum(len(seq) for seq in self.sequences)

    @property
    def largest_seq(self):
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

    @property
    def orfs_called(self):
        """Retrieve whether `prot_orfs_fpath` and `nucl_orfs_fpath` have been called.

        Note: This will check whether the aforementioned paths exist and are not empty.
        In the future, a checksum comparison will be performed to ensure file integrity.

        Returns
        -------
        bool
            Description of returned object.

        """
        # COMBAK: Add checkpointing checksum check here
        for fp in [self.prot_orfs_fpath, self.nucl_orfs_fpath]:
            if not os.path.exists(fp):
                return False
            elif not os.path.getsize(fp) > 0:
                return False
        return True

    @property
    @lru_cache(maxsize=None)
    def nucls(self):
        """Retrieve `assembly` nucleotide ORFs.

        Returns
        -------
        list
            [SeqRecord, SeqRecord, ...]

        """
        return self.orfs(orf_type="nucl")

    @property
    @lru_cache(maxsize=None)
    def prots(self):
        """Retrieve `assembly` amino-acid ORFs.

        Returns
        -------
        list
            [SeqRecord, SeqRecord, ...]

        """
        return self.orfs(orf_type="prot")

    def fragmentation_metric(self, quality_measure=0.50):
        """Describes the quality of assembled genomes that are fragmented in
        contigs of different length.

        For more information see:
            http://www.metagenomics.wiki/pdf/definition/assembly/n50

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
        target_size = self.size * quality_measure
        lengths = []
        for length in sorted([len(seq) for seq in self.sequences], reverse=True):
            lengths.append(length)
            if sum(lengths) > target_size:
                return length

    def describe(self, autometa_details=True):
        """Print `assembly` details.

        Parameters
        ----------
        autometa_details : bool
            Also log Autometa specific information to the terminal (Default is True).

        Returns
        -------
        NoneType

        """
        print(
            "Metagenome Details\n"
            "________________________\n"
            f"Assembly: {self.assembly}\n"
            f"Num. Sequences: {self.nseqs:,}\n"
            f"Size: {self.size:,} bp\n"
            f"N50: {self.fragmentation_metric():,} bp\n"
            f"N10: {self.fragmentation_metric(.1):,} bp\n"
            f"N90: {self.fragmentation_metric(.9):,} bp\n"
            f"Length Weighted Avg. GC content: {self.length_weighted_gc:4.2f}%\n"
            f"Largest sequence: {self.largest_seq}\n"
            "________________________\n"
        )
        if not autometa_details:
            return
        print(
            "Autometa Details\n"
            "________________________\n"
            f"Outdir: {self.outdir}\n"
            f"ORFs called: {self.orfs_called}\n"
            f"Prots filepath: {self.prot_orfs_fpath}\n"
            f"Nucl filepath: {self.nucl_orfs_fpath}\n"
        )

    @timeit
    def length_filter(self, out, cutoff=3000, force=False):
        """Filters sequences by length with provided cutoff.

        Note: A WARNING will be emitted if the length filter is applied *after*
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
        outdir = os.path.dirname(out)
        gunzipped_fname = os.path.basename(self.assembly.rstrip(".gz"))
        gunzipped_fpath = os.path.join(outdir, gunzipped_fname)
        if self.assembly.endswith(".gz"):
            if not os.path.exists(gunzipped_fpath):
                gunzip(self.assembly, gunzipped_fpath)
            self.assembly = gunzipped_fpath
        logger.info(f"Getting contigs greater than or equal to {cutoff:,} bp")
        records = [seq for seq in self.seqrecords if len(seq) >= cutoff]
        if self.orfs_called:
            msg = (
                f"{self.nucl_orfs_fpath} and {self.prot_orfs_fpath} have already been called!"
                "Call orfs again to retrieve only ORFs corresponding to filtered assembly"
            )
            logger.warning(msg)
        SeqIO.write(records, out, "fasta")
        return Metagenome(
            assembly=out,
            outdir=self.outdir,
            nucl_orfs_fpath=self.nucl_orfs_fpath,
            prot_orfs_fpath=self.prot_orfs_fpath,
            fwd_reads=self.fwd_reads,
            rev_reads=self.rev_reads,
        )

    def call_orfs(self, force=False, cpus=0, parallel=True):
        """Calls ORFs on Metagenome assembly.

        (Wrapper using external executable: prodigal).

        Parameters
        ----------
        force : bool
            force overwrite of existing ORFs files (the default is False).
        cpus : int
            Description of parameter `cpus` (the default is 0).
        parallel : bool
            Will parallelize prodigal using GNU parallel (the default is True).

        Returns
        ----------
        2-tuple
            (</path/to/nucls.orfs.fna>, </path/to/prots.orfs.faa>)
        Raises
        -------
        TypeError
            `force`,`parallel` or `cpus` type was incorrectly supplied.
        ChildProcessError
            ORF calling failed.
        """
        for arg, argname in zip([force, parallel], ["force", "parallel"]):
            if not isinstance(arg, bool) and isinstance(arg, numbers.Number):
                raise TypeError(f"{argname} must be a boolean!")
        if not isinstance(cpus, int) or isinstance(cpus, bool):
            raise TypeError(f"cpus:({cpus}) must be an integer!")

        # COMBAK: Add checkpointing checksum check here
        try:
            nucls_fp, prots_fp = prodigal.run(
                assembly=self.assembly,
                nucls_out=self.nucl_orfs_fpath,
                prots_out=self.prot_orfs_fpath,
                force=force,
                cpus=cpus,
                parallel=parallel,
            )
        except FileExistsError as err:
            return self.nucl_orfs_fpath, self.prot_orfs_fpath
        except ChildProcessError as err:
            raise err
        return nucls_fp, prots_fp

    def orfs(self, orf_type="prot", cpus=0):
        """Retrieves ORFs after being called from self.call_orfs.

        Parameters
        ----------
        orf_type : str
            format of ORFs to retrieve choices=['nucl','prot'] either nucleotide
            or amino acids (the default is 'prot').

        Returns
        -------
        list
            [SeqRecord, ...]

        Raises
        -------
        ValueError
            Invalid `orf_type`. Choices=['prot','nucl']

        """
        if not self.orfs_called:
            self.call_orfs(cpus=cpus)
        if orf_type not in {"prot", "nucl"}:
            raise ValueError('orf_type must be "prot" or "nucl"!')
        orfs_fpath = (
            self.prot_orfs_fpath if orf_type == "prot" else self.nucl_orfs_fpath
        )
        return [orf for orf in SeqIO.parse(orfs_fpath, "fasta")]


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
        "assembly", help="Path to metagenome assembly (nucleotide fasta)."
    )
    parser.add_argument(
        "out", help="Path to output length-filtered assembly fasta file.",
    )
    parser.add_argument(
        "--cutoff",
        help="Cutoff to apply to length filter",
        default=3000,
        type=int,
        metavar="<int>",
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
    parser.add_argument(
        "--stats",
        help="Print various metagenome assembly statistics.",
        action="store_true",
        default=False,
    )

    args = parser.parse_args()
    dirpath = os.path.dirname(os.path.realpath(args.assembly))
    raw_mg = Metagenome(
        assembly=args.assembly,
        outdir=dirpath,
        prot_orfs_fpath="",
        nucl_orfs_fpath="",
        force=args.force,
    )

    filtered_mg = raw_mg.length_filter(
        out=args.out, cutoff=args.cutoff, force=args.force
    )
    if args.stats:
        filtered_mg.describe(autometa_details=False)


if __name__ == "__main__":
    main()
