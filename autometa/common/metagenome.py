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

from autometa.common import kmers
from autometa.common import coverage
from autometa.common.external import prodigal
from autometa.common.metabin import MetaBin
from autometa.common.utilities import timeit
from autometa.common.utilities import gunzip
from autometa.taxonomy.majority_vote import majority_vote
from autometa.taxonomy.ncbi import NCBI, NCBI_DIR

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
    taxonomy_fpath : str
        </path/to/taxonomy.tsv>
    taxon_method : str, optional
        method to assign taxonomy (the default is 'majority_vote').
        choices=['majority_vote']
    fwd_reads : list, optional
        [</path/to/forward_reads.fastq>, ...]
    rev_reads : list, optional
        [</path/to/reverse_reads.fastq>, ...]
    se_reads : list, optional
        [</path/to/single_end_reads.fastq>, ...]

    Attributes
    ----------
    taxonomy_fname : str
        basename of `taxonomy_fpath`
    taxonomy : pd.DataFrame
        index=contig cols=[taxid] may also contain lineage of taxid
    taxonomy_assigned : bool
        True if `taxonomy_fpath` exists else False
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
    * self.get_kmers()
    * self.assign_taxonomy()
    * self.get_kingdoms()
    * self.write_ranks()

    """

    def __init__(
        self,
        assembly,
        outdir,
        nucl_orfs_fpath,
        prot_orfs_fpath,
        taxonomy_fpath,
        taxon_method="majority_vote",
        fwd_reads=None,
        rev_reads=None,
        se_reads=None,
    ):
        self.assembly = os.path.realpath(assembly)
        self.fwd_reads = fwd_reads
        self.rev_reads = rev_reads
        self.se_reads = se_reads
        self.outdir = outdir
        self.taxon_method = taxon_method
        self.nucl_orfs_fpath = nucl_orfs_fpath
        self.prot_orfs_fpath = prot_orfs_fpath
        self.taxonomy_fpath = taxonomy_fpath
        self.taxonomy_fname = os.path.basename(self.taxonomy_fpath)
        self.taxonomy = (
            pd.read_csv(self.taxonomy_fpath, sep="\t", index_col="contig")
            if self.taxonomy_assigned
            else None
        )

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
        weights = [len(seq) / self.size for seq in self.sequences]
        gc_counts = [SeqUtils.GC(seq) for seq in self.sequences]
        return np.average(a=gc_counts, weights=weights)

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

    @property
    def taxonomy_assigned(self):
        """Retrieve whether taxonomy has been assigned to `assembly`. This will
        check whether `taxonomy_fpath` exists and is non-empty.

        Note: In the future, a checksum comparison should be performed to ensure
        file integrity.

        Returns
        -------
        type
            Description of returned object.

        Raises
        -------
        ExceptionName
            Why the exception is raised.

        """
        # COMBAK: Add checkpointing checksum check here
        if (
            os.path.exists(self.taxonomy_fpath)
            and os.path.getsize(self.taxonomy_fpath) > 0
        ):
            return True
        return False

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
            f"Taxonomy method: {self.taxon_method}\n"
            f"Taxonomy assigned: {self.taxonomy_assigned}\n"
            f"Taxonomy filepath: {self.taxonomy_fpath}\n"
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
        cutoff : int
            Lengths above or equal to `cutoff` that will be retained (the default is 3000).

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
            taxonomy_fpath=self.taxonomy_fpath,
            fwd_reads=self.fwd_reads,
            rev_reads=self.rev_reads,
            taxon_method=self.taxon_method,
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
        OSError
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
            logger.exception(err)
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

    @timeit
    def get_kmers(
        self,
        kmer_size=5,
        multiprocess=True,
        out=None,
        normalized=None,
        force=False,
        nproc=1,
    ):
        """Counts k-mer frequencies using provided `kmer_size`.

        Parameters
        ----------
        kmer_size : int, optional
            length of k-mer to count (the default is 5).
        normalized : str, optional
            Perform Centered-log ratio normalization on counted k-mers and write
            to provided `normalized` file path (the default is None).
        out : str, optional
            Write counted k-mers to `out` (the default is None).
        force : bool, optional
            Overwrite existing k-mers `out` file (the default is False).

        Returns
        -------
        pandas.DataFrame
            pandas DataFrame

        TODO: get_kmers should handle both files and SeqRecords...
        NOTE: above TODO should be handled in kmers.py not here...

        """
        out_specified = out is not None
        out_exists = os.path.exists(out) if out else False
        case1 = out_specified and out_exists and not force
        case2 = out_specified and out_exists and force
        case3 = out_specified and not out_exists
        if case1:
            logger.warning(f"FileExists: {out} force to overwrite. [retrieving]")
            return pd.read_csv(out, sep="\t", index_col="contig")
        normalize_kmers = True if normalized else False
        logger.info(f"Counting {kmer_size}-mers. Normalize: {normalize_kmers}")
        kmers_df = kmers.count(
            assembly=self.assembly,
            kmer_size=kmer_size,
            multiprocess=multiprocess,
            nproc=nproc,
        )
        if case2 or case3:
            kmers_df.to_csv(out, sep="\t", header=True, index=True)
        if normalize_kmers:
            normalized_df = kmers.normalize(kmers_df)
            normalized_df.to_csv(normalized, sep="\t", header=True, index=True)
            return normalized_df
        else:
            return kmers_df

    @timeit
    def get_coverages(self, out, from_spades=True, **kwargs):
        """Retrieve contig coverages using provided `assembly` and `*_reads` or
        if the metagenome assembly was generated from SPAdes, use the k-mer coverages
        provided in each contig's header.

        Parameters
        ----------
        out : str
            </path/to/write/coverages.tsv>
        from_spades : bool
            Description of parameter `from_spades` (the default is True).
        **kwargs : dict
            May contain the following keys: 'sam', 'bam', 'lengths', 'bed'
            Keys should correspond to their respective alignment files.
            'lengths' is a path to a tab-delimited table of contig and its length.

        Returns
        -------
        pd.DataFrame
            index=contig cols=['coverage']

        """
        if from_spades:
            return coverage.from_spades_names(self.seqrecords)
        return coverage.get(
            fasta=self.assembly,
            out=out,
            fwd_reads=self.fwd_reads,
            rev_reads=self.rev_reads,
            se_reads=self.se_reads,
            sam=kwargs.get("sam"),
            bam=kwargs.get("bam"),
            lengths=kwargs.get("lengths"),
            bed=kwargs.get("bed"),
        )

    @timeit
    def assign_taxonomy(self, method, force=False, *args, **kwargs):
        """Assign taxonomy to each sequence in assembly.

        Parameters
        ----------
        force : bool, optional
            overwrite existing voting method's file (the default is False).
        *args : list
            Description of parameter `*args`.
        **kwargs : dict
            May contain the following keys:

            * cpus : int, num. cpus to use
            * ncbi : str, <path/to/ncbi/databases/directory>
            * usepickle : bool, whether to pickle taxonomy-specific files
            * verbose : bool, Add verbosity to stream to terminal
            * blast : str, </path/to/diamond/output/blast.tsv>
            * hits : str, </path/to/diamond/hits.pkl.gz>
            * force : bool, force overwrite existing results files

        Raises
        -------
        NotImplementedError
            Provided `method` not yet implemented.

        """
        logger.debug(f"assigning taxonomy via {method}")
        if not self.orfs_called:
            cpus = kwargs.get("cpus", 0)
            try:
                self.call_orfs(force=force, cpus=cpus)
            except FileExistsError as err:
                logger.warning(err)
        if self.taxonomy_assigned and not force:
            logger.debug(
                f"FileExistsError: {self.taxonomy_fpath}. Use force to overwrite. skipping..."
            )
            return pd.read_csv(self.taxonomy_fpath, sep="\t", index_col="contig")
        if method == "majority_vote":
            self.taxonomy_fpath = majority_vote(
                fasta=self.prot_orfs_fpath,
                ncbi_dir=kwargs.get("ncbi", NCBI_DIR),
                outdir=self.outdir,
                votes_fname=self.taxonomy_fname,
                *args,
                **kwargs,
            )
        else:
            raise NotImplementedError(
                f"method: {method}\nargs:{args}\nkwargs: {kwargs}"
            )
        return pd.read_csv(self.taxonomy_fpath, sep="\t", index_col="contig")

    @timeit
    def get_kingdoms(self, **kwargs):
        """Separate sequences by kingdom using supplied taxon assignment method.

        Parameters
        ----------
        **kwargs : dict, optional
            May contain the following keys:

            * cpus : int, num. cpus to use
            * ncbi : str, <path/to/ncbi/databases/directory>
            * usepickle : bool, whether to pickle taxonomy-specific files
            * verbose : bool, Add verbosity to stream to terminal
            * blast : str, </path/to/diamond/output/blast.tsv>
            * hits : str, </path/to/diamond/hits.pkl.gz>
            * force : bool, force overwrite existing results files


        Returns
        -------
        dict
            {'bacteria':MetaBin, ...}

        Raises
        -------
        KeyError
            Why the exception is raised.

        """
        if not self.taxonomy_assigned:
            logger.info("Assigning taxonomy. This may take a while...")
            self.taxonomy = self.assign_taxonomy(method=self.taxon_method, **kwargs)
        if self.taxonomy.shape[1] <= 2:
            # taxonomy_fp should only contain contig and taxid columns from voting method
            ncbi = NCBI(kwargs.get("ncbi", NCBI_DIR))
            dff = ncbi.get_lineage_dataframe(self.taxonomy["taxid"].unique().tolist())
            self.taxonomy = pd.merge(
                left=self.taxonomy,
                right=dff,
                how="left",
                left_on="taxid",
                right_index=True,
            )
            self.taxonomy.to_csv(self.taxonomy_fpath, sep="\t", index=True, header=True)
            # COMBAK: Add checkpointing checksum check here
            logger.debug(f"Added canonical rank names to {self.taxonomy_fpath}")
        if "superkingdom" not in self.taxonomy.columns:
            raise KeyError(
                f"superkingdom is not in taxonomy columns {self.taxonomy.columns}"
            )
        kingdoms = dict(list(self.taxonomy.groupby("superkingdom")))
        bins = {}
        for kingdom, df in kingdoms.items():
            bins.update({kingdom: MetaBin(self.assembly, df.index.tolist())})
        return bins

    def write_ranks(self, rank="superkingdom"):
        """Write fastas split by `rank`.

        Parameters
        ----------
        rank : str
            `rank` from canonical ranks (the default is 'superkingdom').

        Returns
        -------
        list
            [rank_name_fpath, ...]

        Raises
        -------
        ValueError
            `rank` not in canonical ranks

        """
        if rank not in NCBI.CANONICAL_RANKS:
            raise ValueError(f"rank: {rank} not in {NCBI.CANONICAL_RANKS}")
        fpaths = []
        for rank_name, dff in self.taxonomy.groupby(rank):
            records = [r for r in self.seqrecords if r.id in dff.index]
            rank_name = rank_name.replace(" ", "_")
            rank_name_fname = ".".join([rank_name.title(), "fna"])
            rank_name_fpath = os.path.join(self.outdir, rank_name_fname)
            if not records:
                logger.warning(f"No records to write to {rank_name_fpath}")
            else:
                n_written = SeqIO.write(records, rank_name_fpath, "fasta")
                logger.debug(f"Wrote {n_written} records to {rank_name_fpath}")
                fpaths.append(rank_name_fpath)
        return fpaths


def taxonomy_entrypoint():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )

    parser = argparse.ArgumentParser(
        description="""
    This script handles filtering by taxonomy and can calculate various metagenome statistics.
    """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "assembly", help="Path to metagenome assembly (nucleotide fasta)."
    )
    parser.add_argument("outdir", help="Output directory to store annotations.")
    parser.add_argument(
        "nucls",
        help="Path to nucleotide ORFs corresponding to `assembly.` "
        "(Will write to path if ORFs do not exist).",
        type=str,
    )
    parser.add_argument(
        "prots",
        help="Path to amino acid ORFs corresponding to `assembly.` "
        "(Will write to path if ORFs do not exist).",
        type=str,
    )

    parser.add_argument(
        "--taxon-method", default="majority_vote", choices=["majority_vote"]
    )
    parser.add_argument(
        "--taxon-fname",
        help="Filename to assign voted taxonomy annotation.",
        default="taxonomy_vote.tsv",
    )
    parser.add_argument(
        "--ncbi", help="Path to NCBI databases directory.", default=NCBI_DIR
    )
    parser.add_argument(
        "--pickle",
        help="Whether to serialize taxonomy-specific files",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--blast",
        help="Path to diamond blast results.tsv. (outfmt=6) "
        "(Will write to path if it does not exist).",
        default=None,
    )
    parser.add_argument(
        "--hits",
        help="Path to diamond blast hits.pkl.gz. "
        "(Will write to path if it does not exist).",
        default=None,
    )
    # Eventually will need to create a subparser for the taxon assignment methods
    # to include help information and required parameters according to method.
    parser.add_argument(
        "--cpus",
        help="Number of cpus to use when performing "
        "annotations (Default will use all possible if `--parallel` is supplied).",
        default=0,
        type=int,
    )
    parser.add_argument(
        "--parallel",
        help="Use GNU parallel when performing annotations.",
        action="store_true",
        default=False,
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

    taxonomy_fpath = os.path.join(args.outdir, args.taxon_fname)
    mg = Metagenome(
        assembly=args.assembly,
        outdir=args.outdir,
        prot_orfs_fpath=args.prots,
        nucl_orfs_fpath=args.nucls,
        taxonomy_fpath=taxonomy_fpath,
        taxon_method=args.taxon_method,
    )

    if args.stats:
        mg.describe(autometa_details=False)

    try:
        mg.call_orfs(
            force=args.force, cpus=args.cpus, parallel=args.parallel,
        )
    except FileExistsError:
        logger.warning(f"{mg.prots_out} already exists. Skipping...")

    logger.info(f"autometa-taxonomy method={args.taxon_method}")
    mg.get_kingdoms(
        cpus=args.cpus,
        ncbi=args.ncbi,
        usepickle=args.pickle,
        verbose=args.verbose,
        blast=args.blast,
        hits=args.hits,
        force=args.force,
    )


def length_filter_entrypoint():
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
        taxonomy_fpath="",
        force=args.force,
    )

    filtered_mg = raw_mg.length_filter(
        out=args.out, cutoff=args.cutoff, force=args.force
    )
    if args.stats:
        filtered_mg.describe(autometa_details=False)


def main():
    import sys

    help_info = """
    autometa-metagenome script to filter assembly by length cutoff or assign taxonomy.

    For taxonomy-based filtering usage information:
    autometa.common.metagenome taxonomy --help

    For length-cutoff filtering usage information:
    autometa.common.metagenome length-filter --help

    """
    # main() is available here if the user wants to run autometa without installing console-scripts.
    if sys.argv and "taxonomy" in sys.argv:
        sys.exit(taxonomy_entrypoint())

    if sys.argv and "length-filter" in sys.argv:
        sys.exit(length_filter_entrypoint())
    else:
        print(help_info)


if __name__ == "__main__":
    main()
