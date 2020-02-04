#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script containing Metagenome class for general handling of metagenome assembly
"""


import logging
import os

import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio import SeqUtils

from autometa.common import kmers
from autometa.common import coverage
from autometa.common.external import prodigal
from autometa.common.mag import MAG
from autometa.common.utilities import timeit
from autometa.common.utilities import gunzip
from autometa.taxonomy.majority_vote import majority_vote
from autometa.taxonomy.ncbi import NCBI,NCBI_DIR

# TODO: Should place all imports of Database paths to a config module so they
# all exist in one place

logger = logging.getLogger(__name__)


class Metagenome:
    """Autometa Metagenome Class.

    Parameters
    ----------
    assembly : str
        </path/to/assembled/metagenome.fasta>
    outdir : type
        </path/to/output/directory> (the default is None)
    taxon_method : str
        method to assign taxonomy (the default is 'majority_vote').
        choices=['majority_vote']

    Attributes
    ----------
    nucl_orfs_fpath : str
        Description of attribute `nucl_orfs_fpath`.
    prot_orfs_fpath : str
        Description of attribute `prot_orfs_fpath`.
    taxonomy_fname : str
        Description of attribute `taxonomy_fname`.
    taxonomy_fpath : str
        Description of attribute `taxonomy_fpath`.
    taxonomy : pd.DataFrame
        index=contig cols=[taxid] may also contain lineage of taxid
    taxonomy_assigned : bool
        `taxonomy_fpath` exists
    orfs_called : bool
        `nucl_orfs_fpath` and `prot_orfs_fpath` exists
    sequences : list
        [SeqRecord,...]
    nseqs : int
        Number of sequences in assembly.
    mean_gc : float
        Mean GC% of assembly.
    size : int
        Total assembly size.
    largest_seq : str
        id of longest sequence in assembly
    nucls : list
        SeqRecords all correspond to retrieved nucleotide ORFs. [SeqRecord, ...].
    prots : list
        SeqRecords all correspond to retrieved amino-acid ORFs. [SeqRecord, ...].

    Methods
    ----------
    - self.fragmentation_metric()
    - self.describe()
    - self.length_filter()
    - self.call_orfs()
    - self.orfs()
    - self.get_kmers()
    - self.assign_taxonomy()
    - self.get_kingdoms()
    - self.write_ranks()
    """
    def __init__(self, assembly, outdir, nucl_orfs_fpath, prot_orfs_fpath,
        taxonomy_fpath, fwd_reads=None, rev_reads=None, taxon_method='majority_vote'):
        self.assembly = os.path.realpath(assembly)
        self.fwd_reads = fwd_reads
        self.rev_reads = rev_reads
        self.outdir = outdir
        self.taxon_method = taxon_method
        self.nucl_orfs_fpath = nucl_orfs_fpath
        self.prot_orfs_fpath = prot_orfs_fpath
        self.taxonomy_fpath = taxonomy_fpath
        self.taxonomy_fname = os.path.basename(self.taxonomy_fpath)
        self.taxonomy = pd.read_csv(self.taxonomy_fpath,sep='\t',index_col='contig') if self.taxonomy_assigned else None

    def __repr__(self):
        return str(self)

    def __str__(self):
        return self.assembly

    @property
    def sequences(self):
        return [seq for seq in SeqIO.parse(self.assembly, 'fasta')]

    @property
    def nseqs(self):
        return len(self.sequences)

    @property
    def mean_gc(self):
        return np.mean([SeqUtils.GC(record.seq) for record in self.sequences])

    @property
    def size(self):
        return sum(len(record) for record in self.sequences)

    @property
    def largest_seq(self):
        max = float('-inf')
        largest = None
        for rec in self.sequences:
            if len(rec) > max:
                largest = rec
                max = len(rec)
        return largest.id

    @property
    def orfs_called(self):
        for fp in [self.prot_orfs_fpath, self.nucl_orfs_fpath]:
            if not os.path.exists(fp):
                return False
            elif not os.stat(fp).st_size > 0:
                return False
        return True

    @property
    def nucls(self):
        return self.orfs(orf_type='nucl')

    @property
    def prots(self):
        return self.orfs(orf_type='prot')

    @property
    def taxonomy_assigned(self):
        if os.path.exists(self.taxonomy_fpath) and os.stat(self.taxonomy_fpath).st_size > 0:
            return True
        return False

    def fragmentation_metric(self, quality_measure=.50):
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
        for length in sorted([len(r) for r in self.sequences], reverse=True):
            lengths.append(length)
            if sum(lengths) > target_size:
                return length

    def describe(self):
        print(f"""
Metagenome Details
________________________
Assembly: {self.assembly}
Num. Sequences: {self.nseqs:,}
Size: {self.size:,} bp
N50: {self.fragmentation_metric():,} bp
N10: {self.fragmentation_metric(.1):,} bp
N90: {self.fragmentation_metric(.9):,} bp
Mean GC: {self.mean_gc:4.2f}%
Largest sequence: {self.largest_seq}
________________________
Autometa Details
________________________
Outdir: {self.outdir}
ORFs called: {self.orfs_called}
Prots filepath: {self.prot_orfs_fpath}
Nucl filepath: {self.nucl_orfs_fpath}
Taxonomy method: {self.taxon_method}
Taxonomy assigned: {self.taxonomy_assigned}
Taxonomy filepath: {self.taxonomy_fpath}
""")

    @timeit
    def length_filter(self, out, cutoff=3000):
        """Filters sequences by length with provided cutoff.

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
        try:
            cutoff = float(cutoff)
        except Exception as err:
            pass
        if not type(cutoff) in [int, float]:
            raise TypeError(f'cutoff: {cutoff} must be a float or int')
        if cutoff <= 0:
            raise ValueError(f'cutoff: {cutoff} must be a positive real number')
        if os.path.exists(out):
            raise FileExistsError(out)
        outdir = os.path.dirname(out)
        gunzipped_fname = os.path.basename(self.assembly.rstrip('.gz'))
        gunzipped_fpath = os.path.join(outdir, gunzipped_fname)
        if self.assembly.endswith('.gz'):
            if not os.path.exists(gunzipped_fpath):
                gunzip(self.assembly, gunzipped_fpath)
            self.assembly = gunzipped_fpath
        records = [seq for seq in self.sequences if len(seq) >= cutoff]
        SeqIO.write(records, out, 'fasta')
        return Metagenome(
            assembly=out,
            outdir=self.outdir,
            nucl_orfs_fpath=self.nucl_orfs_fpath,
            prot_orfs_fpath=self.prot_orfs_fpath,
            taxonomy_fpath=self.taxonomy_fpath,
            fwd_reads=self.fwd_reads,
            rev_reads=self.rev_reads,
            taxon_method=self.taxon_method)

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
        if type(force) is not bool:
            raise TypeError(f'force:({force}) must be a boolean. I.e. True|False')
        if type(parallel) is not bool:
            raise TypeError(f'parallel:({parallel}) must be a boolean. I.e. True|False')
        if type(cpus) is not int:
            raise TypeError(f'cpus:({cpus}) must be an integer')
        try:
            nucls_fp, prots_fp = prodigal.run(
                assembly=self.assembly,
                nucls_out=self.nucl_orfs_fpath,
                prots_out=self.prot_orfs_fpath,
                force=force,
                cpus=cpus,
                parallel=parallel,
            )
        except OSError as err:
            logger.exception(err)
        except FileExistsError as err:
            return self.nucl_orfs_fpath, self.prot_orfs_fpath
        return nucls_fp, prots_fp

    def orfs(self, orf_type='prot', cpus=0):
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
        if orf_type not in {'prot','nucl'}:
            raise ValueError('orf_type must be \'prot\' or \'nucl\'!')
        orfs_fpath = self.prot_orfs_fpath if orf_type == 'prot' else self.nucl_orfs_fpath
        return [orf for orf in SeqIO.parse(orfs_fpath, 'fasta')]

    @timeit
    def get_kmers(self, kmer_size=5, multiprocess=True, out=None, normalized=None,
        force=False, nproc=1):
        """Counts k-mer frequencies using provided `kmer_size`.

        Parameters
        ----------
        kmer_size : int
            length of k-mer to count (the default is 5).
        normalized : str
            Perform Centered-log ratio normalization on counted k-mers (the default is None).
            and write to provided `normalized` path.
        out : str
            Write counted k-mers to `out` (the default is None).
        force : bool
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
            logger.warning(f'FileExistsError: {out} force to overwrite. [retrieving]')
            return pd.read_csv(out, sep='\t', index_col='contig')
        normalize_kmers = True if normalized else False
        logger.info(f'Counting {kmer_size}-mers. Normalize: {normalize_kmers}')
        kmers_df = kmers.count(
            assembly=self.assembly,
            kmer_size=kmer_size,
            multiprocess=multiprocess,
            nproc=nproc,
        )
        if case2 or case3:
            kmers_df.to_csv(out, sep='\t', header=True, index=True)
        if normalize_kmers:
            normalized_df = kmers.normalize(kmers_df)
            normalized_df.to_csv(normalized, sep='\t', header=True, index=True)
            return normalized_df
        return kmers_df

    @timeit
    def get_coverages(self, out, from_spades=True, **kwargs):
        if from_spades:
            return coverage.from_spades_names(self.sequences)
        return coverage.get(
            fasta=self.assembly,
            out=out,
            fwd_reads=self.fwd_reads,
            rev_reads=self.rev_reads,
            sam=kwargs.get('sam'),
            bam=kwargs.get('bam'),
            lengths=kwargs.get('lengths'),
            bed=kwargs.get('bed'),
        )

    @timeit
    def assign_taxonomy(self, method, force=False, *args, **kwargs):
        """Assign taxonomy to each sequence in assembly.

        Parameters
        ----------
        force : bool
            overwrite existing voting method's file (the default is False).
        *args : type
            Description of parameter `*args`.
        **kwargs : type
            Description of parameter `**kwargs`.

        Raises
        -------
        ExceptionName
            Why the exception is raised.

        """
        if not self.orfs_called:
            cpus = kwargs.get('cpus',0)
            try:
                self.call_orfs(force=force, cpus=cpus)
            except FileExistsError as err:
                logger.warning(err)
        if self.taxonomy_assigned and not force:
            logger.debug(f'FileExistsError: {self.taxonomy_fpath}. Use force to overwrite. skipping...')
            return pd.read_csv(self.taxonomy_fpath, sep='\t', index_col='contig')
        if method == 'majority_vote':
            self.taxonomy_fpath = majority_vote(
                fasta=self.prot_orfs_fpath,
                ncbi_dir=kwargs.get('ncbi',NCBI_DIR),
                outdir=self.outdir,
                votes_fname=self.taxonomy_fname,
                *args,
                **kwargs)
        else:
            raise NotImplementedError(f'method: {method}\nargs:{args}\nkwargs: {kwargs}')
        return pd.read_csv(self.taxonomy_fpath,sep='\t',index_col='contig')

    @timeit
    def get_kingdoms(self, **kwargs):
        """Separate sequences by kingdom using supplied taxon assignment method.

        Parameters
        ----------
        **kwargs : dict
            Optional additional keyword arguments to supply to assign_taxonomy
            and supply a separate NCBI_DIR

        Returns
        -------
        type
            Description of returned object.

        Raises
        -------
        KeyError
            Why the exception is raised.

        """
        if not self.taxonomy_assigned:
            self.taxonomy = self.assign_taxonomy(method=self.taxon_method, **kwargs)
        if self.taxonomy.shape[1] <= 2:
            # taxonomy_fp should only contain contig and taxid columns from voting method
            ncbi = NCBI(kwargs.get('ncbi',NCBI_DIR))
            dff = ncbi.get_lineage_dataframe(self.taxonomy['taxid'].unique().tolist())
            self.taxonomy = pd.merge(
                left=self.taxonomy,
                right=dff,
                how='left',
                left_on='taxid',
                right_index=True)
            self.taxonomy.to_csv(self.taxonomy_fpath,sep='\t',index=True,header=True)
            logger.debug(f'Added canonical rank names to {self.taxonomy_fpath}')
        if 'superkingdom' not in self.taxonomy.columns:
            raise KeyError(f'superkingdom is not in taxonomy columns {self.taxonomy.columns}')
        kingdoms = dict(list(self.taxonomy.groupby('superkingdom')))
        bins = {}
        for kingdom, df in kingdoms.items():
            bins.update({kingdom:MAG(self.assembly, df.index.tolist())})
        return bins

    def write_ranks(self, rank='superkingdom'):
        """Write fastas split by `rank`.

        Parameters
        ----------
        rank : str
            `rank` from canonical ranks (the default is 'superkingdom').

        Returns
        -------
        list
            [rank_name_fp, ...]

        Raises
        -------
        ValueError
            `rank` not in canonical ranks

        """
        if rank not in NCBI.CANONICAL_RANKS:
            raise ValueError(f'rank: {rank} not in {NCBI.CANONICAL_RANKS}')
        fpaths = []
        for rank_name,dff in self.taxonomy.groupby(rank):
            records = [r for r in SeqIO.parse(self.assembly, 'fasta') if r.id in dff.index]
            rank_name = rank_name.replace(' ','_')
            rank_name_fname = '.'.join([rank_name.title(),'fna'])
            rank_name_fp = os.path.join(self.outdir, rank_name_fname)
            if not records:
                logger.warning(f'No records to write to {rank_name_fp}')
            else:
                n_written = SeqIO.write(records, rank_name_fp, 'fasta')
                logger.debug(f'Wrote {n_written} records to {rank_name_fp}')
                fpaths.append(rank_name_fp)
        return fpaths

def main(args):
    raw_mags = Metagenome(args.assembly)
    mags = raw_mags.length_filter(cutoff=args.cutoff)
    logger.info(f'{args.cutoff:,}bp length filter {raw_mags.nseqs:,} to {mags.nseqs:,} seqs')
    logger.info(f'CallingORFs: force:{args.force}. cpus:{args.cpus} parallel:{args.noparallel}')
    try:
        mags.call_orfs(
            force=args.force,
            verbose=args.verbose,
            cpus=args.cpus,
            parallel=args.noparallel,
        )
    except FileExistsError as err:
        logger.warning(f'FileExistsError: {mags.prots_out}. Skipping...')
    logger.info(f'TaxonAssignment - method:{args.taxon_method}')
    mags.get_kingdoms(
        method=args.taxon_method,
        fasta=mags.prots_out,
        ncbi_dir=args.ncbi,
        outdir=mags.dirname,
        blast=None,
        usepickle=True,
        verbose=args.verbose,
    )
    logger.info(f'Getting k-mers of size {args.kmer_size}. Normalize: {args.kmer_normalized}')
    kmer_fpath = args.kmer_fpath if args.kmer_fpath else os.path.join(mags.dirname, 'kmers.tsv')
    kmers_df = mags.get_kmers(
        kmer_size=args.kmer_size,
        normalized=args.kmer_normalized,
        out=kmer_fpath,
        force=args.force,
    )
    logger.info('Done')

if __name__ == '__main__':
    import argparse
    import logging as logger
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logger.DEBUG)
    parser = argparse.ArgumentParser('Metagenome class to filter sequences by length')
    parser.add_argument('assembly', help='</path/to/assembly.fasta>')
    parser.add_argument('--ncbi', help='</path/to/ncbi/dir>', required=True)
    parser.add_argument('--cutoff', help='length to filter sequences',default=3000,
        type=int)
    parser.add_argument('--kmer-size', default=5, type=int)
    parser.add_argument(
        '--kmer-normalized',
        help='Perform CLR transform on k-mer frequencies if provided. (</path/to/kmers.normalized.tsv>)',
    )
    parser.add_argument(
        '--kmer-fpath',
        help='</path/to/kmers.tsv>',
    )
    parser.add_argument('--taxon-method', default='majority_vote',
        choices=['majority_vote'])
    parser.add_argument('--vote-fname', help='<taxon-vote.tsv filename>',
        default='taxonomy_vote.tsv')
    # Eventually will need to create a subparser for the taxon assignment methods
    # to include help information and required parameters according to method.
    parser.add_argument('--cpus', help='num cpus to use', default=0)
    parser.add_argument('--noparallel', help="Do not use GNU parallel",
        action='store_false', default=True)
    parser.add_argument('--force', help="overwrite existing files",
        action='store_true', default=False)
    parser.add_argument('--verbose', help="add verbosity",
        action='store_true', default=False)
    args = parser.parse_args()
    main(args)
