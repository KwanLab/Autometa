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

AutometaUser configuration class

"""


import logging
import os

from autometa import config
from autometa.config import environ
from autometa.common import utilities

from autometa.common.metagenome import Metagenome
from autometa.common.metabin import MetaBin
from autometa.config.databases import Databases
from autometa.config.project import Project
from autometa.common.utilities import timeit


logger = logging.getLogger(__name__)


class AutometaUser:
    """AutometaUser Class to handle job submissions.

    Parameters
    ----------
    user_config : str, optional
        </path/to/user/default.config (the default is config.DEFAULT_FPATH).
    nproc : int, optional
        Number of cpus to use when downloading/updating databases. (the default is 2).

    Attributes
    ----------
    home_dir : str
        </path/to/Autometa/>
    config : config.ConfigParser
        User configuration from `user_config`.

    """

    def __init__(self, user_config=config.DEFAULT_FPATH, nproc=2):
        self.home_dir = config.set_home_dir()
        self.nproc = nproc
        self.user_config = user_config
        self.config = config.get_config(self.user_config)
        if not self.config.has_section("common"):
            self.config.add_section("common")
        if not self.config.has_option("common", "home_dir"):
            self.config.set("common", "home_dir", self.home_dir)

    def __str__(self):
        return self.user_config

    def save(self):
        """Saves the current user config to `self.user_config` file path."""
        config.put_config(self.config, self.user_config)

    def configure(self, dryrun=True, update=False):
        """Configure user execution environment and databases.

        Parameters
        ----------
        dryrun : bool
            Log configuration without performing updates (Default is True).

        Returns
        -------
        NoneType

        """

        # Execution env
        self.config, exe_satisfied = environ.configure(self.config)
        logger.info(f"Executable dependencies satisfied: {exe_satisfied}")
        # Database env
        dbs = Databases(self.config, dryrun=dryrun, nproc=self.nproc, update=update)
        no_checksum = not update
        self.config = dbs.configure(no_checksum=no_checksum)
        logger.info(f"Database dependencies satisfied: {dbs.satisfied()}")
        if dryrun:
            return

        if not dbs.satisfied():
            raise LookupError("Database dependencies not satisfied!")
        if not exe_satisfied:
            raise LookupError("Executable dependencies not satisfied!")

        self.save()

    def new_project(self, fpath):
        """Configure new project at `fpath`.

        Parameters
        ----------
        fpath : str
            </path/to/workspace/project_<num>/project.config>

        Returns
        -------
        autometa.config.project.Project
            Project object containing methods used to manipulate autometa project
        """
        dpath = os.path.dirname(fpath)
        if not os.path.exists(dpath):
            os.makedirs(dpath)
        config.put_config(self.config, fpath)
        return Project(fpath)

    def prepare_binning_args(self, fpath):
        """Prepares metagenome binning run using provided `fpath`.

        Notes
        -----
            This method performs a number of configuration checks to ensure the
            binning run will perform without conflicts.

            1. workspace check: Will construct workspace directory if provided does not exist.
            2. Project check: Will configure a new project if project number is not found in workspace directory.
            3. Metagenome check: Will update if existing with edits or resume if existing without edits. Otherwise will add new metagenome to project.


        Example
        -------
        .. code-block:: python

            #:  Generate namespace - mgargs.files.<file>
            mgargs = prepare_binning_args(mg_config)
            #:  Access file from args - mgargs.files.<file>
            mgargs.files.length_filtered
            #:  Access parameter from args - mgargs.parameters.<parameter>
            mgargs.parameters.length_cutoff

        Parameters
        ----------
        fpath : str
            </path/to/metagenome.config>

        Returns
        -------
        argparse.Namespace
            parameters and files parsed from config.

        Raises
        -------
        ExceptionName
            Why the exception is raised.

        """
        # 1. configure user environment
        self.configure()
        # 2. check workspace exists
        mgargs = config.parse_config(fpath)
        workspace = os.path.realpath(mgargs.parameters.workspace)
        if not os.path.exists(workspace):
            os.makedirs(workspace)
        # 3. check project exists
        proj_name = f"project_{mgargs.parameters.project:03d}"
        project_dirpath = os.path.realpath(os.path.join(workspace, proj_name))
        project_config_fp = os.path.join(project_dirpath, "project.config")
        if not os.path.exists(project_dirpath) or not os.path.exists(project_config_fp):
            project = self.new_project(project_config_fp)
        else:
            project = Project(project_config_fp)
        # 4. check whether existing or new run with metagenome_num
        metagenome = f"metagenome_{mgargs.parameters.metagenome_num:03d}"
        if metagenome not in project.metagenomes:
            mgargs = project.add(fpath)
            project.save()
            return mgargs
        # If resuming existing metagenome run. Check whether config file has changed.
        old_config_fpath = project.metagenomes.get(metagenome)
        old_chksum = utilities.get_checksum(old_config_fpath)
        new_chksum = utilities.get_checksum(fpath)
        if old_chksum != new_chksum:
            mgargs = project.update(
                metagenome_num=mgargs.parameters.metagenome_num, fpath=fpath
            )
        project.save()
        return mgargs

    @utilities.timeit
    def run_binning(self, mgargs):
        """Run the autometa metagenome binning pipeline using the provided metagenome args.

        Parameters
        ----------
        mgargs : argparse.Namespace
            metagenome args

        Returns
        -------
        NoneType

        Raises
        -------
        TODO: Need to enumerate all exceptions raised from within binning pipeline.
        I.e. Demarkate new exception (not yet handled) vs. handled exception.
        Subclassing an AutometaException class may be most appropriate use case here.
        """
        mg = Metagenome(
            assembly=mgargs.files.metagenome,
            outdir=mgargs.parameters.outdir,
            nucl_orfs_fpath=mgargs.files.nucleotide_orfs,
            prot_orfs_fpath=mgargs.files.amino_acid_orfs,
            taxonomy_fpath=mgargs.files.taxonomy,
            fwd_reads=mgargs.files.fwd_reads,
            rev_reads=mgargs.files.rev_reads,
            se_reads=mgargs.files.se_reads,
            taxon_method=mgargs.parameters.taxon_method,
        )

        try:
            # Original (raw) file should not be manipulated so return new object
            mg = mg.length_filter(
                out=mgargs.files.length_filtered, cutoff=mgargs.parameters.length_cutoff
            )
            # COMBAK: Checkpoint length filtered
        except FileExistsError as err:
            # COMBAK: Checkpoint length filtered
            logger.debug(f"{mgargs.files.length_filtered} already exists. Continuing..")
            mg = Metagenome(
                assembly=mgargs.files.length_filtered,
                outdir=mgargs.parameters.outdir,
                nucl_orfs_fpath=mgargs.files.nucleotide_orfs,
                prot_orfs_fpath=mgargs.files.amino_acid_orfs,
                taxonomy_fpath=mgargs.files.taxonomy,
                fwd_reads=mgargs.files.fwd_reads,
                rev_reads=mgargs.files.rev_reads,
                se_reads=mgargs.files.se_reads,
                taxon_method=mgargs.parameters.taxon_method,
            )
        # TODO asynchronous execution here (work-queue/Makeflow tasks)
        # TODO kmer counting here should just be employed via the kmers.py file as is done later on with kmers.embed(...)
        # I.e. kmers.count(...)
        mg.get_kmers(
            kmer_size=mgargs.parameters.kmer_size,
            normalized=mgargs.files.kmer_normalized,
            out=mgargs.files.kmer_counts,
            multiprocess=mgargs.parameters.kmer_multiprocess,
            nproc=mgargs.parameters.cpus,
            force=mgargs.parameters.force,
        )
        # COMBAK: Checkpoint kmers
        coverages = mg.get_coverages(
            out=mgargs.files.coverages,
            from_spades=mgargs.parameters.cov_from_spades,
            sam=mgargs.files.sam,
            bam=mgargs.files.bam,
            lengths=mgargs.files.lengths,
            bed=mgargs.files.bed,
        )
        # COMBAK: Checkpoint coverages

        if mgargs.parameters.do_taxonomy:
            # Filter by Kingdom
            kingdoms = mg.get_kingdoms(
                ncbi=mgargs.databases.ncbi,
                usepickle=mgargs.parameters.usepickle,
                blast=mgargs.files.blastp,
                hits=mgargs.files.blastp_hits,
                force=mgargs.parameters.force,
                cpus=mgargs.parameters.cpus,
            )

            if not mgargs.parameters.kingdom in kingdoms:
                recovered_kingdoms = ", ".join(kingdoms.keys())
                raise KeyError(
                    f"{mgargs.parameters.kingdom} not recovered in dataset. Recovered: {recovered_kingdoms}"
                )

            mag = kingdoms.get(mgargs.parameters.kingdom)
        else:
            mag = MetaBin(assembly=mg.assembly, contigs=mg.seqrecords, outdir=mg.outdir)

        # Perform binning
        # TODO: replace MetaBin method out with recursive_dbscan method recursive_dbscan.binning()
        # No reason to have this as a MetaBin method, since almost all parameters are available here.
        bins_df = mag.get_binning(
            method=mgargs.parameters.binning_method,
            kmers=mgargs.files.kmer_normalized,
            embedded=mgargs.files.kmer_embedded,
            do_pca=mgargs.parameters.do_pca,
            pca_dims=mgargs.parameters.pca_dims,
            embedding_method=mgargs.parameters.embedding_method,
            coverage=mgargs.files.coverages,
            domain=mgargs.parameters.kingdom,
            taxonomy=mgargs.files.taxonomy,
            clustering_method=mgargs.parameters.clustering_method,
            starting_rank=mgargs.parameters.starting_rank,
            reverse_ranks=mgargs.parameters.reverse_ranks,
        )
        binning_cols = ["cluster", "completeness", "purity"]
        bins_df[binning_cols].to_csv(
            mgargs.files.binning, sep="\t", index=True, header=True
        )


def main():
    import argparse
    import multiprocessing as mp
    import logging as logger

    parser = argparse.ArgumentParser(
        description="""
    Configures the Autometa user environment/databases.
    Running without args will download and format Autometa database dependencies.
    """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--config",
        help=f"Path to an Autometa default.config file.",
        default=config.DEFAULT_FPATH,
    )
    parser.add_argument(
        "--dryrun",
        help="Log configuration without performing updates.",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--debug", help="Stream debugging information to terminal.", action="store_true"
    )
    parser.add_argument(
        "--cpus",
        help=f"Num. cpus to use when updating/constructing databases.",
        default=mp.cpu_count(),
        type=int,
    )
    args = parser.parse_args()

    level = logger.DEBUG if args.debug else logger.INFO
    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=level,
    )
    user = AutometaUser(user_config=args.config, nproc=args.cpus)
    user.configure(args.dryrun)


if __name__ == "__main__":
    main()
