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

File containing customized AutometaErrors for more specific exception handling
"""


class AutometaError(Exception):
    """
    Base class for other exceptions
    """

    pass


class KmerFormatError(AutometaError):
    """
    Raised when either kmer_fpath file format is invalid or provided kmers or embedded are not formatted correctly for use.

    Parameters
    ----------
    AutometaError : class
        Base class for other exceptions
    """

    def __init__(self, fpath):
        self.fpath = fpath

    def __str__(self):
        """
        Operator overloading to return the path to k-mers frequency tab-delimited table or
        k-mers normalized tab-delimited table

        Returns
        -------
        str
            Path to k-mers frequency tab-delimited table or k-mers normalized tab-delimited table
        """
        return (
            f'{self.fpath} does not contain a "contig" column. '
            "Ensure the k-mer matrix was properly generated."
        )


class KmerEmbeddingError(AutometaError):
    """
    Raisede when either kmers or embedded are not provided.

    Parameters
    ----------
    AutometaError : class
        Base class for other exceptions
    """

    def __init__(self, value):
        self.value = value


class BinningError(AutometaError):
    """
    Raised when no marker information is availble for contigs to be binned or when dataframe is missing kmer/coverage annotations

    Parameters
    ----------
    AutometaError : class
        Base class for other exceptions
    """

    def __init__(self, value):
        self.value = value


class DatabaseOutOfSyncError(AutometaError):
    """
    Raised when NCBI databases nodes.dmp, names.dmp and merged.dmp are out of sync with each other

    Parameters
    ----------
    AutometaError : class
        Base class for other exceptions
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        """
        Operator overloading to return the text message written while raising the error,
        rather than the message of __str__ by base exception

        Returns
        -------
        str
            Message written alongside raising the exception
        """
        message = """

        NCBI databases nodes.dmp, names.dmp merged.dmp, prot.accession2taxid.gz and nr.gz are out of sync with each other

        Up-to-date databases may be downloaded at:

        non-redundant protein database (nr) - ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz

        prot.accession2taxid.gz - ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz

        taxdump.tar.gz - ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

        Required files within taxdump tarball are *nodes.dmp*, *names.dmp* and *merged.dmp*

        """
        return f"{self.value}\n{message}"


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="file containing utilities functions for Autometa pipeline"
    )
    print("file containing utilities functions for Autometa pipeline")
    args = parser.parse_args()
    import sys

    sys.exit(1)
