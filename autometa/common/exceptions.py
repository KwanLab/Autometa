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
    KmerFormatError exception class.

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
    KmerEmbeddingError exception class.

    Parameters
    ----------
    AutometaError : class
        Base class for other exceptions
    """

    def __init__(self, value):
        self.value = value


class BinningError(AutometaError):
    """
    BinningError exception class.

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

        NCBI databases nodes.dmp, names.dmp and merged.dmp are out of sync with each other or not updated or both.

        Up-to-date databases may be downloaded at: ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

        Required files within taxdump tarball are *nodes.dmp*, *names.dmp* and *merged.dmp*

        """
        return message


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="file containing utilities functions for Autometa pipeline"
    )
    print("file containing utilities functions for Autometa pipeline")
    args = parser.parse_args()
    import sys

    sys.exit(1)
