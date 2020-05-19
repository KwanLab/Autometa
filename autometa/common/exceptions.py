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

File containing customized AutometaExceptions for more specific exception handling
"""


class AutometaException(Exception):
    """docstring for AutometaException."""

    def __init__(self, value):
        self.value = value

    issue_request = '''
    An error was encountered!

    Please help us fix your problem!

    You may file an issue with us at https://github.com/KwanLab/Autometa/issues/new
    '''

    def __str__(self):
        return f'{self.value}\n\n{issue_request}'


class KmerFormatError(AutometaException):
    """KmerFormatError exception class."""

    def __init__(self, fpath):
        super(AutometaException, self).__init__(fpath)
        self.fpath = fpath

    def __str__(self):
        return f'{self.fpath} does not contain a \"contig\" column. '\
            'Ensure the k-mer matrix was properly generated.'


class KmerEmbeddingError(AutometaException):
    """KmerEmbeddingError exception class."""

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value


class BinningError(AutometaException):
    """BinningError exception class."""

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='file containing utilities functions for Autometa pipeline')
    print('file containing utilities functions for Autometa pipeline')
    args = parser.parse_args()
    import sys
    sys.exit(1)
