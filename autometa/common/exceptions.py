#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
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

File containing customized AutometaExceptions for more specific exception handling
"""

class AutometaException(Exception):
    """docstring for AutometaException."""
    issue_request = '''
    An error was encountered!

    Please help us fix your problem!

    You may file an issue with us at https://github.com/KwanLab/Autometa/issues/new
    '''
    pass


class KmerFormatError(Exception):
    """KmerFormatError exception class."""

    def __init__(self, fpath):
        self.fpath = fpath

    def __str__(self):
        return f'{self.fpath} does not contain a \"contig\" column. '\
        'Ensure the k-mer matrix was properly generated.'

class KmerEmbeddingError(Exception):
    """KmerEmbeddingError exception class."""

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value

class RecursiveDBSCANError(Exception):
    """RecursiveDBSCANError exception class."""

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value

if __name__ == '__main__':
    print('This file contains custom exceptions for Autometa and should not be run directly')
    import sys;sys.exit(1)
