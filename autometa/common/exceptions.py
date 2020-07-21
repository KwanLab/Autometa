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


class AutometaError(Exception):
    """Base class for Autometa Errors."""

    pass


class TableFormatError(AutometaError):
    """TableFormatError exception class.

    Exception called when Table format is incorrect.

    This is usually a result of a table missing the 'contig' column as this is
    often used as the index.
    """

    pass


class ChecksumMismatchError(AutometaError):
    """ChecksumMismatchError exception class

    Exception called when checksums do not match.

    """

    pass


if __name__ == "__main__":
    print(
        "This file contains Exceptions for the Autometa pipeline and should not be run directly!"
    )
