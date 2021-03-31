#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
COPYRIGHT
Copyright 2021 Ian J. Miller, Evan R. Rees, Kyle Wolf, Siddharth Uppal,
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

Script to test autometa/common/external/samtools.py
"""

import pytest
from autometa.common.external import samtools


def test_sort_missing_file():
    with pytest.raises(FileNotFoundError):
        samtools.sort(sam="sam", bam="bam")


@pytest.mark.parametrize("cpus", [2.9, -2])
def test_sort_invalid_cpu_input(cpus):
    with pytest.raises(TypeError):
        samtools.sort(sam="sam_fpath", bam="bam", cpus=cpus)
