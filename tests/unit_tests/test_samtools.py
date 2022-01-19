#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

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
