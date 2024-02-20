#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

Script to test autometa/common/metagenome.py
"""

import pandas as pd
import pytest
import argparse
import os

from Bio import SeqIO

from autometa.common import metagenome as am_metagenome
from autometa.common.metagenome import Metagenome
from autometa.common import metagenome
from autometa.common.external import prodigal


@pytest.fixture(name="metagenome_testdir", scope="module")
def fixture_metagenome_testdir(tmp_path_factory):
    return tmp_path_factory.mktemp("metagenome")


@pytest.fixture(name="assembly", scope="module")
def fixture_assembly(variables, metagenome_testdir):
    metagenome_test_data = variables["metagenome"]
    records = metagenome_test_data["assembly"]
    lines = ""
    for record, seq in records.items():
        lines += f"{record}\n{seq}\n"
    fpath = metagenome_testdir / "assembly.fna"
    with open(fpath, "w") as fh:
        fh.write(lines)
    return str(fpath)


@pytest.fixture(name="metagenome", scope="module")
def fixture_metagenome(assembly):
    return Metagenome(assembly=assembly)


@pytest.mark.parametrize("quality_measure", [0.1, 0.9])
def test_fragmentation_metric(metagenome, quality_measure):
    metric = metagenome.fragmentation_metric(quality_measure=quality_measure)
    if quality_measure == 0.1:
        expected = 7237
    elif metagenome.nseqs == 3 and quality_measure == 0.9:
        expected = 7231
    elif metagenome.nseqs == 4 and quality_measure == 0.9:
        expected = 7229
    else:
        expected = None

    assert isinstance(metric, int)
    if expected:
        assert metric == expected


def test_fragmentation_metric_wrong_quality_measure_value(metagenome):
    with pytest.raises(ValueError):
        metagenome.fragmentation_metric(quality_measure=1.4)


def test_length_weighted_gc(metagenome):
    if metagenome.nseqs == 3:
        expected = 65.3167472932504
    elif metagenome.nseqs == 4:
        expected = 63.19554848966614
    else:
        expected = None
    assert isinstance(metagenome.length_weighted_gc, float)
    if expected:
        assert metagenome.length_weighted_gc == expected


def test_largest_seq(metagenome):
    expected = "NODE_1501_length_7237_cov_222.15"
    assert metagenome.largest_seq == expected


def test_size(metagenome):
    if metagenome.nseqs == 3:
        expected = 21705
    elif metagenome.nseqs == 4:
        expected = 28934
    else:
        expected = None
    assert isinstance(metagenome.size, int)
    if expected:
        assert metagenome.size == expected


def test_describe(metagenome):
    df = metagenome.describe()
    assert isinstance(df, pd.DataFrame)
    assert df.index.name == "assembly"
    expected_cols = [
        "nseqs",
        "size (bp)",
        "N50 (bp)",
        "N10 (bp)",
        "N90 (bp)",
        "length_weighted_gc_content (%)",
        "largest_seq",
    ]
    for col in expected_cols:
        assert col in df.columns


@pytest.mark.parametrize("force", [False, True])
def test_force_length_filter(metagenome, force):
    existing = str(metagenome)
    if not force:
        with pytest.raises(FileExistsError):
            filtered_metagenome = metagenome.length_filter(
                out=existing, cutoff=7232, force=force
            )
    else:
        filtered_metagenome = metagenome.length_filter(
            out=existing, cutoff=7232, force=True
        )
        assert filtered_metagenome.nseqs == 2


@pytest.mark.parametrize("cutoff", [-6, "12324", 7232, 10000])
def test_length_filter_invalid_cutoff(metagenome, cutoff, tmp_path):
    out = tmp_path / "length_filtered.fna"
    expected_errors = {-6: ValueError, "12324": TypeError}
    if cutoff in expected_errors:
        expected_error = expected_errors[cutoff]
        with pytest.raises(expected_error):
            filtered_metagenome = metagenome.length_filter(
                out=out, cutoff=cutoff, force=False
            )
    else:
        filtered_metagenome = metagenome.length_filter(
            out=out, cutoff=cutoff, force=False
        )
    if cutoff == 7232:
        assert filtered_metagenome.nseqs == 2
    elif cutoff == 10000:
        # No contigs above cutoff. Ensure we return original
        assert filtered_metagenome.nseqs == metagenome.nseqs


@pytest.mark.entrypoint
def test_metagenome_main(monkeypatch, assembly, tmp_path):
    output_fasta = tmp_path / "metagenome.filtered.fna"
    output_stats = tmp_path / "metagenome.stats.tsv"
    output_gc_content = tmp_path / "metagenome.gc_content.tsv"

    class MockArgs:
        def __init__(self):
            self.assembly = assembly
            self.force = True
            self.cutoff = 7232
            self.output_fasta = output_fasta
            self.output_gc_content = output_gc_content
            self.output_stats = output_stats

    # Defining the MockParser class to represent parser
    class MockParser:
        def add_argument(self, *args, **kwargs):
            pass

        def parse_args(self):
            return MockArgs()

    def return_mock_parser(*args, **kwargs):
        return MockParser()

    monkeypatch.setattr(argparse, "ArgumentParser", return_mock_parser, raising=True)

    am_metagenome.main()
    assert output_fasta.exists()
    assert output_gc_content.exists()
    assert output_stats.exists()
    num_expected_seqs = 2
    assert (
        len([record for record in SeqIO.parse(output_fasta, "fasta")])
        == num_expected_seqs
    )
