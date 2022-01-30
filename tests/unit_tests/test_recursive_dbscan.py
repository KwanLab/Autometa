#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

Script to test autometa/binning/recursive_dbscan.py
"""

import pytest
import argparse

import pandas as pd

from autometa.binning import recursive_dbscan
from autometa.common.exceptions import TableFormatError
from autometa.common.markers import load as load_markers


@pytest.fixture(name="binning_testdir", scope="module")
def fixture_binning_testdir(tmp_path_factory):
    return tmp_path_factory.mktemp("binning")


@pytest.fixture(name="test_contigs", scope="module")
def fixture_test_contigs(variables, markers, n=5, seed=42):
    binning_test_data = variables["binning"]
    df = pd.read_json(binning_test_data["kmers_normalized"])
    df.set_index("contig", inplace=True)
    markers_index = (
        df[df.index.isin(markers.index)].sample(n=n, random_state=seed).index
    )
    kmers_index = df.sample(n=n, random_state=seed).index
    return markers_index.union(kmers_index)


@pytest.fixture(name="kmers", scope="module")
def fixture_kmers(variables, test_contigs, binning_testdir):
    binning_test_data = variables["binning"]
    df = pd.read_json(binning_test_data["kmers_embedded"])
    fpath = binning_testdir / "kmers.embed.tsv"
    df.set_index("contig", inplace=True)
    df.loc[test_contigs].to_csv(fpath, sep="\t", index=True, header=True)
    return fpath


@pytest.fixture(name="coverage", scope="module")
def fixture_coverage(variables, test_contigs, binning_testdir):
    binning_test_data = variables["binning"]
    df = pd.read_json(binning_test_data["coverage"])
    df.set_index("contig", inplace=True)
    fpath = binning_testdir / "coverage.tsv"
    df.loc[test_contigs].to_csv(fpath, sep="\t", index=True, header=True)
    return fpath


@pytest.fixture(name="gc_content", scope="module")
def fixture_gc_content(variables, test_contigs, binning_testdir):
    binning_test_data = variables["binning"]
    df = pd.read_json(binning_test_data["gc_content"])
    df.set_index("contig", inplace=True)
    fpath = binning_testdir / "gc_content.tsv"
    df.loc[test_contigs].to_csv(fpath, sep="\t", index=True, header=True)
    return fpath


@pytest.fixture(name="taxonomy", scope="module")
def fixture_taxonomy(variables, test_contigs, binning_testdir):
    binning_test_data = variables["binning"]
    df = pd.read_json(binning_test_data["taxonomy"])
    fpath = binning_testdir / "taxonomy.tsv"
    df.set_index("contig", inplace=True)
    df.loc[test_contigs].to_csv(fpath, sep="\t", index=True, header=True)
    return fpath


@pytest.fixture(name="markers_fpath", scope="module")
def fixture_markers_fpath(variables, binning_testdir):
    binning_test_data = variables["binning"]
    df = pd.read_json(binning_test_data["markers"])
    fpath = binning_testdir / "markers.tsv"
    df.to_csv(fpath, sep="\t", index=False, header=True)
    return str(fpath)


@pytest.fixture(name="disjoint_markers_fpath", scope="module")
def fixture_disjoint_markers_fpath(binning_testdir):
    df = pd.DataFrame(
        {
            "contig": ["invalid_contig", "invalid_contig", "invalid_contig"],
            "sacc": ["sacc", "sacc", "sacc"],
        }
    )
    fpath = binning_testdir / "invalid_markers.tsv"
    df.to_csv(fpath, sep="\t", index=False, header=True)
    return str(fpath)


@pytest.fixture(name="markers", scope="module")
def fixture_markers(markers_fpath):
    return load_markers(markers_fpath)


@pytest.fixture(name="main_df", scope="module")
def fixture_main_df(kmers, coverage, gc_content, taxonomy):
    main_df = pd.read_csv(kmers, sep="\t", index_col="contig")
    for fpath in [coverage, gc_content, taxonomy]:
        df = pd.read_csv(fpath, sep="\t", index_col="contig")
        main_df = pd.merge(main_df, df, how="left", right_index=True, left_index=True)
    main_df = main_df.convert_dtypes()
    return main_df


def test_taxon_guided_binning(main_df, markers):
    num_contigs = main_df.shape[0]
    df = recursive_dbscan.taxon_guided_binning(
        main=main_df,
        markers=markers,
        completeness=10.0,
        purity=80.0,
        coverage_stddev=25.0,
        gc_content_stddev=5.0,
        method="dbscan",
    )
    assert isinstance(df, pd.DataFrame)
    assert "cluster" in df.columns
    assert "purity" in df.columns
    assert "completeness" in df.columns
    assert "coverage_stddev" in df.columns
    assert "gc_content_stddev" in df.columns
    assert df.shape[0] == num_contigs


@pytest.mark.parametrize("method", ["dbscan", "hdbscan"])
def test_get_clusters(main_df, markers, method):
    num_contigs = main_df.shape[0]
    df = recursive_dbscan.get_clusters(
        main=main_df,
        markers_df=markers,
        completeness=20.0,
        purity=95.0,
        coverage_stddev=25.0,
        gc_content_stddev=5.0,
        method=method,
    )
    assert isinstance(df, pd.DataFrame)
    assert "cluster" in df.columns
    assert "purity" in df.columns
    assert "completeness" in df.columns
    assert "coverage_stddev" in df.columns
    assert "gc_content_stddev" in df.columns
    assert df.shape[0] == num_contigs


def test_binning_invalid_clustering_method(main_df, markers):
    with pytest.raises(ValueError):
        recursive_dbscan.get_clusters(
            main=main_df,
            markers_df=markers,
            completeness=20.0,
            purity=95.0,
            coverage_stddev=5.0,
            gc_content_stddev=5.0,
            method="invalid_clustering_method",
        )


@pytest.mark.entrypoint
def test_recursive_dbscan_main(
    monkeypatch,
    kmers,
    coverage,
    gc_content,
    markers_fpath,
    taxonomy,
    tmp_path,
):
    output_binning = tmp_path / "binning.tsv"
    output_main = tmp_path / "binning_main.tsv"

    class MockArgs:
        def __init__(self):
            self.domain = "bacteria"
            self.kmers = kmers
            self.coverages = coverage
            self.gc_content = gc_content
            self.markers = markers_fpath
            self.output_binning = output_binning
            self.output_main = output_main
            self.clustering_method = "dbscan"
            self.completeness = 20.0
            self.purity = 95.0
            self.cov_stddev_limit = 25.0
            self.gc_stddev_limit = 5.0
            self.taxonomy = taxonomy
            self.rank_filter = "superkingdom"
            self.rank_name_filter = "bacteria"
            self.starting_rank = "superkingdom"
            self.reverse_ranks = False
            self.verbose = True
            self.cpus = -1

    class MockParser:
        def add_argument(self, *args, **kwargs):
            pass

        def parse_args(self):
            return MockArgs()

    def return_mock_parser(*args, **kwargs):
        return MockParser()

    monkeypatch.setattr(argparse, "ArgumentParser", return_mock_parser, raising=True)
    recursive_dbscan.main()
    assert output_binning.exists()
    df = pd.read_csv(output_binning, sep="\t")
    assert "contig" in df.columns
    assert "cluster" in df.columns


@pytest.mark.entrypoint
def test_recursive_dbscan_main_tableformaterror(
    monkeypatch,
    kmers,
    coverage,
    gc_content,
    disjoint_markers_fpath,
    taxonomy,
    tmp_path,
):
    output_binning = tmp_path / "binning.tsv"
    output_main = tmp_path / "binning_main.tsv"

    class MockArgs:
        def __init__(self):
            self.domain = "bacteria"
            self.kmers = kmers
            self.coverages = coverage
            self.gc_content = gc_content
            self.markers = disjoint_markers_fpath
            self.output_binning = output_binning
            self.output_main = output_main
            self.clustering_method = "dbscan"
            self.completeness = 20.0
            self.purity = 95.0
            self.cov_stddev_limit = 25.0
            self.gc_stddev_limit = 5.0
            self.taxonomy = taxonomy
            self.rank_filter = "superkingdom"
            self.rank_name_filter = "bacteria"
            self.starting_rank = "superkingdom"
            self.reverse_ranks = False
            self.verbose = True
            self.cpus = -1

    class MockParser:
        def add_argument(self, *args, **kwargs):
            pass

        def parse_args(self):
            return MockArgs()

    def return_mock_parser(*args, **kwargs):
        return MockParser()

    monkeypatch.setattr(argparse, "ArgumentParser", return_mock_parser, raising=True)
    with pytest.raises(TableFormatError):
        recursive_dbscan.main()
