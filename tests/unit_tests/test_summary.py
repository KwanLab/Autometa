#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

Script to test autometa/binning/summary.py
"""

import pytest
import pandas as pd
from autometa.binning import summary
from autometa.taxonomy import majority_vote
import argparse


@pytest.fixture(name="bin_df", scope="module")
def fixture_bin_df(variables):
    binning_test_data = variables["summary"]
    return pd.read_json(binning_test_data["bin_df"]).set_index("contig")


@pytest.fixture(name="assembly", scope="module")
def fixture_assembly(tmp_path_factory, bin_df):
    fpath = tmp_path_factory.mktemp("summary") / "assembly.fna"
    contigs = bin_df.index.unique().tolist()
    lines = ""
    for contig in contigs:
        lines += f">{contig}\nGCGCGGGGCGGCGATACAATATA\n"
    fpath.write_text(lines)
    return fpath


@pytest.fixture(name="markers", scope="module")
def fixture_markers(variables):
    binning_test_data = variables["binning"]
    return pd.read_json(binning_test_data["markers"])


@pytest.fixture(name="markers_fpath")
def fixture_markers_fpath(markers, tmp_path):
    fpath = tmp_path / "markers.tsv"
    markers.to_csv(fpath, sep="\t", index=False, header=True)
    return fpath


@pytest.fixture(name="disjoint_markers_fpath")
def fixture_disjoint_markers_fpath(markers, bin_df, tmp_path):
    fpath = tmp_path / "disjoint_markers.tsv"
    markers[~markers.index.isin(bin_df.index)].iloc[:1, :].to_csv(
        fpath, sep="\t", index=False, header=True
    )
    return fpath


def return_mock_ncbi(*args, **kwargs):
    class MockNCBI:
        def get_lineage_dataframe(self, *args, **kwargs):
            df = pd.DataFrame.from_dict(
                {
                    "superkingdom": {
                        573569: "bacteria",
                        1302151: "bacteria",
                        37326: "bacteria",
                        60890: "bacteria",
                        309887: "bacteria",
                    },
                    "phylum": {
                        573569: "proteobacteria",
                        1302151: "firmicutes",
                        37326: "actinobacteria",
                        60890: "proteobacteria",
                        309887: "deinococcus-thermus",
                    },
                    "class": {
                        573569: "gammaproteobacteria",
                        1302151: "clostridia",
                        37326: "actinobacteria",
                        60890: "alphaproteobacteria",
                        309887: "deinococci",
                    },
                    "order": {
                        573569: "thiotrichales",
                        1302151: "clostridiales",
                        37326: "corynebacteriales",
                        60890: "rhodobacterales",
                        309887: "deinococcales",
                    },
                    "family": {
                        573569: "francisellaceae",
                        1302151: "clostridiaceae",
                        37326: "nocardiaceae",
                        60890: "rhodobacteraceae",
                        309887: "deinococcaceae",
                    },
                    "genus": {
                        573569: "francisella",
                        1302151: "caloranaerobacter",
                        37326: "nocardia",
                        60890: "phaeobacter",
                        309887: "deinococcus",
                    },
                    "species": {
                        573569: "francisella salina",
                        1302151: "caloranaerobacter sp. tr13",
                        37326: "nocardia brasiliensis",
                        60890: "phaeobacter gallaeciensis",
                        309887: "deinococcus maricopensis",
                    },
                }
            )
            df.index.name = "taxid"
            return df

    return MockNCBI()


@pytest.fixture(name="mock_rank_taxids")
def fixture_mock_rank_taxids(monkeypatch):
    def return_metabin_taxonomies(*args, **kwargs):
        return {
            "bin_0001": 1302151,
            "bin_0002": 37326,
            "bin_0003": 573569,
            "bin_0004": 309887,
            "bin_0005": 60890,
        }

    monkeypatch.setattr(
        majority_vote, "rank_taxids", return_metabin_taxonomies, raising=True
    )


@pytest.mark.skip
def test_get_metabin_taxonomies(
    mock_rank_taxids,
    bin_df,
):
    mock_ncbi = return_mock_ncbi()
    df = summary.get_metabin_taxonomies(bin_df=bin_df, ncbi=mock_ncbi)
    assert df.index.name == "cluster"
    assert "taxid" in df.columns
    ranks = [
        "species",
        "genus",
        "family",
        "order",
        "class",
        "phylum",
        "superkingdom",
    ]
    for rank in ranks:
        assert rank in df.columns


@pytest.mark.skip
def test_get_metabin_stats(bin_df, markers_fpath):
    df = summary.get_metabin_stats(bin_df=bin_df, markers_fpath=markers_fpath)
    assert df.index.name == "cluster"
    assert df.shape == (5, 20)


def test_get_metabin_stats_disjoint_markers(disjoint_markers_fpath, bin_df):
    df = summary.get_metabin_stats(bin_df=bin_df, markers_fpath=disjoint_markers_fpath)
    assert df.index.name == "cluster"
    assert df.shape == (5, 20)


@pytest.mark.skip
def test_write_cluster_records(bin_df, assembly, tmp_path):
    dirpath = tmp_path / "summary"
    summary.write_cluster_records(bin_df=bin_df, metagenome=assembly, outdir=dirpath)


@pytest.fixture(name="mock_parser")
def fixture_mock_parser(monkeypatch, bin_df, markers, assembly, tmp_path):
    binning_main = tmp_path / "binning_main.tsv"
    output_stats = tmp_path / "output_stats.tsv"
    output_taxonomy = tmp_path / "output_taxonomy.tsv"
    output_metabins = tmp_path / "metabins"
    bin_df.to_csv(binning_main, sep="\t", index=True, header=True)

    def return_mock_parser(*args, **kwargs):
        return MockParser()

    class MockArgs:
        def __init__(self):
            self.binning_main = binning_main
            self.markers = markers
            self.metagenome = assembly
            self.binning_column = "cluster"
            self.output_stats = output_stats
            self.output_taxonomy = output_taxonomy
            self.output_metabins = output_metabins
            self.ncbi = return_mock_ncbi

    class MockParser:
        def add_argument(self, *args, **kwargs):
            pass

        def parse_args(self):
            return MockArgs()

    # Defining the MockParser class to represent parser
    monkeypatch.setattr(argparse, "ArgumentParser", return_mock_parser, raising=True)
