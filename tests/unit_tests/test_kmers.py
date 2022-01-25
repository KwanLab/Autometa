#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

Script to test autometa/common/kmers.py
"""

import pytest
import argparse

import pandas as pd
import pytest

from autometa.common import kmers
from autometa.common.exceptions import TableFormatError


@pytest.fixture(name="assembly", scope="module")
def fixture_assembly(variables, tmp_path_factory):
    kmer_test_data = variables["metagenome"]
    records = kmer_test_data["assembly"]
    outlines = ""
    for record, seq in records.items():
        outlines += f"{record}\n{seq}\n"
    outdir = tmp_path_factory.mktemp("kmers")
    fpath = outdir / "metagenome.fna"
    with open(fpath, "w") as fh:
        fh.write(outlines)
    return str(fpath)


@pytest.fixture(name="counts", scope="module")
def fixture_counts(variables):
    kmer_test_data = variables["kmers"]
    df = pd.read_json(kmer_test_data["counts"])
    # kmer size is 5 (b/c this is the default).
    df.set_index("contig", inplace=True)
    return df


@pytest.fixture(name="counts_fpath", scope="module")
def fixture_counts_fpath(counts, tmp_path_factory):
    fpath = tmp_path_factory.mktemp("kmers") / "counts.tsv"
    counts.to_csv(fpath, sep="\t", index=True, header=True)
    return str(fpath)


@pytest.fixture(name="norm_df", scope="module")
def fixture_norm_df(variables):
    kmer_test_data = variables["kmers"]
    df = pd.read_json(kmer_test_data["am_clr_normalized_counts"])
    df.set_index("contig", inplace=True)
    return df


@pytest.fixture(name="invalid_df_fpath")
def fixture_invalid_df_fpath(tmp_path):
    invalid_dict = {
        "column1": ["invalid_contig_1", "invalid_contig_2", "invalid_contig_3"],
        "column2": ["invalid_marker1", "invalid_marker2", "invalid_marker3"],
    }
    df = pd.DataFrame(invalid_dict)
    df_fpath = tmp_path / "invalid_df.tsv"
    df.to_csv(df_fpath)
    return str(df_fpath)


def test_kmer_load(counts_fpath):
    df = kmers.load(kmers_fpath=counts_fpath)
    assert not df.empty
    assert df.index.name == "contig"


def test_kmer_load_missing_file():
    with pytest.raises(FileNotFoundError):
        kmers.load(kmers_fpath="Invalid_fpath")


def test_kmer_load_table_missing_contig_column(invalid_df_fpath):
    with pytest.raises(TableFormatError):
        kmers.load(invalid_df_fpath)


def test_count(assembly, tmp_path):
    out = tmp_path / "kmers.tsv"
    size = 5
    force = False
    df = kmers.count(assembly=assembly, size=size, out=out, force=force)
    assert df.shape[1] == 4 ** size / 2
    assert df.index.name == "contig"
    assert out.exists()


@pytest.mark.parametrize("force", [True, False])
def test_count_out_exists(assembly, counts, force, tmp_path):
    out = tmp_path / "kmers.tsv"
    counts.to_csv(out, sep="\t", index=True, header=True)
    size = 5
    df = kmers.count(assembly=assembly, size=size, out=out, force=force)
    assert df.shape[1] == 4 ** size / 2
    assert df.index.name == "contig"
    assert out.exists()


def test_count_wrong_size(assembly):
    size = 5.5
    with pytest.raises(TypeError):
        kmers.count(assembly=assembly, size=size)


@pytest.mark.parametrize("method", ["am_clr", "clr", "ilr"])
def test_normalize(counts, method, tmp_path):
    out = tmp_path / "kmers.norm.tsv"
    force = False
    df = kmers.normalize(df=counts, method=method, out=out, force=force)
    if method in {"am_clr", "clr"}:
        assert df.shape == counts.shape
    else:
        # ILR will reduce the columns by one.
        assert df.shape[1] < counts.shape[1]
    assert out.exists()


@pytest.mark.parametrize("force", [True, False])
def test_normalize_out_exists(counts, norm_df, force, tmp_path):
    out = tmp_path / "kmers.norm.tsv"
    norm_df.to_csv(out, sep="\t", index=True, header=True)
    df = kmers.normalize(df=counts, method="am_clr", out=out, force=force)
    assert df.shape == counts.shape
    assert df.index.name == "contig"


def test_normalize_wrong_method(counts, tmp_path):
    out = tmp_path / "kmers.norm.tsv"
    with pytest.raises(ValueError):
        kmers.normalize(df=counts, method="am_ilr", out=out, force=False)


@pytest.mark.parametrize("method", ["bhsne", "sksne", "umap", "densmap"])
def test_embed_methods(norm_df, method, tmp_path):
    seed = 42
    out = tmp_path / "kmers.embed.tsv"
    force = False
    embed_dimensions = 2
    pca_dimensions = 3
    df = kmers.embed(
        kmers=norm_df,
        out=out,
        force=force,
        embed_dimensions=embed_dimensions,
        pca_dimensions=pca_dimensions,
        method=method,
        seed=seed,
    )
    assert df.shape[1] == embed_dimensions


@pytest.mark.parametrize("embed_dimensions", [2, 3, 4])
def test_embed_dimensions(norm_df, embed_dimensions, tmp_path):
    out = tmp_path / "kmers.embed.tsv"
    df = kmers.embed(
        kmers=norm_df,
        out=out,
        force=False,
        embed_dimensions=embed_dimensions,
        pca_dimensions=5,
        method="bhsne",
        seed=42,
    )
    assert df.shape[1] == embed_dimensions


@pytest.mark.parametrize("force", [False, True])
def test_embed_out_exists(norm_df, force, tmp_path):
    seed = 42
    out = tmp_path / "kmers.embed.tsv"
    method = "bhsne"
    embed_dimensions = 2
    pca_dimensions = 3
    df = kmers.embed(
        kmers=norm_df,
        out=out,
        force=force,
        embed_dimensions=embed_dimensions,
        pca_dimensions=pca_dimensions,
        method=method,
        seed=seed,
    )


def test_embed_TableFormatError(invalid_df_fpath):
    with pytest.raises(TableFormatError):
        kmers.embed(kmers=invalid_df_fpath)


def test_embed_input_not_string_or_dataframe(tmp_path):
    kmer_fpath = tmp_path / "kmers.embed.tsv"
    with pytest.raises(TypeError):
        kmers.embed(kmers=kmer_fpath)


def test_embed_empty_dataframe(tmp_path):
    empty_df = pd.DataFrame({})
    out = tmp_path / "kmers.embed.tsv"
    with pytest.raises(FileNotFoundError):
        kmers.embed(kmers=empty_df, out=out, force=True)


@pytest.mark.entrypoint
def test_kmers_main(monkeypatch, tmp_path, assembly):
    norm_method = "am_clr"
    embed_method = "bhsne"
    counts_out = tmp_path / "kmers.tsv"
    normalized = tmp_path / f"kmers.{norm_method}.tsv"
    embedded = tmp_path / f"kmers.{norm_method}.{embed_method}.tsv"
    embed_dimensions = 2

    class MockArgs:
        def __init__(self):
            self.fasta = assembly
            self.size = 4
            self.kmers = counts_out
            self.force = True
            self.norm_method = norm_method
            self.norm_output = normalized
            self.embedding_method = embed_method
            self.embedding_dimensions = embed_dimensions
            self.embedding_output = embedded
            self.pca_dimensions = 3
            self.cpus = 1
            self.seed = 42

    class MockParser:
        def add_argument(self, *args, **kwargs):
            pass

        def parse_args(self):
            return MockArgs()

    def return_mock_parser(*args, **kwargs):
        return MockParser()

    monkeypatch.setattr(argparse, "ArgumentParser", return_mock_parser, raising=True)
    kmers.main()
    assert embedded.exists()
    df = pd.read_csv(embedded, sep="\t")
    assert "contig" in df.columns
    assert f"x_{embed_dimensions}" in df.columns
    # Make sure we have all of our embedding dimensions and need to account for our contig column
    assert embed_dimensions + 1 == df.shape[1]
