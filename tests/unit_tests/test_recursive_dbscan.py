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
    df = pd.read_json(binning_test_data["kmers_normalized"])
    fpath = binning_testdir / "kmers.norm.tsv"
    df.set_index("contig", inplace=True)
    df.loc[test_contigs].to_csv(fpath, sep="\t", index=True, header=True)
    return str(fpath)


@pytest.fixture(name="embedded_kmers", scope="module")
def fixture_embedded_kmers(variables, test_contigs, binning_testdir):
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


@pytest.fixture(name="markers", scope="module")
def fixture_markers(markers_fpath):
    return load_markers(markers_fpath)


@pytest.fixture(name="master", scope="module")
def fixture_master(embedded_kmers, coverage, taxonomy):
    master = pd.read_csv(embedded_kmers, sep="\t", index_col="contig")
    for fpath in [coverage, taxonomy]:
        df = pd.read_csv(fpath, sep="\t", index_col="contig")
        master = pd.merge(master, df, how="left", right_index=True, left_index=True)
    master = master.convert_dtypes()
    return master


@pytest.mark.parametrize("usetaxonomy", [True, False])
@pytest.mark.parametrize("method", ["dbscan", "hdbscan"])
def test_binning(master, markers, usetaxonomy, method):
    num_contigs = master.shape[0]
    df = recursive_dbscan.binning(
        master=master,
        markers=markers,
        taxonomy=usetaxonomy,
        starting_rank="superkingdom",
        reverse_ranks=False,
        domain="bacteria",
        completeness=20.0,
        purity=90.0,
        method=method,
        verbose=True,
    )
    assert isinstance(df, pd.DataFrame)
    assert "cluster" in df.columns
    assert "purity" in df.columns
    assert "completeness" in df.columns
    assert df.shape[0] == num_contigs


def test_binning_invalid_clustering_method(master, markers):
    with pytest.raises(ValueError):
        recursive_dbscan.binning(
            master=master,
            markers=markers,
            taxonomy=False,
            starting_rank="superkingdom",
            reverse_ranks=False,
            domain="bacteria",
            completeness=20.0,
            purity=90.0,
            method="invalid_clustering_method",
            verbose=False,
        )


def test_binning_empty_markers_table(master):
    invalid_dict = {
        "contig": ["invalid_contig_1", "invalid_contig_2", "invalid_contig_3"],
        "markers": ["invalid_marker1", "invalid_marker2", "invalid_marker3"],
    }
    df = pd.DataFrame(invalid_dict)
    with pytest.raises(TableFormatError):
        recursive_dbscan.binning(
            master=master,
            markers=df,
            domain="bacteria",
            method="hdbscan",
            taxonomy=False,
        )


@pytest.mark.entrypoint
def test_recursive_dbscan_main(
    monkeypatch, kmers, coverage, markers_fpath, embedded_kmers, taxonomy, tmp_path
):
    out = tmp_path / "binning.tsv"

    class MockArgs:
        def __init__(self):
            self.domain = "bacteria"
            self.kmers = kmers
            self.coverage = coverage
            self.markers = markers_fpath
            self.output = out
            self.embedded_kmers = embedded_kmers
            self.embedding_method = "bhsne"
            self.clustering_method = "dbscan"
            self.completeness = 20.0
            self.purity = 90.0
            self.taxonomy = taxonomy
            self.starting_rank = "superkingdom"
            self.reverse_ranks = False
            self.verbose = True

    class MockParser:
        def add_argument(self, *args, **kwargs):
            pass

        def parse_args(self):
            return MockArgs()

    def return_mock_parser(*args, **kwargs):
        return MockParser()

    monkeypatch.setattr(argparse, "ArgumentParser", return_mock_parser, raising=True)
    recursive_dbscan.main()
    assert out.exists()
    df = pd.read_csv(out, sep="\t")
    assert "contig" in df.columns
    assert "cluster" in df.columns
