import pytest
import pandas as pd

from autometa.binning import recursive_dbscan
from autometa.common.markers import load as load_markers
import argparse


@pytest.fixture(name="kmers")
def fixture_kmers(variables, tmp_path):
    binning_test_data = variables["binning"]
    df = pd.read_json(binning_test_data["kmers_normalized"])
    fpath = tmp_path / "kmers.norm.tsv"
    df.to_csv(fpath, sep="\t", index=False, header=True)
    return fpath


@pytest.fixture(name="embedded_kmers")
def fixture_embedded_kmers(variables, tmp_path):
    binning_test_data = variables["binning"]
    df = pd.read_json(binning_test_data["kmers_embedded"])
    fpath = tmp_path / "kmers.embed.tsv"
    df.to_csv(fpath, sep="\t", index=False, header=True)
    return fpath


@pytest.fixture(name="coverage")
def fixture_coverage(variables, tmp_path):
    binning_test_data = variables["binning"]
    df = pd.read_json(binning_test_data["coverage"])
    fpath = tmp_path / "coverage.tsv"
    df.to_csv(fpath, sep="\t", index=False, header=True)
    return fpath


@pytest.fixture(name="taxonomy")
def fixture_taxonomy(variables, tmp_path):
    binning_test_data = variables["binning"]
    df = pd.read_json(binning_test_data["taxonomy"])
    fpath = tmp_path / "taxonomy.tsv"
    df.to_csv(fpath, sep="\t", index=False, header=True)
    return fpath


@pytest.fixture(name="markers_fpath")
def fixture_markers_fpath(variables, tmp_path):
    binning_test_data = variables["binning"]
    df = pd.read_json(binning_test_data["markers"])
    fpath = tmp_path / "markers.tsv"
    df.to_csv(fpath, sep="\t", index=False, header=True)
    return fpath


@pytest.fixture(name="markers")
def fixture_markers(markers_fpath, tmp_path):
    return load_markers(markers_fpath)


@pytest.fixture(name="master")
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
    assert "cluster" in df.columns
    assert "purity" in df.columns
    assert "completeness" in df.columns
    assert df.shape[0] == num_contigs


def test_binning_invalid_clustering_method(master, markers):
    with pytest.raises(ValueError):
        df = recursive_dbscan.binning(
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


@pytest.fixture(name="mock_parser")
def fixture_mock_parser(
    monkeypatch, kmers, coverage, markers_fpath, embedded_kmers, taxonomy, tmp_path
):
    def return_mock_parser(*args, **kwargs):
        return MockParser()

    class MockParseArgs:
        def __init__(self, kmers, coverage, markers, embedded, taxonomy, out):
            self.domain = "bacteria"
            self.kmers = kmers
            self.coverage = coverage
            self.markers = markers_fpath
            self.out = out
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
            out = tmp_path / "binning.tsv"
            return MockParseArgs(
                kmers, coverage, markers_fpath, embedded_kmers, taxonomy, out
            )

    # Defining the MockParser class to represent parser
    monkeypatch.setattr(argparse, "ArgumentParser", return_mock_parser, raising=True)


def test_recursive_dbscan_main(monkeypatch, mock_parser):
    with monkeypatch.context() as m:

        def return_args(*args, **kwargs):
            assert not args
            assert isinstance(kwargs["master"], pd.DataFrame)
            assert isinstance(kwargs["markers"], pd.DataFrame)
            assert isinstance(kwargs["taxonomy"], bool)
            assert kwargs["taxonomy"] == True
            assert isinstance(kwargs["starting_rank"], str)
            assert kwargs["starting_rank"] == "superkingdom"
            assert isinstance(kwargs["reverse_ranks"], bool)
            assert kwargs["reverse_ranks"] == False
            assert isinstance(kwargs["domain"], str)
            assert kwargs["domain"] == "bacteria"
            assert isinstance(kwargs["completeness"], float)
            assert kwargs["completeness"] == 20.0
            assert isinstance(kwargs["purity"], float)
            assert kwargs["purity"] == 90.0
            assert isinstance(kwargs["method"], str)
            assert kwargs["method"] == "dbscan"
            assert isinstance(kwargs["verbose"], bool)
            assert kwargs["verbose"] == True
            return pd.DataFrame(columns=["cluster", "completeness", "purity"])

        m.setattr(recursive_dbscan, "binning", return_args, raising=True)
        recursive_dbscan.main()
