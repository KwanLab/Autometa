import pytest
import argparse
import pandas as pd
from autometa.common import markers
from autometa.common.external import hmmer
from autometa.common.markers import MARKERS_DIR


@pytest.fixture(name="orfs")
def fixture_orfs(variables, tmp_path):
    kmer_test_data = variables["markers"]
    records = kmer_test_data["orfs"]
    outlines = ""
    for record, seq in records.items():
        outlines += f"{record}\n{seq}\n"
    fpath = tmp_path / "orfs.faa"
    with open(fpath, "w") as fh:
        fh.write(outlines)
    return fpath.as_posix()


@pytest.fixture(name="scans_fpath")
def fixture_scans_fpath(variables, tmp_path):
    markers_test_data = variables["markers"]
    df = pd.read_json(markers_test_data["scans"])
    # There should be 24 columns in the hmmscan output file.
    assert df.shape[1] == 24
    out = tmp_path / "hmmscan_output.tsv"
    df.to_csv(out, sep="\t", index=False, header=True)
    return out.as_posix()


@pytest.fixture(name="filtered_markers_fpath")
def fixture_filtered_markers_fpath(variables, tmp_path):
    markers_test_data = variables["markers"]
    df = pd.read_json(markers_test_data["filtered_markers"])
    # There should be 6 columns in the hmmscan output file.
    assert df.shape[1] == 6
    out = tmp_path / "filtered_markers.tsv"
    df.to_csv(out, sep="\t", index=False, header=True)
    return out


@pytest.fixture(name="markers_dbdir")
def fixture_markers_dbdir():
    return MARKERS_DIR


@pytest.fixture(name="mock_hmmer_hmmscan")
def mock_hmmer_hmmscan(monkeypatch, scans_fpath):
    def return_args(*args, **kwargs):
        return scans_fpath

    monkeypatch.setattr(hmmer, name="hmmscan", value=return_args, raising=True)


@pytest.fixture(name="mock_hmmer_filter_markers")
def mock_hmmer_filter_markers(monkeypatch, filtered_markers_fpath):
    def return_args(*args, **kwargs):
        return filtered_markers_fpath

    monkeypatch.setattr(hmmer, name="filter_markers", value=return_args, raising=True)


def test_marker_get(
    markers_dbdir, orfs, mock_hmmer_hmmscan, mock_hmmer_filter_markers, tmp_path
):
    kingdom = "archaea"
    markers_fpath = tmp_path / f"{kingdom}.markers.tsv"
    df = markers.get(
        kingdom=kingdom,
        orfs=orfs,
        dbdir=markers_dbdir,
        out=markers_fpath,
        force=False,
        seed=42,
        format="wide",
    )
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    test_data_shape = (1, 1)
    assert df.shape == test_data_shape
    assert df.index.name == "contig"


@pytest.mark.parametrize(
    "format", ["wide", "long", "list", "counts", "invalid_format", "empty_fpath"]
)
def test_markers_load_formats(filtered_markers_fpath, format):
    if format == "invalid_format":
        with pytest.raises(ValueError):
            loaded_markers = markers.load(fpath=filtered_markers_fpath, format=format)
    elif format == "empty_fpath":
        with pytest.raises(FileNotFoundError):
            loaded_markers = markers.load(fpath="</empty/file/path>", format="wide")
    else:
        loaded_markers = markers.load(fpath=filtered_markers_fpath, format=format)

    if format == "wide":
        assert isinstance(loaded_markers, pd.DataFrame)
        assert loaded_markers.index.name == "contig"
        assert loaded_markers.shape == (1, 1)
    if format == "long":
        assert isinstance(loaded_markers, pd.DataFrame)
        assert loaded_markers.index.name == "contig"
        assert loaded_markers.shape == (1, 2)
        assert loaded_markers.columns.tolist() == ["sacc", "count"]
    if format == "list":
        assert isinstance(loaded_markers, dict)
        assert loaded_markers
        contig_markers_tuple = loaded_markers.popitem()
        # dict should be key="contig" and value=[marker,marker,...]
        assert isinstance(contig_markers_tuple[0], str)
        assert isinstance(contig_markers_tuple[1], list)
    if format == "counts":
        assert isinstance(loaded_markers, dict)
        assert loaded_markers
        contig_markers_tuple = loaded_markers.popitem()
        # dict should be key: str = "contig" and value: int = count
        assert isinstance(contig_markers_tuple[0], str)
        assert isinstance(contig_markers_tuple[1], int)


@pytest.fixture(name="mock_parser")
def fixture_mock_parser(monkeypatch):
    def return_mock_parser(*args, **kwargs):
        return MockParser()

    class MockParseArgs:
        kingdom = "bacteria"
        orfs = "orfs.faa"
        dbdir = "autometa/databases/markers"
        hmmscan = "bacteria.hmmscan.tsv"
        out = "bacteria.markers.tsv"
        force = True
        cpus = 1
        parallel = False
        gnu_parallel = False
        seed = 42

    class MockParser:
        def add_argument(self, *args, **kwargs):
            pass

        def parse_args(self):
            return MockParseArgs()

    # Defining the MockParser class to represent parser
    monkeypatch.setattr(argparse, "ArgumentParser", return_mock_parser, raising=True)


@pytest.mark.entrypoint
def test_markers_main(monkeypatch, mock_parser):
    with monkeypatch.context() as m:

        def return_args(*args, **kwargs):
            assert not args
            assert kwargs["kingdom"] == "bacteria"
            assert kwargs["orfs"] == "orfs.faa"
            assert kwargs["dbdir"] == "autometa/databases/markers"
            assert kwargs["scans"] == "bacteria.hmmscan.tsv"
            assert kwargs["out"] == "bacteria.markers.tsv"
            assert kwargs["force"] == True
            assert kwargs["seed"] == 42

        m.setattr(markers, "get", return_args, raising=True)
        markers.main()
