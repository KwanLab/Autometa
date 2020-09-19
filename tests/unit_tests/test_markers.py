import pytest
import pandas as pd
from Bio import SeqIO
from autometa.common import markers
from autometa.common.markers import MARKERS_DIR

# TODO: from autometa.datasets import download_test_data()


@pytest.fixture(name="orfs")
def fixture_orfs(variables, tmp_path):
    markers_test_data = variables["markers"]
    records = markers_test_data["orfs"]
    outlines = ""
    for record, seq in records.items():
        outlines += f"{record}\n{seq}\n"
    fpath = tmp_path / "orfs.faa"
    with open(fpath, "w") as fh:
        fh.write(outlines)
    return fpath.as_posix()


@pytest.fixture(name="bacteria_markers")
def fixture_bacteria_markers(variables, tmp_path):
    markers_test_data = variables["markers"]
    df = pd.read_json(markers_test_data["bacteria"])
    df.set_index("contig", inplace=True)
    out = tmp_path / "bacteria.markers.tsv"
    df.to_csv(out, sep="\t", index=True, header=True)
    return out.as_posix()


@pytest.fixture(name="archaea_markers")
def fixture_archaea_markers(variables, tmp_path):
    markers_test_data = variables["markers"]
    df = pd.read_json(markers_test_data["archaea"])
    df.set_index("contig", inplace=True)
    out = tmp_path / "archaea.markers.tsv"
    df.to_csv(out, sep="\t", index=True, header=True)
    return out.as_posix()


@pytest.fixture(name="markers_fpath")
def fixture_markers_fpath(request, bacteria_markers, archaea_markers):
    kingdom = request.param
    markers_fpaths = {
        "bacteria_markers": bacteria_markers,
        "archaea_markers": archaea_markers,
    }
    return markers_fpaths[kingdom]


@pytest.fixture(name="dbdir")
def fixture_dbdir():
    return MARKERS_DIR


# TODO: Monkeypatch calls to hmmer?
# Dependency on hmmpressed dbs seems like we need to try something different here now.
@pytest.mark.wip
@pytest.mark.parametrize("kingdom", ["bacteria", "archaea"])
def test_marker_get(kingdom, orfs, dbdir):
    df = markers.get(kingdom, orfs, dbdir)
    assert not df.empty


@pytest.mark.skip
@pytest.mark.wip
@pytest.mark.parametrize(
    "markers_fpath",
    [pytest.param("bacteria_markers"), pytest.param("archaea_markers")],
    indirect=True,
)
def test_markers_load(markers_fpath):
    df = markers.load(fpath=markers_fpath)
    assert df.index.name == "contig"
