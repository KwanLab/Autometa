import pytest
import pandas as pd
from autometa.binning import summary
from autometa.taxonomy import majority_vote
from argparse import Namespace
from Bio import SeqIO
import argparse
from autometa.config import utilities
import glob


@pytest.fixture(name="bin_df", scope="module")
def fixture_bin_df(variables):
    binning_test_data = variables["summary"]
    df = pd.read_json(binning_test_data["bin_df"])
    df.set_index("contig", inplace=True)
    return df


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


@pytest.fixture(name="lengths_fpath", scope="module")
def fixture_lengths_fpath(assembly, tmp_path_factory):
    # Write lengths file for placement of filepath in Namespace
    fpath = tmp_path_factory.mktemp("summary") / "lengths.tsv"
    df = pd.DataFrame(
        [
            {"contig": rec.id, "length": int(rec.id.split("_")[3])}
            for rec in SeqIO.parse(assembly, "fasta")
        ]
    )
    df.to_csv(fpath, sep="\t", index=False, header=True)
    return fpath


@pytest.fixture(name="mgargs", scope="module")
def mock_mgargs(
    variables, tmp_path_factory, assembly, markers, bin_df, lengths_fpath, ncbi_dir
):
    binning_test_data = variables["binning"]
    tmpdir = tmp_path_factory.mktemp("mgargs")
    # Write coverages file for placement of filepath in Namespace
    coverages_fpath = tmpdir / "coverage.tsv"
    df = pd.read_json(binning_test_data["coverage"])
    df.to_csv(coverages_fpath, sep="\t", index=False, header=True)
    # Write taxonomy file for placement of filepath in Namespace
    taxonomy_fpath = tmpdir / "taxonomy.tsv"
    df = pd.read_json(binning_test_data["taxonomy"])
    df.to_csv(taxonomy_fpath, sep="\t", index=False, header=True)
    # Write bacteria binning file for placement of filepath in Namespace
    bacteria_binning = tmpdir / "bacteria.binning.tsv"
    outcols = ["cluster", "completeness", "purity"]
    bin_df[outcols].to_csv(bacteria_binning, sep="\t", index=True, header=True)
    # Write rest as empty files for placement of filepaths in Namespace
    bacteria_markers = tmpdir / "bacteria.markers.tsv"
    markers.to_csv(bacteria_markers, sep="\t", index=False, header=True)
    archaea_binning = tmpdir / "archaea.binning.tsv"
    archaea_markers = tmpdir / "archaea.markers.tsv"
    files = Namespace(
        coverages=coverages_fpath,
        taxonomy=taxonomy_fpath,
        bacteria_binning=bacteria_binning,
        bacteria_markers=bacteria_markers,
        archaea_binning=archaea_binning,
        archaea_markers=archaea_markers,
        metagenome=assembly,
        lengths=lengths_fpath,
    )
    databases = Namespace(ncbi=ncbi_dir)
    parameters = Namespace(outdir=tmpdir)
    return Namespace(files=files, databases=databases, parameters=parameters)


@pytest.mark.skip
def test_merge_annotations(mgargs):
    annotations = summary.merge_annotations(mgargs)
    for domain in ["bacteria", "archaea"]:
        assert domain in annotations
        df = annotations[domain]
        assert df.index.name == "contig"
        if domain == "archaea":
            assert df.empty
        else:
            for annotation in ["coverage", "length", "GC"]:
                assert annotation in df.columns


def test_merge_annotations_lengths_file_not_found(mgargs):
    mgargs.files.lengths = "empty_length_file"
    with pytest.raises(FileNotFoundError):
        annotations = summary.merge_annotations(mgargs)


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
def test_get_metabin_stats(bin_df, markers_fpath, assembly):
    df = summary.get_metabin_stats(
        bin_df=bin_df, markers_fpath=markers_fpath, assembly=assembly
    )
    assert df.index.name == "cluster"
    assert df.shape == (5, 20)


def test_get_metabin_stats_disjoint_markers(disjoint_markers_fpath, bin_df, assembly):
    df = summary.get_metabin_stats(
        bin_df=bin_df, markers_fpath=disjoint_markers_fpath, assembly=assembly
    )
    assert df.index.name == "cluster"
    assert df.shape == (5, 20)


@pytest.mark.skip
def test_write_cluster_records(bin_df, assembly, tmp_path):
    dirpath = tmp_path / "summary"
    summary.write_cluster_records(bin_df=bin_df, metagenome=assembly, outdir=dirpath)


@pytest.fixture(name="mock_parser")
def fixture_mock_parser(monkeypatch):
    def return_mock_parser(*args, **kwargs):
        return MockParser()

    class MockArgs:
        def __init__(self):
            self.workspace = "workspace"
            self.write = True

    class MockParser:
        def add_argument(self, *args, **kwargs):
            pass

        def parse_args(self):
            return MockArgs()

    # Defining the MockParser class to represent parser
    monkeypatch.setattr(argparse, "ArgumentParser", return_mock_parser, raising=True)


@pytest.mark.entrypoint
def test_main(monkeypatch, mock_parser, mgargs, lengths_fpath, mock_rank_taxids):
    mgargs.files.lengths = lengths_fpath
    with monkeypatch.context() as m:

        def return_mgargs(*args, **kwargs):
            return mgargs

        def return_config_fpaths(*args, **kwargs):
            return ["example/config/file/path"]

        m.setattr(glob, "glob", return_config_fpaths, raising=True)
        m.setattr(utilities, "parse_args", return_mgargs, raising=True)
        summary.main()
