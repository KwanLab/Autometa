import argparse
from multiprocessing import Value
from random import seed
from autometa.binning.unclustered_recruitment import Labels
import pandas as pd
import pytest

from autometa.binning import unclustered_recruitment
from autometa.common.markers import load as load_markers


@pytest.fixture(name="recruitment_testdir", scope="module")
def fixture_recruitment_testdir(tmp_path_factory):
    return tmp_path_factory.mktemp("recruitment")


# COMBAK: These are being used in binning, recruitment (and possibly summary)...
# So these should probably be session scoped....
@pytest.fixture(name="kmers", scope="module")
def fixture_kmers(variables, recruitment_testdir):
    binning_test_data = variables["binning"]
    df = pd.read_json(binning_test_data["kmers_normalized"])
    fpath = recruitment_testdir / "kmers.norm.tsv"
    df.to_csv(fpath, sep="\t", index=False, header=True)
    return fpath


@pytest.fixture(name="embedded_kmers", scope="module")
def fixture_embedded_kmers(variables, recruitment_testdir):
    binning_test_data = variables["binning"]
    df = pd.read_json(binning_test_data["kmers_embedded"])
    fpath = recruitment_testdir / "kmers.embed.tsv"
    df.to_csv(fpath, sep="\t", index=False, header=True)
    return fpath


@pytest.fixture(name="coverage", scope="module")
def fixture_coverage(variables, recruitment_testdir):
    binning_test_data = variables["binning"]
    df = pd.read_json(binning_test_data["coverage"])
    fpath = recruitment_testdir / "coverage.tsv"
    df.to_csv(fpath, sep="\t", index=False, header=True)
    return fpath


@pytest.fixture(name="taxonomy", scope="module")
def fixture_taxonomy(variables, recruitment_testdir):
    binning_test_data = variables["binning"]
    df = pd.read_json(binning_test_data["taxonomy"])
    fpath = recruitment_testdir / "taxonomy.tsv"
    df.to_csv(fpath, sep="\t", index=False, header=True)
    return fpath


@pytest.fixture(name="markers_fpath", scope="module")
def fixture_markers_fpath(variables, recruitment_testdir):
    binning_test_data = variables["binning"]
    df = pd.read_json(binning_test_data["markers"])
    fpath = recruitment_testdir / "markers.tsv"
    df.to_csv(fpath, sep="\t", index=False, header=True)
    return fpath


@pytest.fixture(name="assignments_fpath", scope="module")
def fixture_assignments_fpath(assignments, recruitment_testdir):
    fpath = recruitment_testdir / "assignments.tsv"
    assignments.to_csv(fpath, sep="\t", index=True, header=True)
    return fpath


@pytest.fixture(name="assignments", scope="module")
def fixture_assignments(variables):
    recruitment_test_data = variables["recruitment"]
    df = pd.read_json(recruitment_test_data["binning"])
    df.set_index("contig", inplace=True)
    return df


@pytest.fixture(name="markers", scope="module")
def fixture_markers(markers_fpath):
    return load_markers(markers_fpath)


@pytest.mark.wip
@pytest.mark.parametrize("dimensions", [None, 0, 50, 100, 1000])
def test_get_taxa_features(taxonomy, dimensions):
    if dimensions and dimensions > 100:
        with pytest.raises(ValueError):
            df = unclustered_recruitment.get_taxa_features(
                taxonomy, dimensions=dimensions
            )
    else:
        df = unclustered_recruitment.get_taxa_features(taxonomy, dimensions=dimensions)
        assert isinstance(df, pd.DataFrame)
        if dimensions:
            assert df.shape[1] == dimensions


@pytest.mark.parametrize("dimensions", [None, 0, 50, 100, 1000])
def test_get_kmer_features(kmers, dimensions):
    if dimensions and dimensions > 100:
        with pytest.raises(ValueError):
            df = unclustered_recruitment.get_kmer_features(
                filepath=kmers, dimensions=dimensions
            )
    else:
        df = unclustered_recruitment.get_kmer_features(
            filepath=kmers, dimensions=dimensions
        )
        assert isinstance(df, pd.DataFrame)
        if dimensions:
            assert df.shape[1] == dimensions


def test_get_features_no_taxa(
    kmers, coverage,
):
    df = unclustered_recruitment.get_features(kmers=kmers, coverage=coverage)
    assert isinstance(df, pd.DataFrame)


def test_get_features_with_taxa(kmers, coverage, taxonomy):
    df = unclustered_recruitment.get_features(
        kmers=kmers, coverage=coverage, taxonomy=taxonomy
    )


def test_get_labels(assignments):
    labels = unclustered_recruitment.get_labels(assignments)
    assert isinstance(labels, Labels)


@pytest.mark.wip
@pytest.mark.skip
@pytest.mark.parametrize("confidence", [0.1, 0.5, 1.0])
@pytest.mark.parametrize("classifier", ["decision_tree", "random_forest"])
def test_get_confidence_filtered_predictions(
    confidence, classifier, train_data, test_data
):
    df = unclustered_recruitment.get_confidence_filtered_predictions(
        train_data=train_data,
        test_data=test_data,
        num_classifications=2,
        confidence=confidence,
        classifier=classifier,
        seed=10,
    )
    assert isinstance(df, pd.DataFrame)


@pytest.mark.wip
@pytest.mark.skip
def test_filter_contaminating_predictions(markers, assignments):
    preds = pd.DataFrame()
    df = unclustered_recruitment.filter_contaminating_predictions(
        predictions=preds, markers=markers, binning=assignments
    )
    assert isinstance(df, pd.DataFrame)


@pytest.mark.wip
@pytest.mark.skip
def test_add_predictions(assignments):
    preds = pd.DataFrame()
    df = unclustered_recruitment.add_predictions(binning=assignments, predictions=preds)
    assert isinstance(df, pd.DataFrame)


@pytest.mark.entrypoint
def test_main(
    monkeypatch, kmers, coverage, taxonomy, assignments_fpath, markers_fpath, tmp_path
):
    out = tmp_path / "unclustered_recruitment.tsv"

    class MockArgs:
        def __init__(self):
            self.kmers: str = kmers
            self.coverage: str = coverage
            self.binning: str = assignments_fpath
            self.markers: str = markers_fpath
            self.output: str = out
            self.additional_features = []
            self.taxonomy: str = taxonomy
            self.kmer_dimensions: int = 50
            self.taxa_dimensions: int = 0
            self.seed: int = 42
            self.confidence: float = 1.0
            self.num_classifications: int = 10
            self.classifier: str = "decision_tree"

    class MockParser:
        def add_argument(self, *args, **kwds):
            pass

        def parse_args(self):
            return MockArgs()

    def return_mock_parser(*args, **kwds):
        return MockParser()

    monkeypatch.setattr(argparse, "ArgumentParser", return_mock_parser, raising=True)
    unclustered_recruitment.main()
