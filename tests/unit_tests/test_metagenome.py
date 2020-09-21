import pytest

from autometa.common.metagenome import Metagenome
from autometa.common.external import prodigal


@pytest.fixture(name="assembly")
def fixture_assembly(variables, tmp_path):
    kmer_test_data = variables["kmers"]
    records = kmer_test_data["small_metagenome"]
    outlines = ""
    for record, seq in records.items():
        outlines += f"{record}\n{seq}\n"
    fpath = tmp_path / "small_metagenome.fna"
    with open(fpath, "w") as fh:
        fh.write(outlines)
    return fpath.as_posix()


@pytest.fixture(name="metagenome")
def fixture_metagenome(assembly, tmpdir):
    nucl_orfs = tmpdir.join("orfs.fna")
    prot_orfs = tmpdir.join("orfs.faa")
    return Metagenome(
        assembly=assembly,
        outdir=tmpdir,
        nucl_orfs_fpath=nucl_orfs,
        prot_orfs_fpath=prot_orfs,
    )


@pytest.fixture(name="mock_prodigal_run")
def mock_prodigal_run(monkeypatch):
    """Set the prodigal.run method to return filepaths."""

    def return_args(*args, **kwargs):
        return args, kwargs

    monkeypatch.setattr(prodigal, "run", return_args)


@pytest.mark.parametrize(
    "quality_measure,expected", [(0.5, 7004), (0.1, 7227), (0.9, 6764)]
)
def test_fragmentation_metric(metagenome, quality_measure, expected):
    metric = metagenome.fragmentation_metric(quality_measure=quality_measure)
    assert metric == expected


def test_fragmentation_metric_wrong_quality_measure_value(metagenome):
    with pytest.raises(ValueError):
        metagenome.fragmentation_metric(quality_measure=1.4)


def test_length_weighted_gc(metagenome):
    expected = 58.95893030432281
    assert metagenome.length_weighted_gc == expected


def test_largest_seq(metagenome):
    expected = "NODE_1501_length_7237_cov_222.15"
    assert metagenome.largest_seq == expected


@pytest.mark.parametrize("details", [True, False])
def test_describe(metagenome, details):
    metagenome.describe(autometa_details=details)


def test_length_filter(metagenome, tmp_path):
    out = tmp_path / "length_filtered.fna"
    filtered_metagenome = metagenome.length_filter(out=out, cutoff=7000, force=False)
    assert filtered_metagenome.nseqs == 25


def test_length_filter_filtered_exists(metagenome):
    existing = str(metagenome)
    with pytest.raises(FileExistsError):
        filtered_metagenome = metagenome.length_filter(
            out=existing, cutoff=7000, force=False
        )


def test_length_filter_overwrite(metagenome):
    existing = str(metagenome)
    filtered_metagenome = metagenome.length_filter(
        out=existing, cutoff=7000, force=True
    )
    assert filtered_metagenome.nseqs == 25


@pytest.mark.parametrize(
    "cutoff,expected_error", [(-6, ValueError), ("12324", TypeError)]
)
def test_length_filter_invalid_cutoff(metagenome, cutoff, expected_error, tmp_path):
    out = tmp_path / "length_filtered.fna"
    with pytest.raises(expected_error):
        filtered_metagenome = metagenome.length_filter(
            out=out, cutoff=cutoff, force=False
        )


def test_call_orfs(metagenome, mock_prodigal_run):
    force = False
    cpus = 1
    parallel = False
    args, kwargs = metagenome.call_orfs(force=force, cpus=cpus, parallel=parallel)
    assert metagenome.assembly == kwargs.get("assembly")
    assert metagenome.prot_orfs_fpath == kwargs.get("prots_out")
    assert metagenome.nucl_orfs_fpath == kwargs.get("nucls_out")
    assert force == kwargs.get("force")
    assert parallel == kwargs.get("parallel")
    assert cpus == kwargs.get("cpus")


@pytest.mark.parametrize(
    "force,parallel,cpus", [(1, False, 1), (False, 1, 1), (False, False, False)],
)
def test_call_orfs_invalid_params(metagenome, force, parallel, cpus):
    with pytest.raises(TypeError):
        metagenome.call_orfs(force=force, parallel=parallel, cpus=cpus)


def test_orfs_called(metagenome, monkeypatch):
    called = metagenome.orfs_called
    assert not called
    with monkeypatch.context():
        # Here we set the filepaths to the metagenome fixture assembly path
        # (since we know this exists and is non-empty)
        monkeypatch.setattr(metagenome, "prot_orfs_fpath", str(metagenome))
        monkeypatch.setattr(metagenome, "nucl_orfs_fpath", str(metagenome))
        called = metagenome.orfs_called
        assert called
