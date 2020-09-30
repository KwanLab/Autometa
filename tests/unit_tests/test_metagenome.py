import pytest

from autometa.common.metagenome import Metagenome
from autometa.common.external import prodigal


@pytest.fixture(name="assembly")
def fixture_assembly(variables, tmp_path):
    metagenome_test_data = variables["metagenome"]
    records = metagenome_test_data["assembly"]
    lines = ""
    for record, seq in records.items():
        lines += f"{record}\n{seq}\n"
    fpath = tmp_path / "assembly.fna"
    with open(fpath, "w") as fh:
        fh.write(lines)
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


@pytest.mark.parametrize("quality_measure", [0.1, 0.9])
def test_fragmentation_metric(metagenome, quality_measure):
    metric = metagenome.fragmentation_metric(quality_measure=quality_measure)
    if quality_measure == 0.1:
        expected = 7237
    elif metagenome.nseqs == 3 and quality_measure == 0.9:
        expected = 7231
    elif metagenome.nseqs == 4 and quality_measure == 0.9:
        expected = 7229
    else:
        expected = None

    assert isinstance(metric, int)
    if expected:
        assert metric == expected


def test_fragmentation_metric_wrong_quality_measure_value(metagenome):
    with pytest.raises(ValueError):
        metagenome.fragmentation_metric(quality_measure=1.4)


def test_length_weighted_gc(metagenome):
    if metagenome.nseqs == 3:
        expected = 65.3167472932504
    elif metagenome.nseqs == 4:
        expected = 63.195548489666145
    else:
        expected = None
    assert isinstance(metagenome.length_weighted_gc, float)
    if expected:
        assert metagenome.length_weighted_gc == expected


def test_largest_seq(metagenome):
    expected = "NODE_1501_length_7237_cov_222.15"
    assert metagenome.largest_seq == expected


def test_size(metagenome):
    if metagenome.nseqs == 3:
        expected = 21705
    elif metagenome.nseqs == 4:
        expected = 28934
    else:
        expected = None
    assert isinstance(metagenome.size, int)
    if expected:
        assert metagenome.size == expected


@pytest.mark.parametrize("details", [True, False])
def test_describe(metagenome, details):
    metagenome.describe(autometa_details=details)


@pytest.mark.parametrize("force", [False, True])
def test_force_length_filter(metagenome, force):
    existing = str(metagenome)
    if not force:
        with pytest.raises(FileExistsError):
            filtered_metagenome = metagenome.length_filter(
                out=existing, cutoff=7232, force=force
            )
    else:
        filtered_metagenome = metagenome.length_filter(
            out=existing, cutoff=7232, force=True
        )
        assert filtered_metagenome.nseqs == 2


@pytest.mark.parametrize("cutoff", [-6, "12324", 7232, 10000])
def test_length_filter_invalid_cutoff(metagenome, cutoff, tmp_path):
    out = tmp_path / "length_filtered.fna"
    expected_errors = {-6: ValueError, "12324": TypeError}
    if cutoff in expected_errors:
        expected_error = expected_errors[cutoff]
        with pytest.raises(expected_error):
            filtered_metagenome = metagenome.length_filter(
                out=out, cutoff=cutoff, force=False
            )
    else:
        filtered_metagenome = metagenome.length_filter(
            out=out, cutoff=cutoff, force=False
        )
    if cutoff == 7232:
        assert filtered_metagenome.nseqs == 2
    elif cutoff == 10000:
        # No contigs above cutoff. Ensure we return original
        assert filtered_metagenome.nseqs == metagenome.nseqs


def test_call_orfs(metagenome, mock_prodigal_run):
    force = False
    cpus = 1
    parallel = False
    args, kwargs = metagenome.call_orfs(force=force, cpus=cpus, parallel=parallel)
    # Let's ensure we are always passing in keyword arguments
    assert not args
    # Now let's make sure the arguments provided are the
    # same as were called in prodigal.run(...)
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
