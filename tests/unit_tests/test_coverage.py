import pytest

from autometa.common import coverage

# test_coverage_data.json
# Should be test_data.json keyed by file being tested
# e.g. coverage_test_data = variables["coverage"]


@pytest.fixture(name="small_metagenome")
def fixture_metagenome(variables, tmp_path):
    kmer_test_data = variables["kmers"]
    records = kmer_test_data["small_metagenome"]
    outlines = ""
    for record, seq in records.items():
        outlines += f"{record}\n{seq}\n"
    fpath = tmp_path / "small_metagenome.fna"
    with open(fpath, "w") as fh:
        fh.write(outlines)
    return fpath.as_posix()


# TODO: Create fixtures for each set of data.
# Then group fixtures into a group and indirectly parametrize group
# i.e. alignments.sam, alignments.bam, alignments.bed, fwd_reads.fq, rev_reads.fq
# The external tools could be run or we could monkeypatch these.


def test_coverage_get_from_spades(small_metagenome, tmp_path):
    out = tmp_path / "covs_from_spades.tsv"
    df = coverage.get(fasta=small_metagenome, from_spades=True, out=out)
    assert df.index.name == "contig"
    assert "coverage" in df.columns
    assert out.exists()


@pytest.mark.skip
@pytest.mark.wip
def test_coverage_get_from_sam(small_metagenome, tmp_path):
    out = tmp_path / "covs_from_sam.tsv"
    df = coverage.get(fasta=small_metagenome, from_spades=False, out=out)
    assert df.index.name == "contig"
    assert "coverage" in df.columns
    assert out.exists()


@pytest.mark.skip
@pytest.mark.wip
def test_coverage_get_from_bam(small_metagenome, tmp_path):
    out = tmp_path / "covs_from_bam.tsv"
    df = coverage.get(fasta=small_metagenome, from_spades=False, out=out)
    assert df.index.name == "contig"
    assert "coverage" in df.columns
    assert out.exists()
