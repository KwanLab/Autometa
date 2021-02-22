import pytest
import subprocess
import os
import argparse

import pandas as pd


from autometa.common import coverage
from autometa.common.external import samtools


@pytest.fixture(name="coverage_testdir", scope="module")
def fixture_coverage_testdir(tmp_path_factory):
    return tmp_path_factory.mktemp("coverage")


@pytest.fixture(name="metagenome", scope="module")
def fixture_metagenome(variables, coverage_testdir):
    metagenome_test_data = variables["metagenome"]
    records = metagenome_test_data["assembly"]
    outlines = ""
    for record, seq in records.items():
        outlines += f"{record}\n{seq}\n"
    fpath = coverage_testdir / "metagenome.fna"
    with open(fpath, "w") as fh:
        fh.write(outlines)
    return str(fpath)


@pytest.fixture(name="sam_alignment", scope="module")
def fixture_sam_alignment(variables, coverage_testdir):
    coverage_test_data = variables["coverage"]
    outlines = coverage_test_data["sam"]
    fpath = coverage_testdir / "records.sam"
    with open(fpath, "w") as fh:
        fh.write(outlines)
    return str(fpath)


@pytest.fixture(name="bam_alignment", scope="module")
def fixture_bam_alignment(sam_alignment, coverage_testdir):
    bam_fpath = coverage_testdir / "records.bam"
    samtools.sort(sam=sam_alignment, bam=bam_fpath)
    return str(bam_fpath)


@pytest.fixture(name="bed_alignment", scope="module")
def fixture_bed_alignment(variables, coverage_testdir):
    coverage_test_data = variables["coverage"]
    alignment_records = pd.read_json(coverage_test_data["bed"])
    fpath = coverage_testdir / "records.bed"
    alignment_records.to_csv(fpath, sep="\t", header=True, index=False)
    return str(fpath)


@pytest.fixture(name="fwd_reads", scope="module")
def fixture_fwd_reads(variables, coverage_testdir):
    coverage_test_data = variables["coverage"]
    outlines = coverage_test_data["fwd_reads"]
    fpath = coverage_testdir / "fwd_reads.fastq"
    with open(fpath, "w") as fh:
        fh.write(outlines)
    return str(fpath)


@pytest.fixture(name="rev_reads", scope="module")
def fixture_rev_reads(variables, coverage_testdir):
    coverage_test_data = variables["coverage"]
    outlines = coverage_test_data["rev_reads"]
    fpath = coverage_testdir / "rev_reads.fastq"
    with open(fpath, "w") as fh:
        fh.write(outlines)
    return str(fpath)


@pytest.fixture(name="df_exists_fpath")
def fixture_df_without_contig_index_(tmp_path):
    df_dict = {
        "contig": ["contig_1", "contig_2", "contig_3"],
        "coverage": [1, 2, 3],
    }
    df = pd.DataFrame(df_dict)
    df_fpath = tmp_path / "invalid_df.tsv"
    df.to_csv(df_fpath, sep="\t")
    return str(df_fpath)


def test_coverage_get_from_sam(metagenome, sam_alignment, tmp_path):
    out = tmp_path / "covs_from_sam.tsv"
    df = coverage.get(fasta=metagenome, from_spades=False, out=out, sam=sam_alignment)
    assert df.index.name == "contig"
    assert "coverage" in df.columns
    assert out.exists()


def test_coverage_get_from_bam(metagenome, bam_alignment, tmp_path):
    out = tmp_path / "covs_from_bam.tsv"
    df = coverage.get(fasta=metagenome, from_spades=False, out=out, bam=bam_alignment)
    assert df.index.name == "contig"
    assert "coverage" in df.columns
    assert out.exists()


def test_coverage_get_from_bed(metagenome, bed_alignment, tmp_path):
    out = tmp_path / "covs_from_bed.tsv"
    df = coverage.get(fasta=metagenome, from_spades=False, out=out, bed=bed_alignment)
    assert df.index.name == "contig"
    assert "coverage" in df.columns
    assert out.exists()


def test_coverage_get_from_reads(metagenome, fwd_reads, rev_reads, tmp_path):
    out = tmp_path / "covs_from_bed.tsv"
    df = coverage.get(
        fasta=metagenome,
        from_spades=False,
        out=out,
        fwd_reads=fwd_reads,
        rev_reads=rev_reads,
    )
    assert df.index.name == "contig"
    assert "coverage" in df.columns
    assert out.exists()


def test_coverage_get_from_spades(metagenome, tmp_path):
    out = tmp_path / "covs_from_spades.tsv"
    df = coverage.get(fasta=metagenome, from_spades=True, out=out)
    assert df.index.name == "contig"
    assert "coverage" in df.columns
    assert out.exists()


def test_get_not_enough_passed_arguments(metagenome, tmp_path):
    out = tmp_path / "covs.tsv"
    with pytest.raises(ValueError):
        coverage.get(fasta=metagenome, from_spades=False, out=out)


def test_embed_df_already_exists(metagenome, df_exists_fpath, bed_alignment):
    coverage.get(fasta=metagenome, out=df_exists_fpath, bed=bed_alignment)


@pytest.mark.entrypoint
@pytest.mark.parametrize("from_spades", [False, True])
def test_coverage_main(
    monkeypatch,
    tmp_path,
    metagenome,
    sam_alignment,
    bam_alignment,
    bed_alignment,
    fwd_reads,
    rev_reads,
    from_spades,
):
    out = tmp_path / "coverages.tsv"

    class MockArgs:
        def __init__(self):
            self.assembly = metagenome
            self.fwd_reads = fwd_reads
            self.rev_reads = rev_reads
            self.sam = sam_alignment
            self.bam = bam_alignment
            self.lengths = "lengths.tsv"
            self.bed = bed_alignment
            self.cpus = 2
            self.out = out
            self.from_spades = from_spades

    # Defining the MockParser class to represent parser
    class MockParser:
        def add_argument(self, *args, **kwargs):
            pass

        def parse_args(self):
            return MockArgs()

    def return_mock_parser(*args, **kwargs):
        return MockParser()

    monkeypatch.setattr(argparse, "ArgumentParser", return_mock_parser, raising=True)

    coverage.main()
    assert out.exists()
    df = pd.read_csv(out, sep="\t")
    assert "contig" in df.columns
    assert "coverage" in df.columns
