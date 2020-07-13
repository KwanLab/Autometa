#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import pandas as pd
import pytest

from autometa.common.metagenome import Metagenome


@pytest.fixture(scope="session")
def test_metagenome(tmpdir_factory):
    tmpdir = tmpdir_factory.mktemp("data")
    mg = Metagenome(
        assembly="tests/data/metagenome.fna",
        outdir=tmpdir,
        nucl_orfs_fpath=os.path.join(tmpdir, "orfs.fna"),
        prot_orfs_fpath=os.path.join(tmpdir, "orfs.faa"),
        taxonomy_fpath=os.path.join(tmpdir, "taxonomy.tsv"),
        fwd_reads=None,
        rev_reads=None,
        se_reads=None,
        taxon_method="majority_vote",
    )
    assert mg.largest_seq == "NODE_1_length_1389215_cov_225.275"
    assert mg.size == 79468311
    assert mg.nseqs == 3587
    return mg


def test_length_filter(test_metagenome):
    out = os.path.join(test_metagenome.outdir, "metagenome.filtered.fna")
    mg = test_metagenome.length_filter(out=out, cutoff=3000)
    assert mg.nseqs == 2059
    assert mg.size == 78034179
    with pytest.raises(FileExistsError):
        test_metagenome.length_filter(out=out, cutoff=3000)


def test_call_orfs(test_metagenome):
    # TODO: Add mock here for prodigal call
    assert test_metagenome.orfs_called is False
    nucls_fp, prots_fp = test_metagenome.call_orfs(cpus=4, parallel=True)
    assert nucls_fp == test_metagenome.nucl_orfs_fpath
    assert prots_fp == test_metagenome.prot_orfs_fpath
    assert test_metagenome.orfs_called
    assert len(test_metagenome.prots) == 75390
    assert len(test_metagenome.nucls) == 75390


def test_get_kmers(test_metagenome):
    out = os.path.join(test_metagenome.outdir, "kmers.tsv")
    kmers_df = test_metagenome.get_kmers(
        kmer_size=5, multiprocess=True, out=out, nproc=4
    )
    assert kmers_df.shape == (3587, 512)
    normalized = os.path.join(test_metagenome.outdir, "kmers.normalized.tsv")
    normalized_df = test_metagenome.get_kmers(
        kmer_size=5,
        multiprocess=True,
        out=out,
        nproc=4,
        normalized=normalized,
        force=True,
    )
    assert normalized_df.shape == (3587, 512)


def test_get_coverages(test_metagenome):
    out = os.path.join(test_metagenome.outdir, "spades_coverages.tsv")
    cov_df = test_metagenome.get_coverages(out=out, from_spades=True)
    assert isinstance(cov_df, pd.DataFrame)
    assert cov_df.shape == (3587, 1)
    assert "coverage" in cov_df.columns


def test_assign_taxonomy(test_metagenome):
    assert test_metagenome.taxonomy_assigned is False
    params = {
        "cpus": 1,
        "ncbi": "</path/to/mocked/ncbi/dir>",
        "verbose": False,
        "blast": "</path/to/mocked/blastp.tsv",
        "hits": "</path/to/mocked/hits.pkl.gz>",
        "force": True,
    }
    # taxa_df = test_metagenome.assign_taxonomy(force=False, kwargs=params)
    # assert taxa_df.shape == (3587, )
