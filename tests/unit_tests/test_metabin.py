#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import pytest

import numpy as np

from autometa.common.metabin import MetaBin


ASSEMBLY = "tests/data/metagenome.fna"
CONTIGS = {
    "NODE_11_length_590705_cov_223.126",
    "NODE_12_length_584723_cov_225.634",
    "NODE_13_length_533917_cov_222.896",
    "NODE_100_length_142487_cov_223.46",
    "NODE_102_length_138307_cov_223.908",
    "NODE_105_length_136131_cov_224.464",
    "NODE_107_length_135043_cov_223.938",
    "NODE_108_length_134368_cov_224.716",
    "NODE_109_length_132210_cov_222.18",
    "NODE_110_length_132131_cov_223.987",
    "NODE_116_length_123186_cov_224.013",
    "NODE_119_length_119677_cov_223.947",
    "NODE_120_length_119298_cov_224.179",
    "NODE_132_length_112038_cov_223.935",
    "NODE_133_length_111688_cov_223.412",
    "NODE_136_length_111142_cov_224.73",
    "NODE_142_length_107924_cov_223.475",
    "NODE_146_length_105613_cov_224.716",
    "NODE_160_length_97942_cov_222.841",
    "NODE_1007_length_15092_cov_223.109",
    "NODE_1009_length_14996_cov_222.895",
    "NODE_1023_length_14648_cov_224.155",
    "NODE_1026_length_14554_cov_223.897",
    "NODE_1035_length_14277_cov_224.245",
    "NODE_1037_length_14262_cov_223.873",
    "NODE_1039_length_14197_cov_223.624",
    "NODE_1046_length_13906_cov_224.346",
    "NODE_1052_length_13792_cov_224.116",
    "NODE_1057_length_13739_cov_222.742",
    "NODE_1059_length_13673_cov_223.999",
    "NODE_1074_length_13345_cov_224.047",
    "NODE_1076_length_13288_cov_223.745",
    "NODE_1080_length_13246_cov_222.788",
    "NODE_1086_length_13138_cov_224.524",
    "NODE_1126_length_12393_cov_221.06",
    "NODE_1130_length_12333_cov_223.334",
    "NODE_1147_length_12002_cov_224.758",
    "NODE_1154_length_11835_cov_227.452",
    "NODE_1157_length_11817_cov_223.491",
    "NODE_1186_length_11353_cov_224.39",
    "NODE_1191_length_11292_cov_224.194",
    "NODE_1260_length_10179_cov_232.188",
    "NODE_1308_length_9538_cov_217.586",
    "NODE_20_length_436718_cov_222.09",
    "NODE_21_length_435208_cov_224.27",
    "NODE_22_length_427389_cov_223.136",
    "NODE_24_length_372190_cov_224.35",
    "NODE_26_length_356957_cov_224.046",
    "NODE_28_length_331107_cov_224.832",
    "NODE_30_length_320616_cov_227.122",
    "NODE_31_length_318595_cov_223.603",
}


@pytest.fixture(scope="session")
def test_metabin(tmpdir_factory):
    tmpdir = tmpdir_factory.mktemp("data")
    metabin = MetaBin(assembly=ASSEMBLY, contigs=CONTIGS, outdir=tmpdir)
    assert metabin.outdir == tmpdir
    assert metabin.basename == os.path.basename(ASSEMBLY)
    assert metabin.assembly == os.path.realpath(ASSEMBLY)
    return metabin


def test_describe(test_metabin):
    assert test_metabin.nseqs == len(CONTIGS)
    assert test_metabin.nallseqs == 3587
    assert test_metabin.seqs_pct == 1.4218009478672986
    assert test_metabin.size == 6980205
    assert test_metabin.totalsize == 79468311
    assert test_metabin.size_pct == 8.783633264836848
    assert test_metabin.length_weighted_gc == 49.5965089850513


def test_get_orfs(test_metabin):
    with pytest.raises(KeyError):
        test_metabin.get_orfs(orf_type="invalid_type")
    with pytest.raises(FileNotFoundError):
        test_metabin.get_orfs(orf_type="prot")


def test_write_orfs(test_metabin):
    prots_fp = os.path.join(test_metabin.outdir, "metagenome.orfs.faa")
    nucls_fp = os.path.join(test_metabin.outdir, "metagenome.orfs.fna")
    with pytest.raises(KeyError):
        test_metabin.write_orfs(prots_fp, orf_type="invalid_orf_type")
    with pytest.raises(FileNotFoundError):
        test_metabin.write_orfs(prots_fp, orf_type="nucl")
    with pytest.raises(FileNotFoundError):
        test_metabin.write_orfs(prots_fp, orf_type="prot")

    # test_metabin.write_orfs(nucls_fp, orf_type="nucl")


def test_markers(test_metabin):
    # bacteria_markers = test_metabin.markers(kingdom="bacteria", force=True)
    # archaea_markers = test_metabin.markers(kingdom="archaea", force=True)
    pass


def test_subset_df(test_metabin):
    pass


def test_orf_fpath_checking(test_metabin):
    nucls_fp = os.path.join(test_metabin.outdir, "metagenome.orfs.fna")
    prots_fp = os.path.join(test_metabin.outdir, "metagenome.orfs.faa")
    assert test_metabin.nucl_orfs_fpath == nucls_fp
    assert test_metabin.nucl_orfs_exist == False
    assert test_metabin.prepared(nucls_fp) == False
    assert test_metabin.prot_orfs_fpath == prots_fp
    assert test_metabin.prot_orfs_exist == False
    assert test_metabin.prepared(prots_fp) == False
    # assert test_metabin.nnucls == np.nan
    # assert test_metabin.nprots == np.nan
