import pytest

from autometa.taxonomy import vote
from autometa.taxonomy.ncbi import NCBI

import pandas as pd


@pytest.fixture(name="blastp")
def fixture_blastp(variables, tmp_path):
    vote_test_data = variables["taxonomy"]
    blastp_fpath = tmp_path / "blastp.tsv"
    df = pd.read_json(vote_test_data["blastp"])
    df.to_csv(blastp_fpath, sep="\t", index=False, header=True)
    return blastp_fpath.as_posix()


@pytest.fixture(name="prot_orfs")
def fixture_prot_orfs(variables, tmp_path):
    vote_test_data = variables["taxonomy"]
    blastp_fpath = tmp_path / "orfs.faa"
    df = pd.read_json(vote_test_data["prot_orfs"])
    df.to_csv(blastp_fpath, sep="\t", index=False, header=True)
    return blastp_fpath.as_posix()


@pytest.fixture(name="dbdir")
def fixture_dbdir(tmpdir):
    tmpdir.mkdir()
    return tmpdir


@pytest.fixture(name="nodes.dmp")
def fixture_nodes(variables, dbdir):
    vote_test_data = variables["taxonomy"]
    # NOTE: NCBI instance expects file name nodes.dmp
    nodes_fpath = dbdir / "nodes.dmp"
    lines = ""
    for child, info in vote_test_data["nodes"].items():
        # Parsing in ncbi.parse_nodes() follows NCBI nodes.dmp file format.
        # child, parent, rank = line.split("\t|\t")[:3]
        parent = info.get("parent")
        rank = info.get("rank")
        lines += "\t|\t".join([child, parent, rank, "null"]) + "\n"
    with open(nodes_fpath, "w") as fh:
        fh.write(lines)
    return nodes_fpath.as_posix()


@pytest.fixture(name="names.dmp")
def fixture_names(variables, dbdir):
    vote_test_data = variables["taxonomy"]
    # NOTE: NCBI instance expects file name names.dmp
    names_fpath = dbdir / "names.dmp"
    lines = ""
    for taxid, name in vote_test_data["names"].items():
        # Parsing in ncbi.parse_names() follows NCBI names.dmp file format.
        # taxid, name, __, classification = line.strip("\t|\n").split("\t|\t")[:4]
        lines += (
            "\t|\t".join([taxid, name, "null", "scientific name", "null"]) + "\t|\n"
        )
    with open(names_fpath, "w") as fh:
        fh.write(lines)
    return names_fpath.as_posix()


@pytest.fixture(name="merged.dmp")
def fixture_merged(variables, dbdir):
    vote_test_data = variables["taxonomy"]
    # NOTE: NCBI instance expects file name merged.dmp
    merged_fpath = dbdir / "merged.dmp"
    # Parsing in ncbi.parse_merged() follows NCBI merged.dmp file format.
    lines = ""
    for old_taxid, new_taxid in vote_test_data["merged"].items():
        # old_taxid, new_taxid = line.strip("\t|\n").split("\t|\t")
        lines += f"{old_taxid}\t|\t{new_taxid}\t|\n"
    with open(merged_fpath, "w") as fh:
        fh.write(lines)
    return merged_fpath.as_posix()


@pytest.fixture(name="acc2taxid.dmp")
def fixture_acc2taxid(variables, dbdir):
    vote_test_data = variables["taxonomy"]
    # NOTE: NCBI instance expects file name prot.accession2taxid
    acc2taxid_fpath = dbdir / "prot.accession2taxid"
    # Parsing in diamond.add_taxids(...) follows NCBI prot.accession2taxid file format.
    lines = ""
    for acc_num, info in vote_test_data["acc2taxid"].items():
        acc_ver = info.get("acc_ver")
        taxid = info.get("taxid")
        # acc_num, acc_ver, taxid, _ = line.split("\t")
        lines += f"{acc_num}\t{acc_ver}\t{taxid}\tnull\n"
    with open(acc2taxid_fpath, "w") as fh:
        fh.write(lines)
    return acc2taxid_fpath.as_posix()


@pytest.fixture(name="ncbi")
def fixture_ncbi(dbdir):
    dirpath = dbdir.as_posix()
    return NCBI(dirpath=dirpath, verbose=False)


@pytest.fixture(name="ncbi_dirpath")
def fixture_ncbi_dirpath(dbdir):
    return dbdir.as_posix()


def test_add_ranks(ncbi, monkeypatch):
    def mockncbi():
        return ncbi

    monkeypatch.setattr(vote.add_ranks, "ncbi", mockncbi)


@pytest.fixture(name="assigned_vote_fpath")
def fixture_assigned_vote_fpath(tmp_path):
    return tmp_path / "assigned_votes.tsv"


@pytest.mark.wip
def test_vote_assign(blastp, ncbi_dirpath, assigned_vote_fpath, prot_orfs, tmp_path):
    lca_fpath = tmp_path / "lca.tsv"
    vote.assign(
        outfpath=assigned_vote_fpath,
        prot_orfs=prot_orfs,
        blast=blastp,
        ncbi_dir=ncbi_dirpath,
        lca_fpath=lca_fpath,
    )


@pytest.mark.skip
@pytest.mark.wip
def test_vote_get(ncbi):
    fpath = ""
    assembly = ""
    ncbi_dir = ""
    kingdom = ""
    outdir = ""
    vote.get(
        fpath=fpath,
        assembly=assembly,
        kingdom=kingdom,
        ncbi_dir=ncbi_dir,
        outdir=outdir,
    )
