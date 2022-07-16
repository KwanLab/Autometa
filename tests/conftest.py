import pytest

from autometa.taxonomy.ncbi import TAXA_DB


@pytest.fixture(name="ncbi", scope="session", autouse=True)
def fixture_ncbi(taxa_db_dir):
    return TAXA_DB(dirpath=taxa_db_dir, verbose=False)


@pytest.fixture(name="dbdir", scope="session", autouse=True)
def fixture_dbdir(tmp_path_factory):
    """temporary TAXA_DB database directory

    Parameters
    ----------
    tmp_path_factory : Pathlib.LocalPath
        pytest generated temporary directory

    Returns
    -------
    Pathlib.LocalPath
        Path to temporary TAXA_DB database directory
    """
    tmpdir = tmp_path_factory.mktemp("databases")
    return tmpdir


@pytest.fixture(name="nodes", scope="session", autouse=True)
def fixture_nodes(variables):
    """Lines corresponding to nodes.dmp"""
    vote_test_data = variables["taxonomy"]
    lines = ""
    for child, info in vote_test_data["nodes"].items():
        # Parsing in ncbi.parse_nodes() follows TAXA_DB nodes.dmp file format.
        # child, parent, rank = line.split("\t|\t")[:3]
        parent = info.get("parent")
        rank = info.get("rank")
        child, parent, rank = [str(i) for i in [child, parent, rank]]
        lines += "\t|\t".join([child, parent, rank, "null"]) + "\n"
    return lines


@pytest.fixture(name="names", scope="session", autouse=True)
def fixture_names(variables):
    vote_test_data = variables["taxonomy"]
    lines = ""
    for taxid, name in vote_test_data["names"].items():
        # Parsing in ncbi.parse_names() follows TAXA_DB names.dmp file format.
        # taxid, name, __, classification = line.strip("\t|\n").split("\t|\t")[:4]
        # scientific name is checked with classification variable.
        lines += (
            "\t|\t".join([taxid, name, "null", "scientific name", "null"]) + "\t|\n"
        )
        # lines follows structure of names.dmp database.
    return lines


@pytest.fixture(name="merged", scope="session", autouse=True)
def fixture_merged(variables):
    vote_test_data = variables["taxonomy"]
    # Parsing in ncbi.parse_merged() follows TAXA_DB merged.dmp file format.
    lines = ""
    for old_taxid, new_taxid in vote_test_data["merged"].items():
        # old_taxid, new_taxid = line.strip("\t|\n").split("\t|\t")
        lines += f"{old_taxid}\t|\t{new_taxid}\t|\n"
    return lines


@pytest.fixture(name="delnodes", scope="session", autouse=True)
def fixture_delnodes(variables):
    # Parsing in ncbi.parse_delnodes() follows TAXA_DB delnodes.dmp file format.
    deltaxids = [
        2564079,
        2564078,
        2564077,
        2564076,
        2564075,
        2564074,
        2564073,
        2564072,
        2564071,
        2564070,
    ]
    lines = "\t|\n".join([str(taxid) for taxid in deltaxids]) + "\t|\n"
    return lines


@pytest.fixture(name="acc2taxid", scope="session", autouse=True)
def fixture_acc2taxid(variables):
    vote_test_data = variables["taxonomy"]
    # Parsing in diamond.add_taxids(...) follows TAXA_DB prot.accession2taxid file format.
    lines = "accession\taccession.version\ttaxid\tnull\n"
    for acc_num, taxid in vote_test_data["acc2taxid"].items():
        # acc_num, acc_ver, taxid, _ = line.split("\t")
        lines += f"{acc_num}\tnull\t{taxid}\tnull\n"
    return lines


@pytest.fixture(name="taxa_db_dir", scope="session", autouse=True)
def fixture_taxa_db_dir(dbdir, merged, nodes, names, delnodes, acc2taxid):
    taxa_db_dir = dbdir.joinpath("ncbi")
    taxa_db_dir.mkdir()
    # NOTE: TAXA_DB instance expects file name prot.accession2taxid
    acc2taxid_fpath = taxa_db_dir / "prot.accession2taxid"
    # NOTE: TAXA_DB instance expects file name merged.dmp
    merged_fpath = taxa_db_dir / "merged.dmp"
    # NOTE: TAXA_DB instance expects file name names.dmp
    names_fpath = taxa_db_dir / "names.dmp"
    # NOTE: TAXA_DB instance expects file name nodes.dmp
    nodes_fpath = taxa_db_dir / "nodes.dmp"
    # NOTE: TAXA_DB instance expects file name delnodes.dmp
    delnodes_fpath = taxa_db_dir / "delnodes.dmp"
    for fpath, data in zip(
        [merged_fpath, nodes_fpath, names_fpath, delnodes_fpath, acc2taxid_fpath],
        [merged, nodes, names, delnodes, acc2taxid],
    ):
        with open(fpath, "w") as fh:
            fh.write(data)
    return taxa_db_dir
