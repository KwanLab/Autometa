import pytest
import pandas as pd
from Bio import SeqIO
from autometa.common import kmers

# TODO: from autometa.datasets import download_test_data()


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


@pytest.fixture(name="counts")
def fixture_counts(variables):
    kmer_test_data = variables["kmers"]
    df = pd.read_json(kmer_test_data["counts"])
    df.set_index("contig", inplace=True)
    return df


@pytest.fixture(name="counts_fpath")
def fixture_counts_fpath(counts, tmp_path):
    fpath = tmp_path / "counts.tsv"
    counts.to_csv(fpath, sep="\t", index=True, header=True)
    return fpath.as_posix()


@pytest.fixture(name="norm_df")
def fixture_norm_df(variables):
    kmer_test_data = variables["kmers"]
    df = pd.read_json(kmer_test_data["norm_df"])
    df.set_index("contig", inplace=True)
    return df


def test_kmer_load(counts_fpath):
    df = kmers.load(kmers_fpath=counts_fpath)
    assert not df.empty
    assert df.index.name == "contig"


@pytest.mark.parametrize("multiprocess", [True, False])
def test_count(small_metagenome, multiprocess, tmp_path):
    out = tmp_path / "kmers.tsv"
    size = 5
    force = False
    df = kmers.count(
        assembly=small_metagenome,
        size=size,
        out=out,
        force=force,
        multiprocess=multiprocess,
    )
    assert df.shape[1] == 4 ** size / 2
    assert df.index.name == "contig"
    assert out.exists()


@pytest.mark.parametrize("force", [True, False])
def test_count_out_exists(small_metagenome, counts, force, tmp_path):
    out = tmp_path / "kmers.tsv"
    counts.to_csv(out, sep="\t", index=True, header=True)
    size = 5
    df = kmers.count(
        assembly=small_metagenome, size=size, out=out, force=force, multiprocess=True,
    )
    assert df.shape[1] == 4 ** size / 2
    assert df.index.name == "contig"
    assert out.exists()


def test_count_wrong_size(small_metagenome, tmp_path):
    out = tmp_path / "kmers.tsv"
    size = 5.5
    with pytest.raises(TypeError):
        df = kmers.count(assembly=small_metagenome, size=size)


@pytest.mark.parametrize("method", ["am_clr", "clr", "ilr"])
def test_normalize(counts, method, tmp_path):
    out = tmp_path / "kmers.norm.tsv"
    force = False
    df = kmers.normalize(df=counts, method=method, out=out, force=force)
    if method in {"am_clr", "clr"}:
        assert df.shape == counts.shape
    else:
        # ILR will reduce the columns by one.
        assert df.shape[1] < counts.shape[1]
    assert out.exists()


@pytest.mark.parametrize("force", [True, False])
def test_normalize_out_exists(counts, norm_df, force, tmp_path):
    out = tmp_path / "kmers.norm.tsv"
    norm_df.to_csv(out, sep="\t", index=True, header=True)
    df = kmers.normalize(df=counts, method="am_clr", out=out, force=force)
    assert df.shape == counts.shape
    assert df.index.name == "contig"


def test_normalize_wrong_method(counts, tmp_path):
    out = tmp_path / "kmers.norm.tsv"
    with pytest.raises(ValueError):
        df = kmers.normalize(df=counts, method="am_ilr", out=out, force=False)


@pytest.mark.parametrize("method", ["bhsne", "sksne", "umap"])
def test_embed(norm_df, method, tmp_path):
    seed = 42
    out = tmp_path / "kmers.embed.tsv"
    force = False
    embed_dimensions = 2
    do_pca = True
    pca_dimensions = 50
    df = kmers.embed(
        kmers=norm_df,
        out=out,
        force=force,
        embed_dimensions=embed_dimensions,
        do_pca=do_pca,
        pca_dimensions=pca_dimensions,
        method=method,
        seed=seed,
    )
    assert df.shape[1] == embed_dimensions


def test_embed_out_exists(norm_df, tmp_path):
    seed = 42
    out = tmp_path / "kmers.embed.tsv"
    force = False
    method = "bhsne"
    embed_dimensions = 2
    do_pca = True
    pca_dimensions = 50
    df = kmers.embed(
        kmers=norm_df,
        out=out,
        force=force,
        embed_dimensions=embed_dimensions,
        do_pca=do_pca,
        pca_dimensions=pca_dimensions,
        method=method,
        seed=seed,
    )
    assert df.shape[1] == embed_dimensions
    df = kmers.embed(
        kmers=norm_df,
        out=out,
        force=force,
        embed_dimensions=embed_dimensions,
        do_pca=do_pca,
        pca_dimensions=pca_dimensions,
        method=method,
        seed=seed,
    )
