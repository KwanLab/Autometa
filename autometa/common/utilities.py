#!/usr/bin/env python
"""
COPYRIGHT
Copyright 2020 Ian J. Miller, Evan R. Rees, Kyle Wolf, Siddharth Uppal,
Shaurya Chanana, Izaak Miller, Jason C. Kwan

This file is part of Autometa.

Autometa is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Autometa is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with Autometa. If not, see <http://www.gnu.org/licenses/>.
COPYRIGHT

File containing common utilities functions to be used by Autometa scripts.
"""


import gzip
import hashlib
import logging
import os
import pickle
import sys
import tarfile
import time

import numpy as np
import pandas as pd

from datetime import datetime
from functools import wraps
from functools import namedtuple

from autometa.common.exceptions import TableFormatError

logger = logging.getLogger(__name__)


def unpickle(fpath):
    """Load a serialized `fpath` from :func:`make_pickle`.

    Parameters
    ----------
    fpath : str
        </path/to/file.pkl>.

    Returns
    -------
    any
        Python object that was serialized to file via make_pickle()

    Raises
    -------
    ExceptionName
        Why the exception is raised.

    """
    logger.debug(f"unpickling {fpath}")
    if fpath.endswith(".gz"):
        fh = gzip.open(fpath, "rb")
    else:
        fh = open(fpath, "rb")
    obj = pickle.load(file=fh)
    fh.close()
    logger.debug(f"{fpath} object unpickled")
    return obj


def make_pickle(obj, outfpath):
    """Serialize a python object (`obj`) to `outfpath`.
    Note:  Opposite of :func:`unpickle`

    Parameters
    ----------
    obj : any
        Python object to serialize to outfpath.
    outfpath : str
        </path/to/pickled/file>.

    Returns
    -------
    str
        </path/to/pickled/file.pkl>

    Raises
    -------
    ExceptionName
        Why the exception is raised.

    """
    logger.debug(f"pickling object to {outfpath}")
    if outfpath.endswith(".gz"):
        fh = gzip.open(outfpath, "wb")
    else:
        fh = open(outfpath, "wb")
    pickle.dump(obj=obj, file=fh)
    fh.close()
    logger.debug(f"object pickled to {outfpath}")
    return outfpath


def gunzip(infpath, outfpath, delete_original=False, block_size=65536):
    """Decompress gzipped `infpath` to `outfpath` and write checksum of `outfpath` upon successful decompression.

    Parameters
    ----------
    infpath : str
        </path/to/file.gz>
    outfpath : str
        </path/to/file>
    delete_original : bool
        Will delete the original file after successfully decompressing `infpath` (Default is False).
    block_size : int
        Amount of `infpath` to read in to memory before writing to `outfpath` (Default is 65536 bytes).

    Returns
    -------
    str
        </path/to/file>

    Raises
    -------
    FileExistsError
        `outfpath` already exists and is not empty

    """
    logger.debug(
        f"gunzipping {os.path.basename(infpath)} to {os.path.basename(outfpath)}"
    )
    if os.path.exists(outfpath) and os.path.getsize(outfpath) > 0:
        raise FileExistsError(outfpath)
    lines = ""
    with gzip.open(infpath, "rt") as fh, open(outfpath, "w") as out:
        for i, line in enumerate(fh):
            lines += line
            if sys.getsizeof(lines) >= block_size:
                out.write(lines)
                lines = ""
        out.write(lines)
    logger.debug(f"gunzipped {infpath} to {outfpath}")
    write_checksum(outfpath, f"{outfpath}.md5")
    if delete_original:
        os.remove(infpath)
        logger.debug(f"removed original file: {infpath}")
    return outfpath


def untar(tarchive, outdir, member=None):
    """Decompress a tar archive (may be gzipped or bzipped). passing in `member`
    requires an `outdir` also be provided.

    See: https://docs.python.org/3.8/library/tarfile.html#module-tarfile

    Parameters
    ----------
    tarchive : str
        </path/tarchive.tar.[compression]>
    outdir : str
        </path/to/output/directory>
    member : str, optional
        member file to extract.

    Returns
    -------
    str
        </path/to/extracted/member/file> if member else </path/to/output/directory>

    Raises
    -------
    IsADirectoryError
        `outdir` already exists
    ValueError
        `tarchive` is not a tar archive
    KeyError
        `member` was not found in `tarchive`

    """
    if not member and not outdir:
        raise ValueError(
            f"`member` or `outdir` must be passed: member={member} outdir={outdir}"
        )
    logger.debug(f"decompressing tarchive {tarchive} to {outdir}")
    outfpath = os.path.join(outdir, member) if member else None
    if member and os.path.exists(outfpath) and os.path.getsize(outfpath):
        raise FileExistsError(outfpath)
    if not tarfile.is_tarfile(tarchive):
        raise ValueError(f"{tarchive} is not a tar archive")
    if tarchive.endswith(".tar.gz") or tarchive.endswith(".tgz"):
        compression = "gz"
    elif tarchive.endswith(".tar.bz2"):
        compression = "bz2"
    else:
        compression = "*"
    with tarfile.open(tarchive, f"r:{compression}") as tar:
        if member:
            try:
                tar.extract(member=member, path=outdir)
            except KeyError as err:
                raise KeyError(f"member not in tarchive : {member} : {tarchive}")
        else:
            tar.extractall(outdir)
    if member:
        logger.debug(f"{member} extracted to {outdir}")
        return outfpath
    logger.debug(f"{tarchive} decompressed to {outdir}")
    return outdir


def tarchive_results(outfpath, src_dirpath):
    """Generate a tar archive of Autometa Results

    See:
    https://stackoverflow.com/questions/2032403/how-to-create-full-compressed-tar-file-using-python

    Parameters
    ----------
    outfpath : str
        </path/to/output/tar/archive.tar.gz || </path/to/output/tar/archive.tgz
    src_dirpath : str
        </paths/to/directory/to/archive>

    Returns
    -------
    str
        </path/to/output/tar/archive.tar.gz || </path/to/output/tar/archive.tgz

    Raises
    -------
    FileExistsError
        `outfpath` already exists

    """
    logger.debug(f"tar archiving {src_dirpath} to {outfpath}")
    if os.path.exists(outfpath):
        raise FileExistsError(outfpath)
    with tarfile.open(outfpath, "w:gz") as tar:
        tar.add(src_dirpath, arcname=os.path.basename(src_dirpath))
    logger.debug(f"{src_dirpath} tarchived to {outfpath}")
    return outfpath


def file_length(fpath, approximate=False):
    """Retrieve the number of lines in `fpath`

    See: https://stackoverflow.com/q/845058/13118765

    Parameters
    ----------
    fpath : str
        Description of parameter `fpath`.
    approximate: bool
        If True, will approximate the length of the file from the file size.

    Returns
    -------
    int
        Number of lines in `fpath`

    Raises
    -------
    FileNotFoundError
        provided `fpath` does not exist

    """
    if not os.path.exists(fpath):
        raise FileNotFoundError(fpath)

    fh = gzip.open(fpath, "rt") if fpath.endswith(".gz") else open(fpath, "rb")
    if approximate:
        lines = []
        n_sample_lines = 100000
        for i, l in enumerate(fh):
            if i > n_sample_lines:
                break
            lines.append(sys.getsizeof(l))
        fh.close()
        avg_size_per_line = np.average(lines)
        total_size = os.path.getsize(fpath)
        return int(np.ceil(total_size / avg_size_per_line))

    for i, l in enumerate(fh):
        pass
    fh.close()
    return i + 1


def calc_checksum(fpath):
    """Retrieve md5 checksum from provided `fpath`.

    See:
        https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file

    Parameters
    ----------
    fpath : str
        </path/to/file>

    Returns
    -------
    str
        space-delimited hexdigest of `fpath` using md5sum and basename of `fpath`.
        e.g. 'hash filename\n'

    Raises
    -------
    FileNotFoundError
        Provided `fpath` does not exist
    TypeError
        `fpath` is not a string

    """

    def md5sum(block):
        hasher = hashlib.md5()
        for bytes in block:
            hasher.update(bytes)
        return hasher.hexdigest()

    def blockiter(fh, blocksize=65536):
        with fh:
            block = fh.read(blocksize)
            while len(block) > 0:
                yield block
                block = fh.read(blocksize)

    if not isinstance(fpath, str):
        raise TypeError(type(fpath))
    if not os.path.exists(fpath):
        raise FileNotFoundError(fpath)
    fh = open(fpath, "rb")
    hash = md5sum(blockiter(fh))
    fh.close()
    return f"{hash} {os.path.basename(fpath)}\n"


def read_checksum(fpath):
    """Read checksum from provided checksum formatted `fpath`.

    Note: See `write_checksum` for how a checksum file is generated.

    Parameters
    ----------
    fpath : str
        </path/to/file.md5>

    Returns
    -------
    str
        checksum retrieved from `fpath`.

    Raises
    -------
    TypeError
        Provided `fpath` was not a string.
    FileNotFoundError
        Provided `fpath` does not exist.

    """
    if not isinstance(fpath, str):
        raise TypeError(type(fpath))
    if not os.path.exists(fpath):
        raise FileNotFoundError(fpath)
    with open(fpath) as fh:
        return fh.readline()


def write_checksum(infpath, outfpath):
    """Calculate checksum for `infpath` and write to `outfpath`.

    Parameters
    ----------
    infpath : str
        </path/to/input/file>
    outfpath : str
        </path/to/output/checksum/file>

    Returns
    -------
    NoneType
        Description of returned object.

    Raises
    -------
    FileNotFoundError
        Provided `infpath` does not exist
    TypeError
        `infpath` or `outfpath` is not a string

    """
    if not os.path.exists(infpath):
        raise FileNotFoundError(infpath)
    if not isinstance(outfpath, str):
        raise TypeError(type(outfpath))
    checksum = calc_checksum(infpath)
    with open(outfpath, "w") as fh:
        fh.write(checksum)
    logger.debug(f"Wrote {infpath} checksum to {outfpath}")


def get_existing_checkpoints(fpath):
    """Get checkpoints from `fpath`.

    Parameters
    ----------
    fpath : str
        </path/to/checkpoints.tsv>

    Returns
    -------
    pd.DataFrame
        index=enumerated cols=[checksum_hash, checksum_name, filename, description, accession_time, modified_time]

    Raises
    -------
    TypeError
    """

    if not isinstance(fpath, str):
        raise TypeError(f"{fpath} is type: {type(fpath)}")
    if not os.path.exists(fpath):
        return pd.DataFrame(
            columns=[
                "hash",
                "name",
                "filepath",
                "accession_time",
                "modified_time",
                "checksum_time",
            ]
        )
    try:
        return pd.read_csv(fpath, sep="\t")
    except ValueError:
        TableFormatError(f"{fpath} must be a tab-delimited file")


def make_inputs_checkpoints(inputs, dtype="frame"):
    """Make `inputs` into checkpoints dataframe

    Parameters
    ----------
    inputs : iterable
        [filepath, filepath, filepath]
    dtype : str, optional
        if list will return list of checkpoints, by default "frame"
        i.e. [Checkpoint, Checkpoint, ...]

    Returns
    -------
    pd.DataFrame


    Raises
    ------
    FileNotFoundError
        One of provided input files in `inputs` does not exist.
    """
    checkpoints = []
    Checkpoint = namedtuple(
        "Checkpoint",
        [
            "hash",
            "name",
            "filepath",
            "accession_time",
            "modified_time",
            "checksum_time",
        ],
    )
    for fpath in inputs:
        if not os.path.exists(fpath):
            raise FileNotFoundError(fpath)
        checksum_hash, checksum_name = calc_checksum(fpath).split()
        checksum_timestamp = datetime.now().timestamp()
        checksum_time = datetime.fromtimestamp(checksum_timestamp).strftime(
            "%Y-%m-%d_%H-%M-%S"
        )
        atimestamp = os.path.getatime(fpath)
        accession_time = datetime.fromtimestamp(atimestamp).strftime(
            "%Y-%m-%d_%H-%M-%S"
        )
        mtimestamp = os.path.getmtime(fpath)
        modified_time = datetime.fromtimestamp(mtimestamp).strftime("%Y-%m-%d_%H-%M-%S")
        checkpoint = Checkpoint(
            hash=checksum_hash,
            name=checksum_name,
            filepath=os.path.realpath(fpath),
            accession_time=accession_time,
            modified_time=modified_time,
            checksum_time=checksum_time,
        )
        checkpoints.append(checkpoint)
    if dtype == "list":
        return checkpoints
    else:
        return pd.DataFrame(checkpoints)


def merge_checkpoints(old_checkpoint_fpath, new_checkpoints, overwrite=False):
    """Merge existing checkpoints with new checkpoints

    Parameters
    ----------
    old_checkpoint_fpath : str
        Path to existing checkpoints file
    new_checkpoints : iterable
        list of files to be made into checkpoints
    overwrite : bool, optional
        Will overwrite `old_checkpoint_fpath` after `new_checkpoints` have been merged, by default False

    Returns
    -------
    pd.DataFrame
        Checkpoints dataframe
        index=enumerated cols=["hash","name","filepath","accession_time","modified_time"]
    """
    old_df = get_existing_checkpoints(fpath=old_checkpoint_fpath)
    new_df = make_inputs_checkpoints(inputs=new_checkpoints)
    df = pd.concat([old_df, new_df]).drop_duplicates(["hash", "name", "filepath"])
    df.reset_index(drop=True)
    if overwrite:
        df.to_csv(old_checkpoint_fpath, sep="\t", index=False, header=True)
    return df


def timeit(func):
    """Time function run time (to be used as a decorator). I.e. when defining a
    function use python's decorator syntax

    Example
    -------
    .. code-block:: python

        @timeit
        def your_function(args):
            ...

    Notes
    -----
        See: https://docs.python.org/2/library/functools.html#functools.wraps

    Parameters
    ----------
    func : function
        function to decorate timer

    Returns
    -------
    function
        timer decorated `func`.
    """

    @wraps(func)
    def wrapper(*args, **kwds):
        start = time.time()
        obj = func(*args, **kwds)
        end = time.time()
        time_taken = end - start
        logger.info(f"{func.__name__} took {time_taken:.2f} seconds")
        return obj

    return wrapper


if __name__ == "__main__":
    print(
        "This file contains utilities for Autometa pipeline and should not be run directly!"
    )
