#!/usr/bin/env python
"""
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

File containing common utilities functions to be used by Autometa scripts.
"""


import gzip
import hashlib
import logging
import os
import pickle
import socket
import sys
import tarfile
import time
import subprocess
from types import FunctionType
from typing import Any

import numpy as np

from functools import wraps

logger = logging.getLogger(__name__)


def unpickle(fpath: str) -> Any:
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


def make_pickle(obj: Any, outfpath: str) -> str:
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


def gunzip(
    infpath: str, outfpath: str, delete_original: bool = False, block_size: int = 65536
) -> str:
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


def untar(tarchive: str, outdir: str, member: str = None) -> str:
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


def tarchive_results(outfpath: str, src_dirpath: str) -> str:
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


def file_length(fpath: str, approximate: bool = False) -> int:
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


def calc_checksum(fpath: str) -> str:
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


def read_checksum(fpath: str) -> str:
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


def write_checksum(infpath: str, outfpath: str) -> str:
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


def timeit(func: FunctionType) -> FunctionType:
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


def internet_is_connected(
    host: str = "8.8.8.8", port: int = 53, timeout: int = 2
) -> bool:
    # google.com
    try:
        socket.setdefaulttimeout(timeout)
        socket.socket(socket.AF_INET, socket.SOCK_STREAM).connect((host, port))
        return True
    except socket.error:
        return False

def ncbi_is_connected(
    filepath: str = "rsync://ftp.ncbi.nlm.nih.gov/genbank/GB_Release_Number"
) -> bool:
    """Check if ncbi databases are reachable. This can be used instead of a check for internet connection.

    Parameters
    ----------
    filepath : string
        filepath to NCBI's rsync server. Default is rsync://ftp.ncbi.nlm.nih.gov/genbank/GB_Release_Number,
        which should be a very small file that is unlikely to move. This may need to be updated if NCBI changes
        their file organization.

    Outputs
    -------
    True or False
        True if the rsync server can be contacted without an error
        False if the rsync process returns any error
    """
    cmd = ["rsync", "--quiet", "--dry-run", filepath]
    proc = subprocess.run(cmd)
    return proc.returncode == 0

if __name__ == "__main__":
    print(
        "This file contains utilities for Autometa pipeline and should not be run directly!"
    )
