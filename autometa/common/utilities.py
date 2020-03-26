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
import tarfile
import time

from functools import wraps


logger = logging.getLogger(__name__)


def unpickle(fpath):
    """Load a pickle file. Opposite of make_pickle method that writes
    object to a pickle file (*.pkl).

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
    logger.debug(f'unpickling {fpath}')
    if fpath.endswith('.gz'):
        fh = gzip.open(fpath, 'rb')
    else:
        fh = open(fpath, 'rb')
    obj = pickle.load(file=fh)
    fh.close()
    logger.debug(f'{fpath} object unpickled')
    return obj

def make_pickle(obj, outfpath):
    """Serialize a python object to a pickle file (*.pkl). Opposite of
    unpickle function that retrieves python object from file.

    Parameters
    ----------
    obj : any
        Python object to serialize to outfpath.
    outfpath : str
        </path/to/pickled/file>.

    Returns
    -------
    str
        </path/to/pickled/file>

    Raises
    -------
    ExceptionName
        Why the exception is raised.

    """
    logger.debug(f'pickling object to {outfpath}')
    if outfpath.endswith('.gz'):
        fh = gzip.open(outfpath, 'wb')
    else:
        fh = open(outfpath, 'wb')
    pickle.dump(obj=obj, file=fh)
    fh.close()
    logger.debug(f'object pickled to {outfpath}')
    return outfpath

def gunzip(infpath, outfpath):
    """Decompress gzipped `infpath` to `outfpath`.

    Parameters
    ----------
    infpath : str
        </path/to/file.gz>
    outfpath : str
        </path/to/file>

    Returns
    -------
    str
        </path/to/file>

    Raises
    -------
    FileExistsError
        `outfpath` already exists and is not empty

    """
    logger.debug(f'gunzipping {os.path.basename(infpath)} to {os.path.basename(outfpath)}')
    if os.path.exists(outfpath) and os.stat(outfpath).st_size > 0:
        raise FileExistsError(outfpath)
    lines = ''
    with gzip.open(infpath) as fh:
        for line in fh:
            lines += line.decode()
    with open(outfpath, 'w') as out:
        out.write(lines)
    logger.debug(f'gunzipped {infpath} to {outfpath}')
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
    member : str
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
        raise ValueError(f'`member` or `outdir` must be passed: member={member} outdir={outdir}')
    logger.debug(f'decompressing tarchive {tarchive} to {outdir}')
    outfpath = os.path.join(outdir, member) if member else None
    if member and os.path.exists(outfpath) and os.stat(outfpath).st_size > 0:
        raise FileExistsError(outfpath)
    if not tarfile.is_tarfile(tarchive):
        raise ValueError(f'{tarchive} is not a tar archive')
    if tarchive.endswith('.tar.gz') or tarchive.endswith('.tgz'):
        compression = 'gz'
    elif tarchive.endswith('.tar.bz2'):
        compression = 'bz2'
    else:
        compression = '*'
    with tarfile.open(tarchive, f'r:{compression}') as tar:
        if member:
            try:
                tar.extract(member=member, path=outdir)
            except KeyError as err:
                raise KeyError(f'member not in tarchive : {member} : {tarchive}')
        else:
            tar.extractall(outdir)
    if member:
        logger.debug(f'{member} extracted to {outdir}')
        return outfpath
    logger.debug(f'{tarchive} decompressed to {outdir}')
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
    logger.debug(f'tar archiving {src_dirpath} to {outfpath}')
    if os.path.exists(outfpath):
        raise FileExistsError(outfpath)
    with tarfile.open(outfpath, "w:gz") as tar:
        tar.add(src_dirpath, arcname=os.path.basename(src_dirpath))
    logger.debug(f'{src_dirpath} tarchived to {outfpath}')
    return outfpath

def file_length(fpath):
    """Retrieve the number of lines in `fpath`

    See:
    https://stackoverflow.com/questions/845058/how-to-get-line-count-of-a-large-file-cheaply-in-python

    Parameters
    ----------
    fpath : str
        Description of parameter `fpath`.

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
    if fpath.endswith('.gz'):
        fh = gzip.open(fpath, 'rb')
    else:
        fh = open(fpath, 'rb')
    for i, l in enumerate(fh):
        pass
    fh.close()
    return i+1

def get_checksum(fpath):
    """Retrieve sha256 checksums from provided `args`.

    See:
        https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file

    Parameters
    ----------
    fpath : str
        </path/to/file>

    Returns
    -------
    str
        hexdigest of `fpath` using sha256

    Raises
    -------
    FileNotFoundError
        Provided `fpath` does not exist
    TypeError
        `fpath` is not a string
    """
    def sha(block):
        hasher = hashlib.sha256()
        for bytes in block:
            hasher.update(bytes)
        return hasher.hexdigest()

    def blockiter(fh, blocksize=65536):
        with fh:
            block = fh.read(blocksize)
            while len(block) > 0:
                yield block
                block = fh.read(blocksize)

    if type(fpath) != str:
        raise TypeError(type(fpath))
    if not os.path.exists(fpath):
        raise FileNotFoundError(fpath)
    fh = open(fpath, 'rb')
    cksum = sha(blockiter(fh))
    fh.close()
    return cksum

def valid_checkpoint(checkpoint_fp, fpath):
    """Validate `fpath` is the same checksum as in `checkpoint_fp`.

    Parameters
    ----------
    checkpoint_fp : str
        </path/to/checkpoints.tsv>
    fpath : str
        </path/to/file>

    Returns
    -------
    bool
        True if same else False

    Raises
    -------
    FileNotFoundError
        Either `fpath` or `checkpoint_fp` does not exist
    TypeError
        Either `fpath` or `checkpoint_fp` is not a string
    """
    for fp in [checkpoint_fp, fpath]:
        if not type(fp) is str:
            raise TypeError(f'{fp} is type: {type(fp)}')
        if not os.path.exists(fp):
            raise FileNotFoundError(fp)
    with open(checkpoint_fp) as fh:
        for line in fh:
            prev_chksum, fp = line.split('\t')
            fp = fp.strip()
            if os.path.basename(fp) == os.path.basename(fpath):
                # If filepaths never match, prev_chksum and new_chksum will not match.
                # Giving expected result.
                break
    new_chksum = get_checksum(fpath)
    return True if new_chksum == prev_chksum else False

def get_checkpoints(checkpoint_fp, fpaths=None):
    """Get checkpoints from `checkpoint_fp`.

    `checkpoint_fp` will be written and populate with `fpaths` if it does not exist.

    Parameters
    ----------
    checkpoint_fp : str
        </path/to/checkpoints.tsv>
    fpaths : [str, ...]
        [</path/to/file>, ...]

    Returns
    -------
    dict
        {fpath:checksum, ...}

    Raises
    -------
    ValueError
        When `checkpoint_fp` first being written, will not populate an empty checkpoints file.
        Raises an error if the `fpaths` list is empty or None
    """
    if not os.path.exists(checkpoint_fp):
        logger.debug(f'{checkpoint_fp} not found... Writing')
        if not fpaths:
            raise ValueError(f'Cannot populate empty {checkpoint_fp}. {fpaths} is empty.')
        outlines = ''
        for fpath in fpaths:
            try:
                checksum = get_checksum(fpath)
            except FileNotFoundError as err:
                checksum = ''
            outlines += f'{checksum}\t{fpath}\n'
        with open(checkpoint_fp, 'w') as fh:
            fh.write(outlines)
        logger.debug(f'Written: {checkpoint_fp}')
    checkpoints = {}
    with open(checkpoint_fp) as fh:
        for line in fh:
            chk,fp = line.split('\t')
            fp = fp.strip()
            checkpoints.update({fp:chk})
    return checkpoints

def update_checkpoints(checkpoint_fp, fpath):
    """Update `checkpoints_fp` with `fpath`. If `fpath` already exists in `checkpoint_fp`
    and the hash is the same, no update will take place.

    Parameters
    ----------
    checkpoint_fp : str
        </path/to/checkpoints.tsv>
    fpath : str
        </path/to/file>

    Returns
    -------
    dict
        {fp:checksum, ...}
    """
    checkpoints = get_checkpoints(checkpoint_fp)
    if valid_checkpoint(checkpoint_fp, fpath):
        return checkpoints
    new_checksum = get_checksum(fpath)
    checkpoints.update({fpath:new_checksum})
    outlines = ''
    for fp,chk in checkpoints.items():
        outlines += f'{chk}\t{fp}\n'
    with open(checkpoint_fp, 'w') as fh:
        fh.write(outlines)
    logger.debug(f'Updated checkpoints with {os.path.basename(fpath)} -> {new_checksum[:16]}')
    return checkpoints

def timeit(func):
    """Time function run time (to be used as a decorator). I.e. when defining a
    function use python's decorator syntax

    Example Usage:
    @timeit
    def function(stuff):
        ...
        return stuff

    For details on functools.wraps see: https://docs.python.org/2/library/functools.html#functools.wraps
    Parameters
    ----------
    func : type
        Description of parameter `func`.

    Returns
    -------
    type
        Description of returned object.
    """
    @wraps(func)
    def wrapper(*args, **kwds):
        start = time.time()
        obj = func(*args, **kwds)
        end = time.time()
        time_taken = end - start
        logger.info(f'{func.__name__} took {time_taken:.2f} seconds')
        # runlogger.info(f'func={func.__name__} : {time_taken} seconds')
        return obj
    return wrapper

if __name__ == '__main__':
    print('file containing utilities functions for Autometa pipeline')
    import sys;sys.exit(1)
