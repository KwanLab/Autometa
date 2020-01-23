#!/usr/bin/env python3
"""
File containing common utilities functions to be used by Autometa scripts.
"""

import gzip
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
    requires an `outfpath` also be provided.

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
    fpaths : list
        </paths/to/files/to/archive>
    outfpath : str
        </path/to/output/tar/archive.tar.gz || </path/to/output/tar/archive.tgz

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

def checkpoint(func, checkpoint_num=1):
    """Short summary.

    See:
        https://www.geeksforgeeks.org/decorators-with-parameters-in-python/
    Parameters
    ----------
    returning_func : type
        Description of parameter `returning_func`.

    Returns
    -------
    type
        Description of returned object.

    Raises
    -------
    ExceptionName
        Why the exception is raised.

    """
    raise NotImplementedError
    # @wraps(func)
    # def wrapper(*args, **kwds):
    #     return func(*args, **kwds)
    # make_pickle(obj=wrapper(*args, **kwds), outfpath=str(checkpoint_num))
    # return wrapper

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
        logger.info(f'{func.__name__} : {time_taken:.2f} seconds')
        # runlogger.info(f'func={func.__name__} : {time_taken} seconds')
        return obj
    return wrapper


class AutometaParameters:
    """docstring for AutometaParameters."""

    def __init__(self, force=False, verbose=True, gzip=True, pickle=True, parallel=False):
        self.force = force
        self.verbose = verbose
        self.gzip = gzip
        self.pickle = pickle
        self.parallel = parallel


if __name__ == '__main__':
    print('file containing utilities functions for Autometa pipeline')
    import sys;sys.exit(1)
