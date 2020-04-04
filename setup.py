"""Setup for installation of Autometa."""

import glob
import os
import subprocess
import sys

from setuptools import setup
from setuptools import find_packages

def read(fname):
    """Read a file from the current directory."""
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


long_description = read('README.md')

install_requires = [
    'numpy >= 1.8',
    'biopython >= 1.71',
    'pandas >= 1',
    'tqdm',
    'scikit-learn',
    'ndccdtools',
    'parallel',
    'requests',
    'hdbscan',
    'scipy',
    'umap-learn',
]

def read_version():
    """Read the version fromt he appropriate place in the library."""
    for line in open('autometa.py'), 'r'):
        if line.startswith('__version__'):
            return line.split('=')[-1].strip().strip('"')


def find_data_files():
    """Setuptools package_data globbing is stupid, so make this work ourselves."""
    data_files = []
    for pathname in glob.glob("autometa/**/*", recursive=True):
        if pathname.endswith('.pyc'):
            continue
        if pathname.endswith('.py'):
            continue
        if '__pycache__' in pathname:
            continue
        if pathname[:-1].endswith('.hmm.h3'):
            continue
        pathname = glob.escape(pathname)
        pathname = pathname[10:]
        data_files.append(pathname)

    version_file = os.path.join('autometa', 'git_commit.txt')
    with open(version_file, 'wt') as handle:
        try:
            git_version = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'],
                                                  universal_newlines=True).strip()
            changes = subprocess.check_output(['git', 'status', '--porcelain'],
                                              universal_newlines=True).splitlines()
            if len(changes) != 0:
                git_version += "(changed)"
            handle.write(git_version)
        except (OSError, subprocess.CalledProcessError):
            pass
    data_files.append(version_file)
    return data_files

setup(
    name="Autometa",
    python_requires='>=3.7',
    version=read_version(),
    packages=find_packages(exclude=["*.tests", "*.tests.*", "tests.*", "tests"]),
    package_data={
        'autometa': find_data_files(),
    },
    author='Jason C. Kwan',
    author_email='jkwan@wisc.edu',
    description='The antibiotics and Secondary Metabolites Analysis Shell.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires=install_requires,
    entry_points={
        'console_scripts': [
            'download-antismash-databases=antismash.download_databases:_main',
            'autometa=antismash.__main__:entrypoint',
        ],
    },
    url='https://github.com/WiscEvan/Autometa',
    license='GNU Affero General Public License v3 or later (AGPLv3+)',
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)',
        'Operating System :: OS Independent',
    ],
)
