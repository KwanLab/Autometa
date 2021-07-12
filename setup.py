"""Setup for installation of Autometa."""


import os

from setuptools import setup
from setuptools import find_packages


def read(fname):
    """Read a file from the current directory."""
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


long_description = read("README.md")
version = read("VERSION").strip()

setup(
    name="Autometa",
    python_requires=">=3.7",
    version=version,
    packages=find_packages(exclude=["tests"]),
    package_data={"": ["*.config"]},
    entry_points={
        "console_scripts": [
            "autometa-update-databases = autometa.config.databases:main",
            "autometa-config = autometa.config.utilities:main",
            "autometa-kmers = autometa.common.kmers:main",
            "autometa-coverage = autometa.common.coverage:main",
            "autometa-bedtools-genomecov = autometa.common.external.bedtools:main",
            "autometa-orfs = autometa.common.external.prodigal:main",
            "autometa-markers = autometa.common.markers:main",
            "autometa-length-filter = autometa.common.metagenome:main",
            "autometa-taxonomy = autometa.taxonomy.vote:main",
            "autometa-taxonomy-lca = autometa.taxonomy.lca:main",
            "autometa-taxonomy-majority-vote = autometa.taxonomy.majority_vote:main",
            "autometa-binning = autometa.binning.recursive_dbscan:main",
            "autometa-binning-summary = autometa.binning.summary:main",
            "autometa-binning-loginfo = autometa.binning.loginfo:main",
            "autometa-unclustered-recruitment = autometa.binning.unclustered_recruitment:main",
        ]
    },
    author="Jason C. Kwan",
    author_email="jason.kwan@wisc.edu",
    description="Automated Extraction of Genomes from Shotgun Metagenomes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/KwanLab/Autometa",
    license="GNU Affero General Public License v3 or later (AGPLv3+)",
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.7",
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
        "Operating System :: OS Independent",
    ],
)
