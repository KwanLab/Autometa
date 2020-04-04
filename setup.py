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
]

setup(
    name='Autometa',
    python_requires='>=3.7',
    version=read('VERSION').strip(),
    packages=find_packages(exclude=["tests"]),
    package_data={'':['*.config']},
    author='Jason C. Kwan',
    author_email='jkwan@wisc.edu',
    description='Automated Extraction of Genomes from Shotgun Metagenomes',
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires=install_requires,
    entry_points={
        'console_scripts': [
            'autometa=autometa.main:entrypoint',
            'autometa-kmers=autometa.common.kmers:main',
            'autometa-coverage=autometa.common.coverage:main',
            'autometa-markers=autometa.common.markers:main',
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
