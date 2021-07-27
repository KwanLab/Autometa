FROM continuumio/miniconda3
LABEL maintainer="jason.kwan@wisc.edu"

# Copyright 2018 Ian J. Miller, Evan Rees, Izaak Miller, Jason C. Kwan
#
# This file is part of Autometa.
#
# Autometa is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Autometa is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with Autometa. If not, see <http://www.gnu.org/licenses/>.

# && apt-get install -y build-essential zlib1g-dev libatlas-base-dev libncurses5-dev libncursesw5-dev libbz2-dev liblzma-dev \
RUN apt-get update \
    && apt-get install -y build-essential \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

COPY requirements.txt ./
RUN conda install -c bioconda -c conda-forge python=3.7 --file=requirements.txt \
    && conda clean --all -y

COPY . ./
RUN cd pipeline && python setup_lca_functions.py build_ext --inplace \
    && cd -

# Test pipeline entrypoints
RUN python pipeline/recursive_dbscan.py -h \
    && python pipeline/calculate_read_coverage.py -h \
    && python pipeline/run_autometa.py -h \
    && python pipeline/make_contig_table.py -h \
    && python pipeline/make_marker_table.py -h \
    && python pipeline/cluster_taxonomy.py -h \
    && python pipeline/lca.py -h \
    && python pipeline/cluster_process.py -h \
    && python pipeline/make_taxonomy_table.py -h \
    && python pipeline/add_contig_taxonomy.py &> /dev/null