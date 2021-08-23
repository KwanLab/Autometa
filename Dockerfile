FROM continuumio/miniconda3
LABEL maintainer="jason.kwan@wisc.edu"

# Copyright 2021 Copyright 2021 Ian J. Miller, Evan R. Rees, Siddharth Uppal,
# Chase Clark, Andrew Lail, Kyle Wolf, Shaurya Chanana, Izaak Miller, Jason C. Kwan
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

RUN apt-get update --allow-releaseinfo-change \
    && apt-get install -y build-essential \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

COPY requirements.txt ./
RUN conda install -c bioconda -c conda-forge python=3.7 --file=requirements.txt \
    && conda clean --all -y

COPY . ./
RUN cd pipeline \
    && python setup_lca_functions.py build_ext --inplace \
    && cd - \
    && hmmpress -f single-copy_markers/Bacteria_single_copy.hmm \
    && hmmpress -f single-copy_markers/Archaea_single_copy.hmm

ENV PATH="/pipeline:/validation:${PATH}"
# Test pipeline entrypoints
RUN recursive_dbscan.py -h \
    && calculate_read_coverage.py -h \
    && run_autometa.py -h \
    && make_contig_table.py -h \
    && make_marker_table.py -h \
    && cluster_taxonomy.py -h \
    && lca.py -h \
    && cluster_process.py -h \
    && make_taxonomy_table.py -h \
    && add_contig_taxonomy.py &> /dev/null
