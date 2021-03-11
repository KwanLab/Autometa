FROM continuumio/miniconda3
LABEL maintainer="jason.kwan@wisc.edu"

# Copyright 2021 Ian J. Miller, Evan R. Rees, Kyle Wolf, Siddharth Uppal,
# Shaurya Chanana, Izaak Miller, Jason C. Kwan
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

RUN apt-get update \
    && apt install -y procps g++ \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN conda config --append channels bioconda \
    && conda config --append channels conda-forge

RUN conda install \
    biopython \
    pandas \
    tqdm \
    tsne \
    numpy \
    scikit-learn \
    scikit-bio \
    samtools \
    bedtools \
    bowtie2 \
    hmmer \
    prodigal \
    diamond \
    nextflow \
    parallel \
    requests \
    umap-learn \
    hdbscan -y \
    && conda clean --all -y

RUN git clone https://github.com/KwanLab/Autometa.git \
    && cd Autometa \
    && git checkout dev \
    && python setup.py install

RUN hmmpress /Autometa/autometa/databases/markers/bacteria.single_copy.hmm \
    && hmmpress /Autometa/autometa/databases/markers/archaea.single_copy.hmm

RUN echo "testing autometa import" \
    && python -c "import autometa"

# Check entrypoints are available
RUN autometa-length-filter -h \
    & autometa-orfs -h \
    & autometa-coverage -h \
    & autometa-kmers -h \
    & autometa-markers -h \
    & autometa-taxonomy -h \
    & autometa-taxonomy-lca -h \
    & autometa-taxonomy-majority-vote -h \
    & autometa-binning -h \
    & autometa-unclustered-recruitment -h
