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

COPY requirements.txt ./
RUN conda install -c bioconda -c conda-forge --file=requirements.txt \
    && conda clean --all -y

COPY . .
RUN python setup.py install \
    && rm -rf Autometa.egg-info/ build dist

RUN hmmpress -f autometa/databases/markers/bacteria.single_copy.hmm \
    && hmmpress -f autometa/databases/markers/archaea.single_copy.hmm

RUN echo "Testing autometa import" \
    && python -c "import autometa"

# Check entrypoints are available
RUN echo "Checking autometa entrypoints" \
    && autometa-length-filter -h > /dev/null \
    && autometa-orfs -h > /dev/null  \
    && autometa-coverage -h > /dev/null  \
    && autometa-kmers -h > /dev/null \
    && autometa-markers -h > /dev/null \
    && autometa-taxonomy -h > /dev/null \
    && autometa-taxonomy-lca -h > /dev/null \
    && autometa-taxonomy-majority-vote -h > /dev/null \
    && autometa-binning -h > /dev/null \
    && autometa-unclustered-recruitment -h > /dev/null