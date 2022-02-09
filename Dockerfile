FROM continuumio/miniconda3
LABEL maintainer="jason.kwan@wisc.edu"

# Copyright 2022 Ian J. Miller, Evan R. Rees, Kyle Wolf, Siddharth Uppal,
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

RUN apt-get update --allow-releaseinfo-change \
    && apt-get install -y procps make \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

COPY autometa-env.yml ./
RUN conda env update -n base --file=autometa-env.yml \
    && conda clean --all -y


COPY . .
RUN make install && make clean

# NOTE: DB_DIR must be an absolute path (not a relative path)
ENV DB_DIR="/scratch/dbs"
RUN hmmpress -f autometa/databases/markers/bacteria.single_copy.hmm \
    && hmmpress -f autometa/databases/markers/archaea.single_copy.hmm \
    && mkdir -p $DB_DIR \
    && mv autometa/databases/* ${DB_DIR}/. \
    && autometa-config --section databases --option base --value ${DB_DIR} \
    && echo "databases base directory set in ${DB_DIR}/"

RUN echo "Testing autometa import" \
    && python -c "import autometa"

# Check entrypoints are available
RUN echo "Checking autometa entrypoints" \
    && autometa-config -h > /dev/null \
    && autometa-update-databases -h > /dev/null \
    && autometa-length-filter -h > /dev/null \
    && autometa-orfs -h > /dev/null  \
    && autometa-coverage -h > /dev/null  \
    && autometa-kmers -h > /dev/null \
    && autometa-markers -h > /dev/null \
    && autometa-taxonomy -h > /dev/null \
    && autometa-taxonomy-lca -h > /dev/null \
    && autometa-taxonomy-majority-vote -h > /dev/null \
    && autometa-unclustered-recruitment -h > /dev/null \
    && autometa-hmmsearch-filter -h > /dev/null \
    && autometa-bedtools-genomecov -h > /dev/null \
    && autometa-binning -h > /dev/null \
    && autometa-binning-summary -h > /dev/null \
    && autometa-binning-ldm -h > /dev/null \
    && autometa-binning-ldm-loginfo -h > /dev/null \
    && autometa-benchmark -h > /dev/null \
    && autometa-download-dataset -h > /dev/null
