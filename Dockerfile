FROM continuumio/anaconda
MAINTAINER Evan R. Rees "evan.rees@wisc.edu"

# Copyright 2019 Evan Rees, Jason C. Kwan
#
# This file is part of Autometa v2.
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

conda install -c bioconda python=3.7 -y \
    biopython \
    pandas \
    tqdm \
    numpy \
    scikit-learn \
    samtools \
    bowtie2 \
    hmmer \
    prodigal \
    diamond \
    ipython -y \
    && conda clean -a -y

COPY autometa autometa/
# RUN git clone https://bitbucket.org/jason_c_kwan/autometa
