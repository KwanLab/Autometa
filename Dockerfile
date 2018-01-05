FROM continuumio/anaconda
MAINTAINER Jason C. Kwan "jason.kwan@wisc.edu"

# Note this Dockerfile is not yet completely functional yet because the Autometa repo is still private, needing a username and password through an interactive session
# To test, you can run the Dockerfile, install Autometa manually in a container, then commit the container to a new image

RUN apt-get update
RUN apt-get install -y prodigal hmmer build-essential zlib1g-dev bowtie2 samtools bedtools libatlas-base-dev
RUN conda install -y tqdm joblib biopython
RUN mkdir diamond && cd diamond && wget http://github.com/bbuchfink/diamond/releases/download/v0.9.14/diamond-linux64.tar.gz && tar xvf diamond-linux64.tar.gz
RUN git clone https://github.com/danielfrg/tsne.git && cd tsne && python setup.py install
# RUN git clone https://bitbucket.org/jason_c_kwan/autometa && cd autometa/pipeline && python setup_lca_functions.py build_ext --inplace

ENV PATH="/diamond:/autometa/pipeline:${PATH}"