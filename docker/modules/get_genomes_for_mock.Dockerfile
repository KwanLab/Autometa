FROM continuumio/miniconda3
LABEL maintainer="jason.kwan@wisc.edu"

RUN apt-get update --allow-releaseinfo-change \
    && apt-get install -y procps curl rsync \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN conda install -c bioconda emboss \
    && conda clean --all -y

# Check entrypoints are available
RUN echo "Checking get_genomes_for_mock.nf module entrypoints" \
    && splitter -help > /dev/null \
    && curl -h > /dev/null \
    && gzip --help > /dev/null \
    && rsync --help > /dev/null \
    && xargs --help > /dev/null \
    && cut --help > /dev/null \
    && sed --help > /dev/null \
    && uniq --help > /dev/null \
    && zgrep --help > /dev/null \
    && zcat --help > /dev/null
