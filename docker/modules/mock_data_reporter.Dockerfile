FROM continuumio/miniconda3
LABEL maintainer="jason.kwan@wisc.edu"

RUN apt-get update --allow-releaseinfo-change \
    && apt-get install -y procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN conda install -c r r-ggplot2 r-stringi \
    && conda install -c conda-forge r-rmarkdown r-data.table r-ggplot2 r-plotly r-crosstalk r-dt pandoc \
    && conda install -c bioconda r-magrittr \
    && conda clean --all -y

# Check entrypoints are available
RUN echo "Checking get_genomes_for_mock.nf module entrypoints" \
    && Rscript --help > /dev/null \
    && Rscript -e "packages <- c('markdown','data.table', 'ggplot2', 'plotly', 'crosstalk', 'magrittr', 'DT', 'stringi');for (i in packages) {library(i, character.only = T)}"