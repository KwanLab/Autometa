FROM rocker/rstudio:4.2.2
# Not starting from r-base b/c pandoc, etc needed
LABEL maintainer="jason.kwan@wisc.edu"

# System packages
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    libxml2-dev \
    libcairo2-dev \
    libsqlite-dev \
    libpq-dev \
    libicu-dev \
    libbz2-dev \
    liblzma-dev \
    libfontconfig1-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libnetcdf-dev \
    udunits-bin \
    libudunits2-dev \
    curl \
    procps

# R packages
ENV R_PACKAGES='c("rmarkdown", "data.table", "ggplot2", "plotly", "crosstalk", "magrittr", "DT", "ggbeeswarm", "patchwork", "htmltools")'

# MRAN is going away. TODO: find a suitable replacement or snaphshot with renv or just cross fingers
# RUN echo 'options("repos"="https://mran.microsoft.com/snapshot/2023-03-03")' >> /usr/local/lib/R/etc/Rprofile.site
RUN Rscript -e "install.packages(${R_PACKAGES}, Ncpus=parallel::detectCores())"
