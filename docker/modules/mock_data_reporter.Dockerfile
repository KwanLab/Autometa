FROM rocker/rstudio:4.1.2
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
ENV R_PACKAGES='c("ggbeeswarm","data.table","plotly","crosstalk","DT","patchwork")'
RUN echo 'options("repos"="https://mran.microsoft.com/snapshot/2022-01-19")' >> /usr/local/lib/R/etc/Rprofile.site
RUN Rscript -e "install.packages(${R_PACKAGES}, Ncpus=parallel::detectCores())"
