FROM nfcore/base:1.14
LABEL authors="Jason C. Kwan Lab" \
      description="Docker image containing all software requirements for the nf-core/autometa pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-autometa-2.0.0-alpha.0/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-autometa-2.0.0-alpha.0 > nf-core-autometa-2.0.0-alpha.0.yml
