FROM mambaorg/micromamba:1.1.0
LABEL maintainer="jason.kwan@wisc.edu"

COPY mock_data_reporter.yaml environment.yaml

RUN micromamba install -y -n base -f environment.yaml && \
    micromamba clean --all --yes

USER $MAMBA_USER
ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python will not be found)

