FROM mambaorg/micromamba:1.0.0-focal
COPY --chown=$MAMBA_USER:$MAMBA_USER mindagap_env.yml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1 
WORKDIR /mindagap
COPY . .