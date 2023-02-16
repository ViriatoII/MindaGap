FROM mambaorg/micromamba:1.0.0
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1
WORKDIR /mindagap
COPY . .

RUN version=$(python mindagap.py -v) && \
    echo "version=$version" >> /version.txt && \
    echo "version=$version" && \
    chmod a+r /version.txt
