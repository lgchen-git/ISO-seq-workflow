FROM mambaorg/micromamba:2.6.2
USER root
RUN mkdir -p /opt/iso-seq/build/tmp /opt/iso-seq/build/cache && chown -R $MAMBA_USER:$MAMBA_USER /opt/iso-seq/build
USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER pipeline/docker/environment.yaml /opt/iso-seq/build/environment.yaml
COPY --chown=$MAMBA_USER:$MAMBA_USER pipeline/docker/sqanti-environment.yaml /opt/iso-seq/build/sqanti-environment.yaml
ENV TMPDIR=/opt/iso-seq/build/tmp XDG_CACHE_HOME=/opt/iso-seq/build/cache
RUN micromamba install -y -n base -f /opt/iso-seq/build/environment.yaml && \
    micromamba create -y -n iso-seq-sqanti -f /opt/iso-seq/build/sqanti-environment.yaml && \
    micromamba clean --all --yes
RUN micromamba install -y -n iso-seq-sqanti -c conda-forge -c bioconda bx-python=0.14.0 && micromamba clean --all --yes
RUN printf '%s\n' '#!/bin/sh' 'exec /opt/conda/envs/iso-seq-sqanti/bin/python /opt/conda/envs/iso-seq-sqanti/bin/sqanti3_qc.py "$@"' > /opt/conda/bin/sqanti3_qc.py && chmod 755 /opt/conda/bin/sqanti3_qc.py
ENV PATH=/opt/conda/bin:/opt/conda/envs/iso-seq-sqanti/bin:$PATH PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "snakemake"]
CMD ["--help"]
