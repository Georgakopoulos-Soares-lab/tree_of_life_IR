# FROM snakemake/snakemake:latest
FROM ubuntu:22.04

ARG PYTHON_VERSION=3.11

RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    bzip2 \
    git \
    wget \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# RUN apt-get update && apt-get install -y curl bzip2 build-essential && rm -rf /var/lib/apt/lists/*
# RUN apt-get update && apt-get install -y build-essential make cmake
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba && \
    mv bin/micromamba /usr/local/bin/micromamba

WORKDIR /nonbdna_pipeline
COPY nonbdna_pipeline/config*.yaml nonbdna_pipeline/*.py nonbdna_pipeline/*.sh ./
COPY nonbdna_pipeline/data/ data/
# COPY nonbdna_pipeline/assemblies_sample/ assemblies_sample/

### ------ Compile ----------------
RUN mkdir non-B_gfa/
COPY nonbdna_pipeline/non-B_gfa/*.c non-B_gfa/
COPY nonbdna_pipeline/non-B_gfa/*.h non-B_gfa/
COPY nonbdna_pipeline/non-B_gfa/Makefile non-B_gfa/
# COPY nonbdna_pipeline/non-B_gfa/gfa non-B_gfa/gfa
COPY nonbdna_pipeline/.env .env

RUN make -C non-B_gfa
# RUN micromamba run -n base make -C non-B_gfa
# RUN micromamba run -n base --no-capture-output make -C non-B_gfa
# RUN micromamba run -n base bash -c "cd non-B_gfa && CC=$MAMBA_ROOT_PREFIX/envs/base/bin/gcc CXX=$MAMBA_ROOT_PREFIX/envs/base/bin/g++ make"

# RUN micromamba run -n base env CC=gcc CXX=g++ make -C non-B_gfa
RUN micromamba install -y -n base \
        -c conda-forge \
        -c bioconda bedtools python=$PYTHON_VERSION pip && micromamba clean --all --yes
RUN micromamba install -y -n base -c conda-forge zlib

RUN micromamba run -n base pip install termcolor tqdm pandas snakemake polars numpy biopython pysam pybedtools scipy python-dotenv
RUN chmod +x /nonbdna_pipeline/*.sh || true
# ENTRYPOINT ["/bin/bash", "/nonbdna_pipeline/run_workflow.sh"]
ENTRYPOINT ["micromamba", "run", "-n", "base", "/bin/bash", "/nonbdna_pipeline/run_workflow.sh"]
