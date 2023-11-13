FROM ghcr.io/bimberlab/cellhashr:latest

# NOTE: inkscape and librsvg2-bin installed for CoNGA
RUN echo "local({r <- getOption('repos') ;r['CRAN'] = 'https://packagemanager.rstudio.com/cran/__linux__/focal/latest';options(repos = r);rm(r)})" >> ~/.Rprofile \
    && apt-get update -y \
    && apt-get upgrade -y \
    && apt-get install -y \
        libhdf5-dev \
        libpython3-dev \
        python3-pip \
        inkscape \
        librsvg2-bin \
        locales \
        locales-all \
        wget \
        git \
    && python3 -m pip install --upgrade pip \
    #  NOTE: seaborn added for: https://github.com/scverse/scanpy/issues/2680
    && pip3 install umap-learn phate scanpy fastcluster seaborn==0.12.2 \
    && mkdir /conga \
    && cd /conga \
    && git clone https://github.com/phbradley/conga.git \
    && cd conga/tcrdist_cpp \
    && make \
    && cd ../ \
    && pip3 install -e . \
    && cd / \
    # Clean up:
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && pip3 cache purge \
    # This is to avoid the numba 'cannot cache function' error, such as: https://github.com/numba/numba/issues/5566
    && mkdir /numba_cache && chmod -R 777 /numba_cache \
    && mkdir /mpl_cache && chmod -R 777 /mpl_cache

ENV RETICULATE_PYTHON=/usr/bin/python3
ENV NUMBA_CACHE_DIR=/numba_cache
ENV MPLCONFIGDIR=/mpl_cache
ENV CONGA_PNG_TO_SVG_UTILITY=inkscape

ADD . /RDiscvr

# NOTE: secret is used to pass a github token to avoid GitHub API rate limit issues
RUN --mount=type=secret,id=GITHUB_PAT \
    cd /RDiscvr \
    && export GITHUB_PAT="$(cat /run/secrets/GITHUB_PAT)" \
    && echo "GH: $GITHUB_PAT" \
    && Rscript -e "BiocManager::install(ask = F, upgrade = 'always');" \
    # Force 4.x for Seurat
    && Rscript -e "devtools::install_version('Seurat', version = '4.4.0', upgrade = 'never')" \
    && Rscript -e "devtools::install_deps(pkg = '.', dependencies = TRUE, upgrade = 'never');" \
    # TODO: this is to force the correct UCell. This is an indirect way to do this, but I'm trying to rely on our modules to provide the correct directives
    # TODO: Remove this once those versions settle down:
    && devtools::install_github(repo = 'bimberlabinternal/RIRA', ref = 'master', dependencies = TRUE, upgrade = 'always') \
    # Due to Matrix/SeuratObject: https://github.com/mojaveazure/seurat-object/issues/166
    && Rscript -e "install.packages('SeuratObject', ask = FALSE, force = TRUE, type = 'source', repos = 'https://cloud.r-project.org')" \
    && R CMD build . \
    && R CMD INSTALL --build *.tar.gz \
    && rm -Rf /tmp/downloaded_packages/ /tmp/*.rds
