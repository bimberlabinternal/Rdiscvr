from ghcr.io/bimberlab/cellhashr:latest

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
    && pip3 install umap-learn phate scanpy[leiden] \
    && mkdir /conga \
    && cd /conga \
    && git clone -b rhesus2 https://github.com/bbimber/conga.git \
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
    # TODO: drop this once main CRAN repo contains version 0.2.1:
    && Rscript -e "install.packages('aplot', repos = 'https://cran.wustl.edu/')" \
	&& Rscript -e "devtools::install_deps(pkg = '.', dependencies = TRUE, upgrade = 'always');" \
    # NOTE: Related to: https://github.com/satijalab/seurat/issues/7328. Should revert to a release once patched.
    && Rscript -e "remotes::install_github('satijalab/seurat', ref='443ab86684253d9a7290c3d38c2bc1d8db021776');" \
    && R CMD build . \
	&& R CMD INSTALL --build *.tar.gz \
	&& rm -Rf /tmp/downloaded_packages/ /tmp/*.rds
