FROM ghcr.io/bimberlabinternal/discvr-base:latest

ADD . /RDiscvr

# NOTE: secret is used to pass a github token to avoid GitHub API rate limit issues
RUN --mount=type=secret,id=GITHUB_PAT \
    cd /RDiscvr \
    && export GITHUB_PAT="$(cat /run/secrets/GITHUB_PAT)" \
    && Rscript -e "BiocManager::install(ask = F, upgrade = 'always');" \
    && Rscript -e "devtools::install_deps(pkg = '.', dependencies = TRUE, upgrade = 'never');" \
    # This should ultimately be removed. It is being used to fix the merge() bug in SeuratObject 5.3.0
    && Rscript -e "remotes::install_version('SeuratObject', version = '5.2.0')" \
    && R CMD build . \
    && R CMD INSTALL --build *.tar.gz \
    && rm -Rf /tmp/downloaded_packages/ /tmp/*.rds
