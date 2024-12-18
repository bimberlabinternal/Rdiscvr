FROM ghcr.io/bimberlabinternal/discvr-base:latest

ADD . /RDiscvr

# NOTE: secret is used to pass a github token to avoid GitHub API rate limit issues
RUN --mount=type=secret,id=GITHUB_PAT \
    cd /RDiscvr \
    && export GITHUB_PAT="$(cat /run/secrets/GITHUB_PAT)" \
    && Rscript -e "BiocManager::install(ask = F, upgrade = 'always');" \
    # TODO: added to avoid fgsea installation issues:
    && Rscript -e "remotes::install_version('BH', '1.84.0-0');" \
    && Rscript -e "devtools::install_deps(pkg = '.', dependencies = TRUE, upgrade = 'never');" \
    && R CMD build . \
    && R CMD INSTALL --build *.tar.gz \
    && rm -Rf /tmp/downloaded_packages/ /tmp/*.rds
