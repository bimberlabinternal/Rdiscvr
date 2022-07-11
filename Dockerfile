from ghcr.io/bimberlab/cellhashr:latest

ADD . /RDiscvr

ENV RETICULATE_PYTHON=/usr/bin/python3

# NOTE: secret is used to pass a github token to avoid GitHub API rate limit issues
RUN --mount=type=secret,id=GITHUB_PAT \
    cd /RDiscvr \
    && export GITHUB_PAT="$(cat /run/secrets/GITHUB_PAT)" \
    && echo "GH: $GITHUB_PAT" \
	&& R CMD build . \
	&& Rscript -e "BiocManager::install(ask = F, upgrade = 'always');" \
	&& Rscript -e "devtools::install_deps(pkg = '.', dependencies = TRUE, upgrade = 'always');" \
	&& R CMD INSTALL --build *.tar.gz \
	&& rm -Rf /tmp/downloaded_packages/ /tmp/*.rds
