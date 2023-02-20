from ghcr.io/bimberlab/cellhashr:latest

ADD . /RDiscvr

ENV RETICULATE_PYTHON=/usr/bin/python3

# NOTE: ggplot2 added to force version 3.4.0, which is needed by ggtree. Otherwise this container is pegged to ./focal/2022-10-28
RUN echo "local({r <- getOption('repos') ;r['CRAN'] = 'https://packagemanager.rstudio.com/cran/__linux__/focal/latest';options(repos = r);rm(r)})" >> ~/.Rprofile

# NOTE: secret is used to pass a github token to avoid GitHub API rate limit issues
RUN --mount=type=secret,id=GITHUB_PAT \
    cd /RDiscvr \
    && export GITHUB_PAT="$(cat /run/secrets/GITHUB_PAT)" \
    && echo "GH: $GITHUB_PAT" \
	&& R CMD build . \
	# NOTE: Matrix install has been added to force 1.5.3. Can be removed after main repos updated. \
    && Rscript -e "install.packages('Matrix', repos = 'https://packagemanager.posit.co/cran/__linux__/focal/latest', update = TRUE, ask = FALSE)" \
	&& Rscript -e "BiocManager::install(ask = F, upgrade = 'always');" \
	&& Rscript -e "devtools::install_deps(pkg = '.', dependencies = TRUE, upgrade = 'always');" \
	&& R CMD INSTALL --build *.tar.gz \
	&& rm -Rf /tmp/downloaded_packages/ /tmp/*.rds
