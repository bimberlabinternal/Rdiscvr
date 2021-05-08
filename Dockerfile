from ghcr.io/bimberlab/cellhashr:latest

ADD . /RDiscvr

ENV RETICULATE_PYTHON=/usr/bin/python3

RUN cd /RDiscvr \
	&& R CMD build . \
	&& Rscript -e "BiocManager::install(ask = F, upgrade = 'always');" \
    # NOTE: related to: https://github.com/satijalab/seurat/issues/4436. Should remove this once Matrix issue is fixed.
    && Rscript -e "devtools::install_version('Matrix', version = '1.3-2', dependencies=TRUE, ask = FALSE)" \
	&& Rscript -e "devtools::install_deps(pkg = '.', dependencies = TRUE, upgrade = 'always');" \
	&& R CMD INSTALL --build *.tar.gz \
	&& rm -Rf /tmp/downloaded_packages/ /tmp/*.rds