# Bioconductor devel with annotatr's dependencies preinstalled, for docker/check.sh.
# Rebuild when Bioconductor devel moves to a new version (after each release), and
# periodically within a cycle so dependencies match the Bioconductor builders:
#   docker/check.sh --rebuild
#
# The image comes from GitHub's registry; the Docker Hub copy is
# bioconductor/bioconductor_docker:devel.
ARG BIOC_IMAGE=ghcr.io/bioconductor/bioconductor:devel
FROM ${BIOC_IMAGE}

COPY docker/install_deps.R /opt/annotatr/install_deps.R
COPY DESCRIPTION /opt/annotatr/DESCRIPTION
RUN Rscript /opt/annotatr/install_deps.R /opt/annotatr/DESCRIPTION
