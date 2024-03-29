# Dockerfile with development tools for:
# easy testing, documentation generation and benchmarking
FROM rocker/r-ver:4.1.1

RUN apt-get update && apt-get -y install \
  texlive-latex-base \
  texlive-fonts-extra \
  texinfo \
  libcurl4-gnutls-dev \
  libssl-dev \
  libxml2-dev \
  libgit2-dev \
  pandoc \
  pandoc-citeproc \
  valgrind \
  time

RUN R -e "install.packages(c('devtools', 'roxygen2', 'testthat'))"

COPY DESCRIPTION /home/docker/msio/

RUN R -e "library('remotes'); install_deps('/home/docker/msio', dependencies = TRUE)"
