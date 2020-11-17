### copyright 2017-2020 Regents of the University of California and the Broad Institute. All rights reserved.

FROM satijalab/seurat:3.2.0
MAINTAINER Edwin Juarez <ejuarez@ucsd.edu>

ENV LANG=C LC_ALL=C
USER root

RUN R -e "library('Seurat');sessionInfo()"
RUN R -e 'install.packages("BiocManager",repos = "http://cran.us.r-project.org")'

RUN R -e 'requireNamespace("BiocManager", quietly = TRUE); BiocManager::install("scater")'
RUN R -e "library('scater');sessionInfo()"

RUN R -e 'install.packages("optparse",repos = "http://cran.us.r-project.org")'
RUN R -e "library('optparse')"

RUN apt-get update && apt-get install -y time

RUN R -e 'install.packages("log4r", repos = "http://cran.us.r-project.org")'

RUN R -e 'requireNamespace("BiocManager", quietly = TRUE); BiocManager::install("Signac")'
RUN R -e "library('Signac');sessionInfo()"

RUN mkdir /module
ADD sample_data /module/sample_data

RUN mkdir /temp
COPY seurat_preprocess.R /module/

# build using this:
# docker build -t genepattern/seurat-suite:2.4 .
