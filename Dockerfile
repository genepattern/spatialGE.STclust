FROM bioconductor/bioconductor:3.18-R-4.3.1

RUN \
  apt-get update && \
  apt-get install -y --no-install-recommends libxml2-dev libcurl4-openssl-dev libssl-dev && \
  Rscript -e "install.packages('optparse')" && \
#   Rscript -e "devtools::install_github('fridleylab/spatialGE')" && \
  ## hand edited
  Rscript -e "install.packages('msigdbr')" && \
  Rscript -e "install.packages('tidyverse')"


RUN Rscript -e "remotes::install_github('fridleylab/spatialGE')"