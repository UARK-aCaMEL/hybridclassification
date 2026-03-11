FROM rocker/r-ver:4.5.2

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates curl git \
    build-essential gfortran pkg-config cmake \
    libcurl4-openssl-dev libssl-dev libxml2-dev zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Use Posit Package Manager (more likely to give you binaries on Debian)
ENV R_REPO="https://packagemanager.posit.co/cran/__linux__/bookworm/latest"

# Install dependencies (Rcpp/Rstan toolchain etc.)
RUN Rscript -e 'options(repos=c(CRAN=Sys.getenv("R_REPO")), download.file.method="libcurl"); \
  install.packages(c("Rcpp","BH","RcppEigen","RcppParallel","StanHeaders","rstan","rstantools"))'

# Install bgc-hm without GitHub API calls
RUN git clone --depth 1 https://github.com/zgompert/bgc-hm.git /tmp/bgc-hm \
 && R CMD INSTALL /tmp/bgc-hm \
 && Rscript -e 'library(bgchm); cat("bgchm loaded OK\n")'

WORKDIR /work
