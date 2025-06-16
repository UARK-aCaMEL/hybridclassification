# Use official R environment
FROM rocker/r-base:4.2.3

# Install build tools, Git, and system libraries needed for R packages
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      build-essential \
      gfortran \
      git \
      libcurl4-gnutls-dev \
      libssl-dev \
      libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Clone & compile NewHybrids (no-GUI)
RUN git clone https://github.com/eriqande/newhybrids.git /opt/newhybrids && \
    cd /opt/newhybrids && \
    git submodule update --init --recursive && \
    chmod +x Compile-with-no-gui-Linux.sh && \
    ./Compile-with-no-gui-Linux.sh && \
    mv newhybrids-no-gui-linux.exe /usr/local/bin/newhybrids && \
    ln -s /usr/local/bin/newhybrids /usr/local/bin/newhybs && \
    chmod +x /usr/local/bin/newhybrids

# Install required R packages from CRAN
RUN Rscript -e "install.packages(c('adegenet','vcfR','plyr','stringr','tidyr','dartR'), repos='https://cloud.r-project.org')"

# Install remotes and GitHub packages
RUN Rscript -e "install.packages('remotes', repos='https://cloud.r-project.org'); \
                remotes::install_github(c('bwringe/parallelnewhybrid', 'bwringe/hybriddetective', 'rystanley/genepopedit'))"

# Set working directory for Nextflow runs
WORKDIR /data

# Default to bash
CMD ["bash"]
