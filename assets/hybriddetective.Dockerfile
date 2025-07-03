###############################################################################
# hybridclassification / hybriddetective debug image
###############################################################################

FROM rocker/r-base:4.2.3

# ───── system libraries ─────────────────────────────────────────────────────────
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      build-essential gfortran git procps \
      libcurl4-gnutls-dev libssl-dev libxml2-dev libgit2-dev \
      libpng-dev libjpeg-dev libtiff5-dev libfreetype6-dev \
      libharfbuzz-dev libfribidi-dev libfontconfig1-dev \
      libgdal-dev libproj-dev libudunits2-dev libgsl-dev \
      libgeos-dev libsqlite3-dev ca-certificates \
      autoconf automake libtool pkg-config file m4 perl \
    && rm -rf /var/lib/apt/lists/*

# ───── shared site-library for every UID ────────────────────────────────────────
ENV R_LIBS_SITE=/usr/local/lib/R/site-library
ENV R_LIBS_USER=$R_LIBS_SITE

RUN mkdir -p "$R_LIBS_SITE" && chmod -R a+rX "$R_LIBS_SITE" && \
    mkdir -p /usr/lib/R/etc && \
    echo '.libPaths(c("/usr/local/lib/R/site-library", .libPaths()))' \
        >> /usr/lib/R/etc/Rprofile.site

# ───── NewHybrids binary (no GUI) ───────────────────────────────────────────────
RUN git clone https://github.com/eriqande/newhybrids.git /opt/newhybrids && \
    cd /opt/newhybrids && \
    git submodule update --init --recursive && \
    chmod +x Compile-with-no-gui-Linux.sh && \
    ./Compile-with-no-gui-Linux.sh && \
    mv newhybrids-no-gui-linux.exe /usr/local/bin/newhybrids && \
    ln -s /usr/local/bin/newhybrids /usr/local/bin/newhybs && \
    chmod +x /usr/local/bin/newhybrids

# ───── CRAN packages ────────────────────────────────────────────────────────────
RUN Rscript -e "install.packages(c( \
      'optparse','adegenet','vcfR','plyr','stringr','tidyr','dartR', \
      'dplyr','ggplot2','ggrepel','ggridges','ggpubr','cowplot','randomcoloR', \
      'reshape2','snowfall','data.table','VGAM','Rmisc','gridExtra', \
      'sf','rgdal','raster'), \
      repos='https://cloud.r-project.org', dependencies=TRUE)"

# ───── GitHub packages (visible output, log captured) ──────────────────────────
RUN Rscript -e "install.packages('remotes', repos='https://cloud.r-project.org')"
RUN Rscript -e "\
  remotes::install_github('bwringe/parallelnewhybrid', \
      dependencies = TRUE, upgrade = 'never', build_vignettes = FALSE, force = TRUE)"

# ───── genepopedit with full log capture ────────────────────────────────────────────
RUN set -e; \
    echo '::: installing genepopedit :::'; \
    Rscript -e "\
      tryCatch( \
        remotes::install_github('rystanley/genepopedit', \
            dependencies=TRUE, upgrade='never', build_vignettes=FALSE, force=TRUE), \
        error=function(e){ \
          writeLines(conditionMessage(e), '/tmp/genepopedit.log'); \
          stop(e) \
        })" \
    || { echo '==== genepopedit build log ===='; cat /tmp/genepopedit.log || true; exit 1; }; \
    test -d /usr/local/lib/R/site-library/genepopedit

# ───── sanity check ─────────────────────────────────────────────────────────────────
RUN Rscript -e "stopifnot(requireNamespace('genepopedit', quietly=TRUE))"

# ---- install hybriddetective with full log ----
RUN set -e; \
    echo '::: installing hybriddetective :::'; \
    Rscript -e "\
      tryCatch( \
        remotes::install_github('bwringe/hybriddetective', \
            dependencies = TRUE, upgrade = 'never', build_vignettes = FALSE, \
            force = TRUE, keep_outputs = TRUE), \
        error = function(e) { \
            writeLines(conditionMessage(e), '/tmp/hybriddetective.log'); \
            stop(e) \
        })" \
    || { echo '==== build log ===='; cat /tmp/hybriddetective.log || true; exit 1; }; \
    test -d /usr/local/lib/R/site-library/hybriddetective


# ───── sanity check – must load for root at build time ─────────────────────────
RUN Rscript -e "stopifnot(requireNamespace('hybriddetective', quietly = TRUE))"

WORKDIR /data
CMD [\"bash\"]

