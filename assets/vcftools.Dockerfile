FROM debian:bookworm-slim

ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies (including htslib headers for cyvcf2) and procps
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      build-essential \
      gcc \
      curl \
      procps \
      ca-certificates \
      zlib1g-dev \
      libbz2-dev \
      liblzma-dev \
      libcurl4-gnutls-dev \
      libssl-dev \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        bash \
        build-essential \
        gcc \
        curl \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-gnutls-dev \
        libssl-dev \
        coreutils \
        gawk \
        grep \
        sed \
        procps \
        vcftools \
        tabix \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app
