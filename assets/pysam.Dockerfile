# Use official Python 3.12 slim image as base
FROM python:3.12-slim

# Avoid Python writing .pyc files and buffer issues
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

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

# Set working directory
WORKDIR /app

# Install Python packages
RUN pip install --no-cache-dir \
      cyvcf2 \
      numpy \
      pandas \
      scipy \
      pysam

