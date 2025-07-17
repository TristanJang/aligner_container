# syntax=docker/dockerfile:1

# Use an official lightweight Python image as base
FROM python:3.10-slim

# Set metadata
LABEL maintainer="Tristan Jang"
LABEL description="Standalone sequencing aligner and reporting system for SOPHiA GENETICS assignment"

# Avoid interactive prompts during installs
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies (bwa, samtools) and cleanup to minimize image size
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        bwa \
        samtools \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Set the working directory
WORKDIR /app

# Copy Python dependency definitions first (for Docker layer caching)
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy remaining source code into the image
COPY . .

# Use unbuffered output (useful for logs)
ENV PYTHONUNBUFFERED=1

# Default command (can be overridden)
CMD ["bash"]
