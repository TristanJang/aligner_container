# syntax=docker/dockerfile:1

# Use official lightweight Python image as base
FROM python:3.10-slim

# Install bwa and samtools
RUN apt-get update && \
    apt-get install -y bwa samtools time && \
    apt-get clean

# Set working directory
WORKDIR /app

# Copy your project files into the container
COPY . .

# Set default command (can be overridden)
CMD ["bash"]
