# InfectioGIT Docker Image
# Multi-source infectious disease data analysis platform
# Supports: BioModels, NCBI, Zenodo, SQLite database

FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    wget \
    git \
    ca-certificates \
    r-base \
    r-base-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements first for better Docker layer caching
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy entire project
COPY . .

# Set environment variables
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

# Create output directories
RUN mkdir -p /app/outputs \
    && mkdir -p /app/BioModels_Database_Final \
    && mkdir -p /app/BioModels_Stats \
    && mkdir -p /app/downloaded_ncbi_models

# Default command: interactive Python shell
# To run specific scripts:
#   docker run infectio python src/Biomodels/biomodels_pipeline.py
#   docker run infectio python src/NCBI/ncbi_extraction.py
#   docker run infectio python src/Database/database_creation.py
CMD ["python", "-i"]
