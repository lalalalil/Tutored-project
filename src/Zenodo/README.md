# Zenodo pipeline extractor

## Overview

This project provides a complete Python pipeline for downloading, organizing, filtering, and preparing computational models of infectious diseases from the Zenodo repository.

The pipeline automatically connects to the Zenodo API to search for target pathogens while filtering out general epidemiological noise. It systematically separates modeling scripts from metadata, creates structured local directories, isolates valid executable files, and eliminates empty repositories.

The script is specifically designed for large-scale comparative analysis of infectious disease computational models, serving as the data ingestion step before downstream statistical analysis.

---

# Features

* **Automated Zenodo API Crawling**: Performs multi-page queries using specific pathogen keywords and automated pagination.
* **Targeted Noise Filtering**: Uses negative API queries to exclude unrelated articles focusing on epidemics, climate, transmission, or hospital bed capacities.
* **Format-Specific Extraction**: Detects and isolates specific computational files including Python, R, MATLAB, Julia, C/C++, and systems biology standards.
* **Hierarchical Repository Storage**: Automatically builds organized subdirectories for data and metadata based on extracted clean DOIs.
* **Self-Cleaning Architecture**: Automatically removes temporary or empty virus directories that do not contain valid executable files.
* **Rate-Limit Safeguards**: Integrates strategic pauses to comply with Zenodo API limitations and ensure data transfer stability.

---

# Supported Diseases

The pipeline currently searches models related to:

* Dengue (DENV)
* Chikungunya (CHIKV)
* Mpox (monkeypox)
* West Nile Virus (WNV)
* Influenza (including avian influenza and H5N1 strains)
* Tuberculosis (mycobacterium)
* HIV (human immunodeficiency virus)
* COVID-19 (SARS-CoV-2)

---

# Project Structure

The extraction pipeline dynamically builds and manages the following file system layout:

```text
├──Models
└── [Pathogen_Name]/
    └── Zenodo/
        └── DOI[Extracted_ID]/
            ├── data/
            │   └── model_file.py
            └── metadata/
                └── metadata_[Extracted_ID].json
