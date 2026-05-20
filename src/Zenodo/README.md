# Zenodo Pathogen Model Extractor

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
└── [Pathogen_Name]/
    └── Zenodo/
        └── DOI[Extracted_ID]/
            ├── data/
            │   └── model_file.py
            └── metadata/
                └── metadata_[Extracted_ID].json
```

---

# Installation

## Requirements

* Python 3.x

## Python Dependencies

Install the required network package:

```bash
pip install requests
```

---

# Usage

Run the extraction script directly from your terminal:

```bash
python zenodo_extraction.py
```

---

# Pipeline Steps

## 1. Automated Search & Pagination
The script queries the Zenodo API with a targeted keyword payload. It loops through up to five pages per query variant to gather relevant software entries while actively filtering out duplicate versions.

## 2. Extension Filtering
Every retrieved record undergoes a strict validation check. The pipeline inspects the file keys inside the metadata and ensures they match the target computational biology and programming formats before initiating any download.

## 3. Directory Creation & Metadata Isolation
For every valid record, a clean numeric identifier is extracted from the DOI. The script generates specialized data and metadata folders, writing the complete, unedited Zenodo JSON response directly into the metadata directory.

## 4. Source File Download
The script connects to the specific download links provided by the API context, downloading the raw script assets directly into the localized data subdirectories using binary writing modes.

## 5. Post-Extraction Clean Up
Once all keyword variants for a specific pathogen are processed, the system evaluates the target folder contents. If a pathogen directory contains no valid downloaded model files, the entire directory tree is wiped to keep the workspace clean.

---

# Metadata Structure

Example of the structured metadata JSON saved by the script:

```json
{
    "id": 1234567,
    "links": {
        "doi": "[https://doi.org/10.5281/zenodo.1234567](https://doi.org/10.5281/zenodo.1234567)",
        "self": "[https://zenodo.org/api/records/1234567](https://zenodo.org/api/records/1234567)"
    },
    "metadata": {
        "title": "Pathogen Computational Model",
        "publication_date": "2026-05-20",
        "description": "Source code repository for target disease dynamics."
    },
    "files": [
        {
            "key": "simulation.py",
            "links": {
                "self": "[https://zenodo.org/api/files/.../simulation.py](https://zenodo.org/api/files/.../simulation.py)"
            }
        }
    ]
}
```

---

# Recognized File Types

The script isolates computational modeling files and code scripts using a synchronized extension list:

* **Python**: `.py`, `.ipynb`
* **R Script**: `.r`
* **MATLAB**: `.m`
* **Julia**: `.jl`
* **C / C++**: `.c`, `.cpp`
* **Systems Biology / Copasi**: `.sbml`, `.sedml`, `.omex`, `.cps`, `.cellml`

---
