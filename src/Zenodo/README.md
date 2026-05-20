# Pathogen Model Extractor (Zenodo API)

This Python script runs a data extraction pipeline designed to systematically retrieve and filter mathematical and computational models of infectious diseases from the Zenodo repository. It queries the Zenodo API, applies strict metadata filters, isolates records containing executable model code, and downloads both the raw files and their structured metadata.

---

## Features

* **Targeted Filtering**: Automatically excludes general epidemiological noise (e.g., climate, bed capacity, transmission rates) using negative API queries to focus strictly on structural or computational models.
* **Format-Specific Extraction**: Detects and downloads specific formats: Python (`.py`, `.ipynb`), R (`.r`), MATLAB (`.m`), Julia (`.jl`), C/C++ (`.c`, `.cpp`), and systems biology standards (`.sbml`, `.sedml`, `.omex`, `.cps`, `.cellml`).
* **Automated Architecture**: Structures downloaded files hierarchically by target pathogen, creating organized `metadata` and `data` subdirectories.
* **Rate-Limit Friendly**: Embedded pauses (`time.sleep`) prevent API overload, connection drops, or rate-limiting blocks.

---

## Directory Structure

The script builds the following architecture. This exact layout is required for subsequent downstream analysis scripts to run properly:

```text
Models/
└── [Pathogen_Name]/
    └── Zenodo/
        └── DOI[Extracted_ID]/
            ├── data/
            │   └── model_file.py
            └── metadata/
                └── metadata_[Extracted_ID].json
