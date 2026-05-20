# Pathogen Model Extractor & Analyzer

This project is a two-part pipeline designed to systematically retrieve, filter, and analyze mathematical and computational models of infectious diseases from the Zenodo repository. 

1. **Python Script**: Queries the Zenodo API for specific pathogens, filters for software records containing executable or modeling file extensions, and downloads both the metadata and the source files.
2. **R Script**: Parses the downloaded local JSON metadata and files to generate statistical insights and structured visualizations.

---

## Features

* **Targeted Filtering**: Automatically excludes general epidemiological noise (e.g., climate, bed capacity, transmission rates) using negative API queries to focus strictly on structural or computational models.
* **Format-Specific Extraction**: Detects and downloads specific formats: Python (`.py`, `.ipynb`), R (`.r`), MATLAB (`.m`), Julia (`.jl`), C/C++ (`.c`, `.cpp`), and systems biology standards (`.sbml`, `.sedml`, `.omex`, `.cps`, `.cellml`).
* **Automated Architecture**: Structures downloaded files hierarchically by target pathogen, creating organized `metadata` and `data` subdirectories.
* **Rate-Limit Friendly**: Embedded pauses (`time.sleep`) prevent API overload and connection drops.

---

## Directory Structure

The Python script builds the following architectur. The R script relies on this layout to aggregate statistics:

```text
Models/
└── [Pathogen_Name]/
    └── Zenodo/
        └── DOI[Extracted_ID]/
            ├── data/
            │   └── model_file.py
            └── metadata/
                └── metadata_[Extracted_ID].json
