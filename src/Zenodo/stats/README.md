# Pathogen Model Analytics & Visualizer

This R script is a data processing and analytics tool designed to parse, aggregate, and visualize downloaded modeling assets from the Zenodo extraction pipeline. It recursively scans directory structures, counts computational modeling files based on their extensions, reads localized metadata files, and automatically exports clean statistical plots.

---

## Features

* **Automated Data Harvesting**: Recursively reads nested directories on disk without needing a central database to rebuild dataset summaries dynamically.
* **Smart Context Parsing**: Automatically extracts metadata records and determines target pathogen categories directly from the file system's folder hierarchy.
* **Extension Frequency Tracking**: Computes exact counts for target script types: Python (`py`, `ipynb`), R (`r`), MATLAB (`m`), Julia (`jl`), C/C++ (`c`, `cpp`), and systems biology modeling files (`sbml`, `sedml`, `omex`, `cps`, `cellml`).
* **Visual Export Suite**: Pre-configured thematic engine utilizing `ggplot2` that exports production-ready analytics charts instantly.

---

## Directory Context Requirement

For the script to work properly, it must be placed and executed in the parent directory where the data extraction pipeline stored its outputs. It assumes the following file tree layout:

```text
. (Current Working Directory)
└── [Pathogen_Name]/
    └── Zenodo/
        └── DOI[Extracted_ID]/
            ├── data/
            │   └── [Files with targeted extensions]
            └── metadata/
                └── metadata_[Extracted_ID].json
