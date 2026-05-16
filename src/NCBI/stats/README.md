# NCBI Analysis: File Extensions and SBML Models

Analysis script for examining file types in NCBI article downloads and identifying computational models (SBML files).

## Overview

This R script processes NCBI article downloads across multiple disease categories and generates:

1. A visualization of file type distributions showing that NCBI downloads contain mostly metadata (PDFs, images, documents) rather than computational models
2. A detailed list of SBML organisms found in a specific article (PMC6620235) that contains 20 SBML files

## Requirements

- R 4.0 or later
- tidyverse
- ggplot2

Install dependencies:
```r
install.packages("tidyverse")
```

## Usage

```bash
Rscript ncbi_analysis.R
```

The script will read from:
- `downloaded_ncbi_models 2/` - subdirectories containing disease summary CSVs
- `data_manifest.csv` - manifest of all extracted files

Output files are saved to `ncbi_analysis_outputs/`:
- `01_file_extensions.png` - Bar chart of top 15 file types
- `02_pmc6620235_sbml_organisms.csv` - Table of SBML organisms in PMC6620235

## Data Structure

### Input: disease_summary.csv
One file per disease category (influenza_files, tuberculosis_files, covid_files, etc.)

Columns include counts for different file types:
- data_xml, data_json, data_csv, etc.
- metadata_pdf, metadata_png, metadata_docx, etc.

### Input: data_manifest.csv
Complete file manifest with columns:
- pmcid: PubMed Central article ID
- filename: Name of the file
- category: File type (excel, sbml_xml, pdf, etc.)
- source: Original file path

### Output: 01_file_extensions.png
Horizontal bar chart showing count of each file type found across all downloads. Types are colored by category (DATA vs METADATA).

### Output: 02_pmc6620235_sbml_organisms.csv
CSV with three columns:
- No.: Row number
- Organism: Organism or model name extracted from SBML filename
- Type: Classification (Tuberculosis, Influenza, or Other)

## Key Finding

NCBI article downloads contain mostly metadata files (PDFs, images, documents). Computational SBML models are rare. PMC6620235 is an example article that contains 20 SBML files, including models for tuberculosis and influenza organisms.

## Script Details

The script:
1. Reads disease_summary.csv files from each disease folder
2. Aggregates file extension counts across all diseases
3. Creates a visualization of the top 15 file types
4. Extracts SBML organism data from the manifest for PMC6620235
5. Identifies organism type based on filename pattern matching
6. Exports results as PNG and CSV files

## Files

- `ncbi_analysis.R` - Main analysis script
- `README.md` - This file
