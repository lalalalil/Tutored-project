# InfectioGIT: A Repository of Computational Models and Metadata for Infectious Diseases

## Description

InfectioGIT is a project dedicated to the collection, organization, standardization, and analysis of computational models related to infectious diseases. The project integrates models and associated metadata retrieved from multiple scientific repositories in order to create a centralized, reusable, and FAIR-compliant resource.

The repository gathers models from several public scientific platforms, including BioModels, Zenodo, and National Center for Biotechnology Information.

InfectioGIT was developed within the broader context of the *Digital Twin* initiative, aiming to contribute to the development of infectious disease and immune-system-oriented digital twins.

---

# Table of Contents

* [Project Objectives](#project-objectives)
* [Key Features](#key-features)
* [Repositories and Data Sources](#repositories-and-data-sources)
* [Supported Diseases](#supported-diseases)
* [Pipeline Overview](#pipeline-overview)
* [Folder Structure](#folder-structure)
* [Metadata and FAIR Principles](#metadata-and-fair-principles)
* [Statistical Analysis](#statistical-analysis)
* [Installation and Usage](#installation-and-usage)
* [Outputs](#outputs)
* [Contributing](#contributing)
* [Documentation](#documentation)
* [License](#license)
* [Credits](#credits)

---

# Project Objectives

The main objectives of InfectioGIT are:

* Collect computational models related to infectious diseases
* Centralize models from multiple scientific repositories
* Standardize metadata and repository structures
* Improve accessibility and reusability of modelling resources
* Apply FAIR principles to infectious disease computational models
* Facilitate comparative analyses across diseases and repositories
* Support future infectious disease digital twin initiatives

---

# Key Features

## Automated Model Collection

Python pipelines were developed to automatically retrieve models and metadata from different repositories.

### Features

* repository querying
* automated downloads
* metadata extraction
* DOI retrieval
* repository organization

---

## Metadata Standardization

Metadata are either:

* directly retrieved from repositories
* automatically generated and enriched

The project centralizes metadata into structured JSON and CSV formats.

---

## Repository Organization

Models are classified according to:

* disease category
* source repository
* modelling approach
* curation status
* publication year

---

## Statistical Analysis

Dedicated R scripts generate:

* publication trend analyses
* modelling approach distributions
* repository statistics
* file format analyses
* curation quality visualizations

---

## FAIR and Reproducibility

To improve reproducibility and accessibility:

* pipelines were containerized using Docker
* metadata were standardized
* repository structures were normalized
* models and associated resources were centralized in GitHub

---

# Repositories and Data Sources

The project integrates resources from:

* BioModels
* Zenodo
* National Center for Biotechnology Information

Each repository provides complementary types of resources and modelling files.

---

# Supported Diseases

The repository currently includes models related to:

* COVID-19
* HIV
* Tuberculosis
* Influenza
* Dengue
* Chikungunya
* Mpox
* West Nile Virus

---

# Pipeline Overview

The global workflow followed during the project is summarized below:

```text id="5r2qja"
Scientific Literature Review
            ↓
Selection of Priority Infectious Diseases
            ↓
Automatic Querying of Scientific Repositories
(BioModels / Zenodo / NCBI)
            ↓
Model and Metadata Retrieval
            ↓
Metadata Enrichment and Standardization
            ↓
Classification and Filtering
            ↓
Statistical Analysis and Visualization
            ↓
Integration into InfectioGIT Repository
```

---

# Folder Structure

```text id="8fww8z"
InfectioGIT/
├── Results/
│   ├── COVID/
│   │   ├── BioModels/
│   │       ├── metadata
│   │       └── model
│   │   ├── Zenodo/
│   │   └── NCBI/
│   │
│   ├── HIV/
│   ├── Tuberculosis/
│   ├── Influenza/
│   ├── Dengue/
│   ├── Chikungunya/
│   ├── Mpox/
│   ├── West_Nile/
│   │
│   ├── Statistics/
│   └── Outputs/
│
├── scripts/
│   ├── Python/
│   └── R/
│
├── docker/
├── docs/
└── README.md
```

---

# Metadata and FAIR Principles

The project follows FAIR principles:

* **Findable**
* **Accessible**
* **Interoperable**
* **Reusable**

To support these objectives:

* metadata are standardized
* DOI information is preserved
* repository structures are harmonized
* reproducible workflows are provided through Docker

---

# Statistical Analysis

The repository includes statistical analyses such as:

* publication trends over time
* modelling approach distributions
* disease-specific repository comparisons
* file extension analyses
* curation quality evaluations

Generated figures are exported in publication-ready PNG format.

---

# Installation and Usage

## Quick Start with Docker (Recommended)

Docker provides a complete, reproducible environment with all dependencies pre-configured.

### Prerequisites
- [Docker Desktop](https://www.docker.com/products/docker-desktop)

### Setup

```bash
# Navigate to docker folder
cd docker

# Build the Docker image (takes ~3 minutes)
docker-compose build

# Verify installation
docker-compose run infectio python -c "
from Bio import Entrez
import bioservices
print('✓ All dependencies installed successfully!')
"
```

### Running Pipelines with Docker

#### BioModels Analysis Pipeline
Interactive 8-step pipeline for downloading and analyzing BioModels:
```bash
cd docker
docker-compose run infectio python src/Biomodels/biomodels_pipeline.py
```

**Menu Options:**
1. Download & Sort (metadata vs model files)
2. Clean Empty Models & Generate Extension Stats
3. Separate by Curation Status
4. Classify by Modelling Approach
5. Generate Publication Dates CSV
6. Delete Models Published Before 2015
7. Clean Metadata Directories
8. Run R Statistical Analysis
9. Exit

#### NCBI PubMed Central Extraction
Download open-access articles and code from PubMed Central:
```bash
cd docker
docker-compose run infectio python src/NCBI/ncbi_extraction.py
```

#### Database Creation
Scan GitHub repository and populate SQLite database:
```bash
cd docker
docker-compose run infectio python src/Database/database_creation.py
```

### Interactive Python Shell

```bash
cd docker
docker-compose run infectio python -i
```

Then import and use the modules:
```python
from Bio import Entrez
import bioservices
import pandas as pd
# ... your code here
```

---

## Local Installation (Alternative)

### Python Requirements

```bash
pip install biopython bioservices biomodels requests pandas numpy scipy xmltodict pyyaml tqdm
```

### R Requirements

```r
install.packages(c(
  "tidyverse",
  "janitor",
  "RColorBrewer"
))
```

### Running the Pipelines

```bash
# BioModels Pipeline
python src/Biomodels/biomodels_pipeline.py

# NCBI Extraction
python src/NCBI/ncbi_extraction.py

# Database Creation
python src/Database/database_creation.py
```

---

# Outputs

The project generates:

* organized model repositories
* standardized metadata
* CSV statistical summaries
* publication-ready visualizations
* comparative repository analyses

Example outputs:

```text id="az96hk"
publication_dates_summary.csv
modelling_approaches_summary.csv
summary_avg_per_category.csv
01_publications_global.png
05_approaches_by_disease.png
```

---

# Contributing

Contributions are welcome for:

* new disease categories
* additional repositories
* metadata enrichment
* FAIR standardization
* pipeline optimization
* statistical analyses

---

# Documentation

Additional documentation includes:

* pipeline descriptions
* repository organization guides
* statistical analysis documentation
* usage examples
* reproducibility instructions

---

# License

This project is distributed under the MIT License.

---

# Credits

## Supervisors

* Anna Niarakis
* Virginie Jouffret

## Students

* Aya Sebbah
* Diana Bravais
* Thalie Holmiere
* [Additional contributors]

---

# Affiliations

* University of Toulouse
* Centre de Biologie Intégrative Toulouse (CBI)
* Computational Systems Biology teams
* Digital Twin research initiatives

---

# References

1. [BioModels](https://www.ebi.ac.uk/biomodels/?utm_source=chatgpt.com)
2. [Zenodo](https://zenodo.org/?utm_source=chatgpt.com)
3. [NCBI](https://www.ncbi.nlm.nih.gov/?utm_source=chatgpt.com)
