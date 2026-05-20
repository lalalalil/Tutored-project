# InfectioGIT Database Creator

## Overview

This project provides an automated Python pipeline designed to build a comprehensive SQLite database (`infectio_git.db`) from biological computational models. 

The script dynamically fetches JSON metadata files directly from the `thalieh/InfectioGIT` GitHub repository. It parses complex metadata from BioModels, Zenodo, and NCBI to extract and infer critical relationships between mathematical models, infectious diseases, pathogens, and biological hosts. The result is a fully structured relational database ready for complex SQL queries and Pandas data analysis.

---

# Features

* **Automated GitHub Scraping**: Uses the GitHub API to recursively locate and download JSON metadata files without manual local cloning.
* **Multi-Source Parsing**: Intelligently detects and processes varying JSON structures from BioModels, Zenodo, and NCBI.
* **Semantic Inference**: Automatically deduces disease names, pathogen classifications, and biological scales from free text and titles when formal ontology annotations are missing.
* **Relational SQLite Schema**: Automatically generates a complex schema including tables for Models, Diseases, Pathogens, Hosts, Tissues, Cell Types, and Datasets with strict foreign key constraints.
* **Data Integration & Analytics**: Built-in Pandas integration to instantly verify database integrity and perform multi-table SQL joins.

---

# Supported Diseases

The pipeline is pre-configured to infer and map the following targets to their respective DOID and taxonomy IDs:

* COVID-19 (SARS-CoV-2)
* Influenza (Influenza A virus)
* Tuberculosis (Mycobacterium tuberculosis)
* Dengue (Dengue virus)
* HIV
* Zika
* Chikungunya

---

# Database Schema

The script generates the following primary entities:

* `Model`: Stores DOIs, titles, mathematical formalisms, scales, years, and software.
* `Disease`: Maps to DOID/Mondo ontologies.
* `Pathogen`: Stores NCBI Taxonomy IDs and pathogen classes.
* `Host`: Stores NCBI Taxonomy IDs (e.g., Homo sapiens) and in vitro/in vivo status.
* `Tissue` & `Cell_types`: Maps to UBERON and CL ontologies.
* `Dataset`: Tracks related GEO identifiers or clinical datasets.

---

# Installation

## Requirements

* Python 3.9+

## Python Dependencies

Install the required packages:

```bash
pip install requests pandas
```
---

# Usage

Run the script directly from your terminal to construct and populate the database:

```bash
python database_creation.py
```

*Note: If you encounter GitHub API rate-limiting restrictions during large-scale scans, paste a valid personal access token into the `GITHUB_TOKEN` variable inside the script.*

---

# Pipeline Steps

## 1. Schema Initialization
The script establishes a connection to SQLite, enforces `PRAGMA foreign_keys = ON`, and runs a comprehensive DDL script. This builds the foundational structure for models, pathogens, hosts, tissues, cells, and datasets alongside their corresponding junction tables.

## 2. Remote Repository Scanning
The pipeline queries the GitHub tree API for the `thalieh/InfectioGIT` repository. It filters the recursive file system list to isolate every `.json` file located within the designated target results folder.

## 3. Dynamic Source Detection
As each JSON file is pulled into memory, the script evaluates its root keys. Records containing specific biological fields are categorized as BioModels, while records with custom programming languages or DOI strings are routed to Zenodo or NCBI handlers.

## 4. Relational Ingestion & Deduplication
Extracted entities are inserted sequentially into the database. The system utilizes `INSERT OR IGNORE` logic across all primary and junction tables to ensure duplicate models, cross-references, or taxons do not corrupt database integrity.

## 5. Diagnostic Validation
Once the ingestion loop terminates, the pipeline executes a suite of predefined analytical SQL queries via Pandas. It outputs row counts for every table and displays sample joins combining models with their corresponding diseases, hosts, and pathogen classes.

---

# Database Schema

The script manages a fully normalized relational structure containing the following core tables:

* **Model**: Tracks the primary DOI, title, mathematical formalism, biological scale, publication year, software environment, and source database origin.
* **Disease**: Stores standardized DOID/Mondo ontology identifiers and matching common disease names.
* **Pathogen**: Maps NCBI Taxonomy IDs to species names, life cycle stages, and classifications (e.g., virus, bacterium).
* **Host**: Tracks taxonomic information and maps experimental environments to in vitro, in vivo, or in silico tags.
* **Tissue & Cell_types**: Stores structural anatomical context mapped to UBERON and CL ontologies.
* **Dataset**: Logs linked external accessions such as Gene Expression Omnibus (GEO) identifiers.

---

# Outputs

The execution pipeline generates a physical relational database asset alongside comprehensive console reporting:

| Asset | Description |
| ----------------------------------- | --------------------------------- |
| `infectio_git.db` | The populated SQLite database containing all relational tables and constraints. |
| Table Row Summaries | Direct console printout displaying the exact number of populated records per entity. |
| Pandas DataFrame Highlights | Analytical matrix logs demonstrating valid relational joins across pathogens and models. |

---

# Authors

Developed for large-scale infectious disease database consolidation, cross-repository standardization, and relational model discovery in computational biology and bioinformatics.
