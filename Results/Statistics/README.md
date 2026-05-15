# Results Directory

## Overview

This directory contains the final datasets and statistical outputs generated during the large-scale collection and analysis of infectious disease computational models from multiple public repositories.

The results are organized by:

1. disease category
2. source database

The datasets originate from three major repositories:

* BioModels
* Zenodo
* National Center for Biotechnology Information

The directory also contains statistical summaries and visualization outputs generated from the collected datasets.

---

# Directory Structure

```text id="u0vh8l"
Results/
├── COVID-19/
│   ├── BioModels/
│   ├── Zenodo/
│   └── NCBI/
│
├── HIV/
│   ├── BioModels/
│   ├── Zenodo/
│   └── NCBI/
│
├── Tuberculosis/
│   ├── BioModels/
│   ├── Zenodo/
│   └── NCBI/
│
├── Influenza/
│   ├── BioModels/
│   ├── Zenodo/
│   └── NCBI/
│
├── Dengue/
│   ├── BioModels/
│   ├── Zenodo/
│   └── NCBI/
│
├── Chikungunya/
│   ├── BioModels/
│   ├── Zenodo/
│   └── NCBI/
│
├── Mpox/
│   ├── BioModels/
│   ├── Zenodo/
│   └── NCBI/
│
├── West_Nile/
│   ├── BioModels/
│   ├── Zenodo/
│   └── NCBI/
│
├── Statistics/
└── Outputs/
```

---

# Disease Directories

Each disease directory contains computational models and associated metadata collected from the three repositories.

## Included Repositories

### BioModels

Contains:

* curated computational models
* mathematical model files
* standardized metadata
* publication information
* modelling approaches

---

### Zenodo

Contains:

* archived datasets
* supplementary modelling files
* repository metadata
* deposited research resources

---

### NCBI

Contains:

* biological datasets
* supplementary computational resources
* metadata records
* associated repository files

---

# Statistics Directory

```text id="j26l3y"
Results/Statistics/
```

Contains statistical summary tables generated during the analysis workflow.

### Example Files

```text id="x3oq6d"
publication_dates_summary.csv
modelling_approaches_summary.csv
statistiques_modeles_nettoyes.csv
summary_avg_per_category.csv
```

### Included Analyses

* publication statistics
* modelling approach distributions
* repository composition summaries
* file extension analyses
* curation quality statistics

---

# Outputs Directory

```text id="wq9b7m"
Results/Outputs/
```

Contains generated visualizations and final analysis figures.

### Included Figures

* publication timelines
* disease comparison plots
* modelling approach distributions
* file extension analyses
* heatmaps
* curation quality visualizations

### File Formats

* `.png`
* `.csv`

---

# Objectives

The datasets and analyses in this directory were generated to support:

* infectious disease modelling studies
* repository comparison
* FAIR model integration
* metadata harmonization
* computational biology research
* bioinformatics analyses
* reproducible workflows

---

# Workflow Summary

```text id="s7w8ew"
Public Repositories
(BioModels / Zenodo / NCBI)
            ↓
Disease-Specific Collection
            ↓
Metadata Cleaning & Standardization
            ↓
Model Classification & Filtering
            ↓
Statistical Analysis
            ↓
Final Results & Visualizations
```

---

# Notes

* Repository structures may vary depending on the source database.
* Some models include supplementary metadata and publication files.
* Statistical analyses and figures were generated automatically using dedicated Python and R workflows.
* Outputs are formatted for downstream analyses, reports, and publications.

---

# Authors

Generated as part of a large-scale infectious disease computational model integration and comparative analysis project.