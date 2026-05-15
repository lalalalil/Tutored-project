# BioModels Analysis Pipeline

## Overview

This project provides a complete Python pipeline for downloading, organizing, cleaning, filtering, and analyzing computational models from BioModels.

The pipeline automatically:

* Downloads disease-related models from BioModels
* Separates modeling files from metadata files
* Consolidates metadata into structured JSON files
* Retrieves publication DOIs
* Cleans invalid or incomplete models
* Organizes models by curation status and modeling approach
* Filters models by publication year
* Generates CSV statistics files
* Launches an R statistical analysis workflow

The script is designed for large-scale comparative analysis of infectious disease computational models.

---

# Features

* Automated download of BioModels entries
* DOI extraction from publication metadata
* Separation of:

  * model files
  * metadata files
* Consolidated metadata generation (`*_metadata.json`)
* Removal of incomplete/non-model repositories
* Classification by:

  * curation status
  * modelling approach
* Publication year filtering
* CSV export for downstream analysis
* Integration with an R visualization/statistics script

---

# Supported Diseases

The pipeline currently searches models related to:

* Dengue
* Chikungunya
* Mpox
* West Nile Virus
* Influenza
* Tuberculosis
* HIV
* COVID-19

---

# Project Structure

```text
.
├── biomodels_pipeline.py
├── stats_biomodels.R
├── BioModels_Database_Final/
│   ├── dengue_files/
│   ├── covid_files/
│   ├── hiv_files/
│   └── ...
│
├── BioModels_Stats/
│   ├── statistiques_modeles_nettoyes.csv
│   ├── modelling_approaches_summary.csv
│   ├── publication_dates_summary.csv
│   └── outputs/
│       ├── figures/
│       └── plots/
```

---

# Installation

## Requirements

* Python 3.9+
* R (optional, for statistical analysis)

## Python Dependencies

Install required packages:

```bash
pip install bioservices biomodels requests
```

---

# Usage

Run the script:

```bash
python biomodels_pipeline.py
```

The script launches an interactive menu:

```text
--- BIOMODELS ANALYSIS PIPELINE ---

1. Download & Sort (Metadata vs Model)
2. Clean Empty Models & Generate Extension Stats CSV
3. Separate by Curation Status
4. Classify Folders by Modelling Approach
5. Delete Models Published Before 2015
6. Generate Publication Dates CSV
7. Clean Metadata Directories
8. Run R Statistical Analysis
9. Exit
```

---

# Pipeline Steps

## 1. Download & Sort

Downloads models from BioModels and:

* creates disease-specific directories
* separates:

  * model files
  * metadata files
* generates a consolidated metadata JSON file
* retrieves DOI information from publication records

Generated metadata file:

```text
MODEL_ID_metadata.json
```

---

## 2. Clean Empty Models & Generate Statistics

This step:

* removes repositories without valid modeling files
* counts file extensions
* generates statistics CSV files

Output:

```text
BioModels_Stats/statistiques_modeles_nettoyes.csv
```

---

## 3. Separate by Curation Status

Models are reorganized into folders according to:

* curated
* non_curated
* unknown

Example:

```text
covid_files/
├── curated/
├── non_curated/
```

---

## 4. Classify by Modelling Approach

Models are reorganized using their BioModels modelling approach metadata.

Examples:

* Ordinary Differential Equation
* Flux Balance Analysis
* Stochastic Modeling

Output:

```text
BioModels_Stats/modelling_approaches_summary.csv
```

---

## 5. Delete Models Published Before 2015

Creates comparative CSV summaries before filtering and optionally removes:

* models published before 2015

Generated files:

```text
stats_BEFORE_2015_filter.csv
stats_AFTER_2015_filter.csv
```

---

## 6. Generate Publication Dates CSV

Exports publication metadata into:

```text
publication_dates_summary.csv
```

Including:

* model ID
* disease category
* publication year
* journal
* DOI

---

## 7. Clean Metadata Directories

Deletes all files in metadata directories except:

```text
*_metadata.json
```

Useful for reducing storage and standardizing repositories.

---

## 8. Run R Statistical Analysis

Runs:

```text
stats_biomodels.R
```

The R script uses generated CSV files to produce:

* plots
* summary figures
* statistical analyses

Outputs are saved to:

```text
BioModels_Stats/outputs/
```

---

# Metadata Structure

Example of generated metadata JSON:

```json
{
    "model_id": "MODEL1234567890",
    "doi": "10.xxxx/xxxxx",
    "web_metadata": {},
    "files_list": []
}
```

---

# Recognized File Types

## Model Files

Examples:

* `.xml`
* `.sbml`
* `.omex`
* `.sedml`
* `.py`
* `.m`
* `.ode`
* `.graphml`

## Metadata Files

Examples:

* `.json`
* `.rdf`
* `.owl`
* `.pdf`
* `.csv`
* `.xlsx`

---

# Outputs

The pipeline generates:

| File                                | Description                       |
| ----------------------------------- | --------------------------------- |
| `statistiques_modeles_nettoyes.csv` | Extension statistics              |
| `modelling_approaches_summary.csv`  | Modelling approach classification |
| `publication_dates_summary.csv`     | Publication metadata              |
| `stats_BEFORE_2015_filter.csv`      | Models before filtering           |
| `stats_AFTER_2015_filter.csv`       | Models after filtering            |

---

# Notes

* DOI extraction uses both:

  * BioModels metadata
  * EBI Publications API fallback
* The script includes safety confirmations before destructive operations.
* Invalid or incomplete repositories are automatically removed.

---

# Example Workflow

Recommended execution order:

```text
1 → 2 → 3 → 4 → 5 → 6 → 8
```

---

# Authors

Developed for large-scale infectious disease model integration and comparative analysis workflows in computational biology and bioinformatics.

---

# License

This project is distributed under the MIT License.