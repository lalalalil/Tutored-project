# BioModels Data Processing Scripts

This repository provides two Python scripts designed to **retrieve, organize, clean, and analyze computational models** from the BioModels database, with a focus on infectious diseases.

These scripts aim to:
- Structure downloaded data in a consistent way
- Separate **model files** from **non-model (metadata) files**
- Filter out irrelevant or unusable models
- Generate **statistics on file types and model content**

---

## About BioModels

BioModels is a curated repository of **quantitative biological models**, commonly used in systems biology and epidemiology. Models are typically stored in formats such as:
- SBML
- SED-ML
- OMEX archives

---

## Script 1 — Data Retrieval & Organization

This script:
- Queries BioModels using disease-related keywords
- Downloads associated files for each model
- Organizes them into a structured directory
- Separates **model files** from **metadata files**

### Targeted diseases
- Dengue
- Chikungunya
- Lyme disease
- Mpox
- West Nile virus
- Influenza
- Tuberculosis
- HIV
- COVID-19

### Output structure

```
BioModels_Data_cleaned_redirection_metadata_json_dico_/
│
├── <disease>/
│ ├── <model_id>/
│ │ ├── model/
│ │ └── metadata/
│
└── summary_extensions.csv
```


### Key features

- Search via BioModels API (`bioservices`, `biomodels`)
- Automatic directory creation
- File classification:
  - **model/** → modeling-related files
  - **metadata/** → non-model files

Files are classified as metadata if:
- Their extension is in: *.png, .jpg, .jpeg, .pdf, .txt, .docx, .doc, .xlsx, .xls, .csv*

- OR their name contains: *metadata, .json, .rdf, .owl*


- Retrieval of:
- Web JSON metadata
- Internal file metadata

- Output:
- `summary_extensions.csv` → number of files per extension per model

---

## Script 2 — Data Cleaning & Filtering

This script processes the downloaded data to **retain only valid computational models**.

### Model validation

A model is kept only if its `model/` folder contains at least one file with a recognized modeling extension.

### Recognized modeling extensions

*.xml, .sbml, .omex, .sedml, .cps, .m, .ode, .py, .f, .java, .vcml, .zip, .tsv, .yaml, .graphml*


### Removal criteria

A model directory is deleted if:
- No modeling files are found
- Only metadata or auxiliary files are present

### Output

- `statistiques_modeles_nettoyes.csv`:
  - Number of files per extension
  - Per model and per disease category

---

## Usage

### 1. Install dependencies

```bash
pip install bioservices biomodels requests
```

### 2. Run the scripts

```
python biomodels_model_extraction.py
python sort_metadata_files_and_delete_models_whithout_modelling_files.py
```


---

## Outputs

| File | Description |
|------|------------|
| `summary_extensions.csv` | File distribution across all downloaded models |
| `statistiques_modeles_nettoyes.csv` | Statistics for cleaned and valid models |

---

## Notes

- API requests are rate-limited (delays are included)
- Some files may fail to download
- Scripts focus on robustness but may require adaptation depending on use case