# NCBI Extraction (PMC Open Access) — Infectious Disease Model Retrieval

This folder contains a notebook/script used to retrieve **open-access supplementary files** from **NCBI PubMed Central (PMC)** for multiple infectious diseases, with a focus on **computational model assets** (SBML/OMEX/SED-ML) and **repository links** (GitHub/GitLab/Bitbucket/DOI) when available.

## What it does
For each disease query:
1. Search **PubMed** using the disease term + model/file keywords (SBML, OMEX, COMBINE archive, etc.) and a minimum publication year.
2. Convert **PMIDs → PMCIDs** (only PMC entries can provide OA packages).
3. For each PMCID:
   - Fetch PMC **Open Access package** links (ZIP / TAR.GZ / TGZ)
   - Download and **extract** archives (including nested archives)
   - Separate extracted files into:
     - `data/` → model/executable-like files (SBML/XML with `<sbml>`, OMEX, SED-ML, scripts, notebooks, etc.)
     - `metadata/` → PDFs, images, text docs, manifests, JSON, etc.
   - Extract and save repository links into `metadata/repository_links.csv`
4. Generate manifests and summaries.

## Diseases / Queries
Default disease queries are defined in the notebook/script (example keys):
- dengue
- chikungunya
- lyme
- mpox
- west_nile
- influenza
- tuberculosis
- hiv
- covid

## Output structure
Example output (one folder per disease):
downloaded_ncbi_models/
  chikungunya_files/
    pubmed_records.csv
    pubmed_records.json
    pmc_mapping.csv
    pmc_mapping.json
    disease_summary.csv
    disease_summary.json
    PMCxxxxxxx/
      data/
        data_manifest.csv
        data_manifest.json
      metadata/
        metadata_manifest.csv
        metadata_manifest.json
        repository_links.csv
        pubmed_linked_records.csv

## Requirements
- Python 3
- biopython
- requests

Optional (only if Excel conversion is enabled):
- pandas
- openpyxl

## NCBI Entrez email (required)
NCBI Entrez requests require an email address.
**Do not commit your personal email to GitHub.**  
Set it locally in the notebook before running:

```python
EMAIL = "your.email@domain.com"
run_all_diseases(email=EMAIL)
