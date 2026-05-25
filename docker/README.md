# InfectioGIT Docker Setup

## Overview
Docker configuration for the InfectioGIT multi-source infectious disease data analysis platform.

### Supported Data Sources
- **BioModels**: Download and analyze computational models from BioModels database
- **NCBI**: Extract open-access articles and supplementary data from PubMed Central
- **Zenodo**: Retrieve research datasets and models from Zenodo
- **SQLite Database**: Unified database schema for disease, pathogen, host, tissue, and cell type relationships

### Supported Diseases
- COVID-19 (SARS-CoV-2)
- Influenza
- Dengue
- Chikungunya
- West Nile Virus
- Tuberculosis
- HIV
- Mpox
- Lyme Disease (in NCBI module)

---

## Quick Start

### Prerequisites
- [Docker Desktop](https://www.docker.com/products/docker-desktop)

### Build the Docker image
```bash
docker-compose build
```

### Run interactive Python shell
```bash
docker-compose run infectio python -i
```

### Run specific scripts

#### BioModels Pipeline
Interactive menu-driven pipeline for downloading, sorting, and analyzing BioModels:
```bash
docker-compose run infectio python src/Biomodels/biomodels_pipeline.py
```

#### NCBI Extraction
Download open-access articles and code from PubMed Central:
```bash
docker-compose run infectio python src/NCBI/ncbi_extraction.py
```

#### Database Creation
Scan GitHub repository and populate SQLite database:
```bash
docker-compose run infectio python src/Database/database_creation.py
```

---

## Project Structure

```
InfectioGIT/
├── docker/                          # Docker configuration
│   ├── Dockerfile                   # Multi-stage build for all dependencies
│   ├── docker-compose.yml           # Service orchestration
│   ├── requirements.txt             # Python dependencies
│   └── README.md                    # This file
├── src/                             # Source code
│   ├── Biomodels/
│   │   └── biomodels_pipeline.py   # BioModels analysis (8-step pipeline)
│   ├── NCBI/
│   │   ├── ncbi_extraction.py      # PubMed Central extraction
│   │   └── cleanup_files.py        # File cleanup utility
│   ├── Zenodo/
│   │   └── zenodo_extraction.py    # Zenodo data retrieval
│   └── Database/
│       └── database_creation.py    # SQLite database schema & population
├── Results/                         # Downloaded models and results
├── BioModels_Database_Final/        # BioModels data (created at runtime)
├── BioModels_Stats/                 # BioModels statistics outputs
└── downloaded_ncbi_models/          # NCBI extracted files (created at runtime)
```

---

## Key Dependencies

### Python Libraries
- **biopython** (1.81): NCBI Entrez queries, sequence analysis
- **bioservices** (1.11.2): BioModels API access
- **pandas** (2.0.3): Data manipulation and CSV output
- **requests** (2.31.0): HTTP operations for API calls
- **numpy** (1.24.3): Numerical computing

### System Dependencies
- **R 4.x**: Statistical analysis (used by BioModels pipeline)
- **curl, wget**: File downloads
- **build-essential**: C/C++ compilation for biopython
- **git**: Version control (already in project)

---

## Usage Examples

### Example 1: Download BioModels and Generate Statistics

```bash
docker-compose run infectio python src/Biomodels/biomodels_pipeline.py

# Interactive menu:
# 1. Download & Sort (metadata vs model)
# 2. Clean Empty Models & Extension Stats CSV
# 3. Separate by Curation Status
# 4. Classify by Modelling Approach
# 5. Generate Publication Dates CSV
# 6. Delete Models Before 2015
# 7. Clean Metadata Dirs
# 8. Run R Statistical Analysis
# 9. Exit
```

### Example 2: Extract NCBI Articles with Code

```bash
docker-compose run infectio python src/NCBI/ncbi_extraction.py

# Downloads open-access articles from PubMed Central for:
# - COVID-19
# - Influenza
# - Dengue
# - Tuberculosis
# - HIV
# etc.
```

### Example 3: Populate Database from GitHub

```bash
docker-compose run infectio python src/Database/database_creation.py

# Scans Results/ folder in GitHub repo
# Creates SQLite database with:
# - Model metadata (DOI, title, formalism, year)
# - Disease associations (DOID/MONDO IDs)
# - Pathogen/Host taxonomy relationships
# - Tissue/Cell type ontologies
# - Dataset links
```

---

## Volume Mounts

The docker-compose configuration provides these mounts:

| Mount Point | Purpose |
|-------------|---------|
| `.:/app` | Project root (source code, data folders) |
| `~/:` | User home directory (for Desktop access) |
| `~/Desktop/infectio_outputs:/app/outputs` | Output files sent to Desktop |

---

## Output Locations

### BioModels Pipeline Outputs
- `./BioModels_Database_Final/` - Downloaded models organized by disease/status/approach
- `./BioModels_Stats/` - CSV statistics and R-generated PNG plots
- Files: `publication_dates_summary.csv`, `modelling_approaches_summary.csv`, etc.

### NCBI Extraction Outputs
- `~/Desktop/downloaded_ncbi_models/{disease}_files/{PMCID}/`
  - `data/` - Executable code and models
  - `metadata/` - Documentation and images
  - `*_manifest.csv` - File inventory
  - `*_manifest.json` - Detailed metadata

### Database Outputs
- `./infectio_git.db` - SQLite database with full schema

---

## Environment Variables

```bash
PYTHONUNBUFFERED=1          # Unbuffered Python output
PYTHONDONTWRITEBYTECODE=1   # No .pyc file generation
ENTREZ_EMAIL=...            # Email for NCBI Entrez API (optional)
```

---

## Troubleshooting

### NCBI Rate Limiting
If you encounter PubMed Central rate limits:
- Add delays between requests (already implemented)
- Use smaller `MAX_PMIDS_PER_DISEASE` in ncbi_extraction.py

### R Script Issues
If R statistical analysis fails:
- Ensure R is installed in the container (included in Dockerfile)
- Place `stats_biomodels.R` in the project root
- Check R script dependencies (ggplot2, etc.)

### Memory Issues
For large downloads:
- Run Docker with increased memory: `docker-compose run --memory=8g infectio ...`
- Process by disease category instead of all at once

---

## Development

### Install Additional Packages
```bash
docker-compose run infectio pip install <package-name>
```

### Debug Mode
```bash
docker-compose run infectio python -i
# Then import and test:
# >>> import src.NCBI.ncbi_extraction as ncbi
# >>> ncbi.build_pubmed_query("covid")
```

### View Logs
```bash
docker-compose logs -f infectio
```

---

## References

- BioModels: https://www.ebi.ac.uk/biomodels/
- NCBI PubMed Central: https://www.ncbi.nlm.nih.gov/pmc/
- Zenodo: https://zenodo.org/
- DisGeNET: https://www.disgenet.org/ (disease/gene associations)
- DOID/MONDO: Disease ontologies

---

## Contact
For issues or questions: sebbahaya03@gmail.com
