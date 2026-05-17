# Docker Setup

This project includes a complete Docker configuration for running the NCBI/BioModels/Zenodo extraction pipeline in a containerized environment.

## Quick Start

### Prerequisites
- [Docker Desktop](https://www.docker.com/products/docker-desktop)

### Build & Run

```bash
# Build the Docker image
docker-compose build

# Start interactive Python shell
docker-compose up

# Run NCBI extraction
docker-compose run infectio python ncbi_extraction.py

# Run Zenodo extraction
docker-compose run infectio python InfectioGIT/src/Zenodo/zenodo_extraction.py

# Run BioModels extraction
docker-compose run infectio python InfectioGIT/src/Biomodels/biomodels_pipeline.py

# Run cleanup
docker-compose run infectio python cleanup_files.py --delete
```

## Files

- **Dockerfile** — Container image definition with Python 3.11 + system dependencies
- **docker-compose.yml** — Container orchestration and volume mounting
- **requirements.txt** — Python package dependencies

## Environment

- Python 3.11 (slim image)
- System tools: build-essential, curl, wget, git
- Python packages: biopython, requests, bioservices, biomodels, numpy

## Output

Container outputs go to:
- `/app/outputs` (mounted to `~/Desktop` on your machine)
- Docker volume mounts: `.:/app` (current directory)

## Included Modules

-  **NCBI extraction** — `ncbi_extraction.py`
-  **Zenodo extraction** — `InfectioGIT/src/Zenodo/zenodo_extraction.py`
-  **BioModels extraction** — `InfectioGIT/src/Biomodels/biomodels_pipeline.py`
-  **File cleanup** — `cleanup_files.py`

## Notes

- Container runs with `PYTHONUNBUFFERED=1` for real-time logging
- All downloads save to `~/Desktop/downloaded_ncbi_models/`
- Remove `downloaded_ncbi_models 2/` folder before building to save space
- All source code in `InfectioGIT/` is included and ready to use
