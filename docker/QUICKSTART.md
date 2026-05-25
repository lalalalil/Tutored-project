# Docker Quick Start Guide

## First Time Setup

```bash
# Navigate to docker folder
cd docker

# Build the Docker image (takes ~3 minutes)
docker-compose build
```

## Verify Installation

```bash
docker-compose run infectio python -c "
from Bio import Entrez
import bioservices
import pandas as pd
print('✓ All dependencies installed successfully!')
"
```

---

## Running the Pipelines

### Option 1: BioModels Analysis (Interactive Menu)
```bash
docker-compose run infectio python src/Biomodels/biomodels_pipeline.py
```
Select menu option 1-8 to run different pipeline steps.

### Option 2: NCBI PubMed Central Extraction
```bash
docker-compose run infectio python src/NCBI/ncbi_extraction.py
```
Downloads open-access articles for specified diseases.

### Option 3: Database Creation
```bash
docker-compose run infectio python src/Database/database_creation.py
```
Creates SQLite database from GitHub repository data.

---

## Interactive Python Shell

```bash
docker-compose run infectio python -i
```

Then use Python normally:
```python
from Bio import Entrez
import bioservices
import pandas as pd
# ... your code
exit()  # to exit
```

---

## Common Tasks

### Run Python Script
```bash
docker-compose run infectio python your_script.py
```

### Install Additional Package
```bash
docker-compose run infectio pip install package-name
```

### View Installed Packages
```bash
docker-compose run infectio pip list
```

### Check Python Version
```bash
docker-compose run infectio python --version
```

### Clean Up Orphan Containers
```bash
docker-compose down --remove-orphans
```

### View Container Logs
```bash
docker-compose logs -f infectio
```

### Increase Memory for Large Datasets
```bash
docker-compose run --memory=8g infectio python src/Biomodels/biomodels_pipeline.py
```

---

## Troubleshooting

**Issue:** "no configuration file provided: not found"
- **Solution:** Make sure you're in the `docker/` folder when running docker-compose commands

**Issue:** "ModuleNotFoundError: No module named 'biopython'"
- **Solution:** Use `from Bio import Entrez` instead of `import biopython`

**Issue:** Docker image takes long to build
- **Solution:** This is normal for first build (~3 minutes). Subsequent runs will be faster due to caching.

**Issue:** NCBI rate limiting errors
- **Solution:** The script includes delays. Reduce `MAX_PMIDS_PER_DISEASE` in `src/NCBI/ncbi_extraction.py` or run disease categories separately.

---

## Files in This Folder

- **Dockerfile** - Container image definition
- **docker-compose.yml** - Service configuration and volume mounts
- **requirements.txt** - Python package dependencies
- **README.md** - Comprehensive documentation
- **QUICKSTART.md** - This file (quick command reference)

---

## Environment Variables

Set in `docker-compose.yml`:

```yaml
PYTHONUNBUFFERED=1          # Real-time output
PYTHONDONTWRITEBYTECODE=1   # No .pyc files
ENTREZ_EMAIL=your@email.com # NCBI API email (optional)
```

---

## Volume Mounts

- `..:/app` - Project root directory
- `~/:/home/user` - User home directory
- `~/Desktop/infectio_outputs:/app/outputs` - Output directory

---

## Key Dependencies

- **Python 3.11** - Programming language
- **R 4.x** - Statistical analysis
- **biopython** - NCBI/sequence analysis
- **bioservices** - BioModels API
- **pandas** - Data manipulation
- **requests** - HTTP operations

See `requirements.txt` for full list with versions.

---

## Tips

✅ Always run from the `docker/` folder  
✅ Use `--remove-orphans` occasionally to clean up  
✅ Check `docker-compose.yml` to customize environment variables  
✅ Read `README.md` for detailed documentation  
✅ Large downloads may require significant disk space  

---

## Getting Help

1. Check `README.md` for detailed documentation
2. Review pipeline source code in `src/`
3. Check container logs: `docker-compose logs -f infectio`
4. Test with: `docker-compose run infectio python -c "import biopython; from Bio import Entrez; print('OK')"`

---

**Last Updated:** May 25, 2026  
**Docker Compose Version:** 3.x+  
**Docker Desktop Required:** Yes
