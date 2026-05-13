# NCBI PMC Article Extraction Pipeline

Download open access articles from PubMed Central, extract SBML models and executable code, and organize files by disease.

## Features

- 🔍 **Search PubMed** for articles by disease keyword
- 📥 **Download OA packages** from PMC
- 📦 **Extract files** (ZIP, TAR.GZ) including nested archives
- 🧬 **Identify SBML models** and executable code (.py, .r, .m, etc.)
- 🔗 **Extract repository links** (GitHub, GitLab, DOI)
- 📁 **Organize files** into data/ and metadata/ folders
- 🧹 **Cleanup script** to keep only SBML and JSON files

## Installation

```bash
pip install -r requirements.txt
```

## Setup

Edit the email address in `ncbi_extraction.py`:

```python
EMAIL = "your.email@ncbi.nlm.nih.gov"  # Required by NCBI
```

## Usage

### Download Articles

```bash
python ncbi_extraction.py
```

This will:
- Search PubMed for all diseases in `DISEASE_QUERIES`
- Download open access packages
- Extract all files
- Organize by disease and article (PMCID)

Output: `~/Desktop/downloaded_ncbi_models/`

### Preview Cleanup (Dry Run)

```bash
python cleanup_files.py
```

Shows what will be deleted (does not delete anything).

### Actually Cleanup

```bash
python cleanup_files.py --delete
```

Keeps ONLY `.sbml` and `.json` files. Deletes everything else.

## Customization

### Change Diseases

Edit `DISEASE_QUERIES` in `ncbi_extraction.py`:

```python
DISEASE_QUERIES = {
    "dengue_files": '(dengue OR DENV)',
    "your_disease_files": '(your disease OR synonym)',
}
```

### Change File Types to Keep

Edit `KEEP_EXTENSIONS` in `cleanup_files.py`:

```python
KEEP_EXTENSIONS = {".sbml", ".json", ".py"}  # Add .py to keep Python files
```




## Dependencies

- `biopython` - NCBI Entrez queries
- `requests` - HTTP downloads

## License

MIT

## Notes

- NCBI requires a valid email address in Entrez queries
- Not all articles have downloadable OA packages (404s are expected)
- The script respects NCBI rate limits (0.34s between requests)
- Files are organized by disease and PMCID for easy traceability
