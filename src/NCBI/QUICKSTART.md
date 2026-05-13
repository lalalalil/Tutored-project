# Quick Start

## 1. Install

```bash
pip install -r requirements.txt
```

## 2. Set Your Email

Edit `ncbi_extraction.py` and change:

```python
EMAIL = "sebbahaya03@gmail.com"
```

To your email address.

## 3. Download Articles

```bash
python ncbi_extraction.py
```

Wait 1-3 hours depending on number of articles.

Output goes to: `~/Desktop/downloaded_ncbi_models/`

## 4. Preview What Will Be Deleted

```bash
python cleanup_files.py
```

This shows what files will be removed (keeps only `.sbml` and `.json`).

## 5. Actually Delete Unnecessary Files

```bash
python cleanup_files.py --delete
```

This removes all files except SBML models and JSON metadata.

## Done!

Your extracted SBML models are in `~/Desktop/downloaded_ncbi_models/`

Each disease folder contains articles organized by PMCID (unique identifier).
