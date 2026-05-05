#!/usr/bin/env python3
"""
NCBI PMC Open Access Article Extraction Pipeline
=================================================

Downloads open access articles from PubMed Central for specified diseases,
extracts all files, identifies SBML models and executable code,
and organizes files for further analysis.

Usage:
    python ncbi_extraction.py

Output:
    ~/Desktop/downloaded_ncbi_models/
        {disease}_files/
            {PMCID}/
                data/              # Executable code + models
                metadata/          # Documentation + images
"""

from __future__ import annotations

import csv, json, logging, re, shutil, tarfile, time, zipfile
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import requests
from Bio import Entrez

# =============================================================================
# CONFIGURATION
# =============================================================================

OA_URL = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi"
REQUEST_TIMEOUT = 60
CHUNK_SIZE = 1024 * 256
ENTREZ_SUMMARY_BATCH = 100

logging.basicConfig(level=logging.INFO, format="%(levelname)s | %(message)s")
logger = logging.getLogger("ncbi_pipeline")

MIN_PUBLICATION_YEAR = 2015
MAX_PMIDS_PER_DISEASE: Optional[int] = None

# File type definitions
EXECUTABLE_EXTENSIONS = {
    ".py", ".r", ".m", ".ipynb",
    ".sbml", ".sedml", ".omex", ".cps", ".cellml",
    ".jl", ".cpp", ".c"
}

SUPPORT_DATA_EXTENSIONS = {
    ".csv", ".tsv", ".tab", ".txt", ".xls", ".xlsx", ".xlsm", ".ods",
    ".mat", ".h5", ".hdf5", ".sqlite", ".db",
    ".fasta", ".fa", ".fna", ".fastq", ".fq",
    ".gb", ".gbk", ".gff", ".gff3", ".bed", ".vcf",
    ".json", ".yaml", ".yml", ".graphml", ".gml", ".sif", ".xgmml",
    ".nwk", ".newick"
}

METADATA_EXTENSIONS = {
    ".pdf", ".png", ".jpg", ".jpeg", ".svg", ".gif", ".tif", ".tiff",
    ".md", ".rst", ".bib", ".ris", ".nfo", ".doc", ".docx", ".ppt", ".pptx"
}

METADATA_NAME_HINTS = {
    "readme", "license", "licence", "authors", "author", "manifest",
    "metadata", "supplement", "description", "legend", "notes", "protocol"
}

URL_REGEX = re.compile(r'https?://[^\s"<>\]\)]+', flags=re.IGNORECASE)

# Disease queries - EDIT THESE TO CUSTOMIZE
DISEASE_QUERIES = {
    "dengue_files": '(dengue OR DENV)',
    "chikungunya_files": '(chikungunya OR CHIKV)',
    "lyme_files": '(lyme OR borrelia OR borreliosis)',
    "mpox_files": '(mpox OR monkeypox)',
    "west_nile_files": '("west nile" OR WNV)',
    "influenza_files": '(influenza OR "influenza virus" OR "avian influenza" OR H5N1)',
    "tuberculosis_files": '(tuberculosis OR TB OR mycobacterium)',
    "hiv_files": '(HIV OR "human immunodeficiency virus")',
    "covid_files": '("covid" OR "SARS-CoV-2")',
}

DEFAULT_KEYWORDS = [
    '"supplementary material"', '"supplementary data"', '"source code"',
    "python", '"R script"', "matlab", "notebook",
    "github", "gitlab", "bitbucket",
    "SBML", "OMEX", '"COMBINE archive"',
]

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def ensure_dir(path: Path) -> Path:
    """Create directory if missing."""
    path.mkdir(parents=True, exist_ok=True)
    return path


def get_desktop_output_root() -> Path:
    """Get output root (Desktop or Bureau)."""
    home = Path.home()
    for candidate in [home / "Desktop", home / "Bureau"]:
        if candidate.exists():
            return ensure_dir(candidate / "downloaded_ncbi_models")
    return ensure_dir(home / "downloaded_ncbi_models")


def normalize_oa_url(url: str) -> str:
    """Convert FTP to HTTPS."""
    return url.replace("ftp://", "https://", 1) if url.startswith("ftp://") else url


def safe_read_text(path: Path, max_bytes: int = 200_000) -> str:
    """Read file as text, ignoring encoding errors."""
    try:
        return path.read_bytes()[:max_bytes].decode("utf-8", errors="ignore")
    except Exception:
        return ""


def is_sbml_xml(path: Path) -> bool:
    """Check if XML contains SBML."""
    try:
        return b"<sbml" in path.read_bytes()[:4000].lower()
    except Exception:
        return False


def is_archive_file(path: Path) -> bool:
    """Check if file is ZIP or TAR.GZ."""
    name = path.name.lower()
    return name.endswith(".zip") or name.endswith(".tar.gz") or name.endswith(".tgz")


def classify_file(path: Path) -> str:
    """Classify file: 'data', 'metadata', or 'ignore'."""
    suffix = path.suffix.lower()
    lower_name = path.name.lower()

    if suffix in EXECUTABLE_EXTENSIONS:
        return "data"
    if suffix == ".xml" and is_sbml_xml(path):
        return "data"
    if suffix in SUPPORT_DATA_EXTENSIONS:
        return "data"
    if suffix in METADATA_EXTENSIONS or any(h in lower_name for h in METADATA_NAME_HINTS):
        return "metadata"
    return "ignore"


def file_category_label(path: Path) -> str:
    """Get human-readable category."""
    suffix = path.suffix.lower()
    mapping = {
        ".py": "python", ".r": "r", ".m": "matlab", ".ipynb": "notebook",
        ".sbml": "sbml", ".sedml": "sedml", ".omex": "omex", ".cellml": "cellml",
        ".csv": "csv", ".json": "json", ".xml": "xml",
    }
    return mapping.get(suffix, suffix[1:] if suffix else "other")


def has_minimum_useful_content(root_dir: Path, extra_links: Optional[List[str]] = None) -> bool:
    """Keep if has executable/model, repo link, or 2+ data files."""
    strong_found = False
    support_count = 0

    for path in root_dir.rglob("*"):
        if not path.is_file():
            continue
        suffix = path.suffix.lower()

        if suffix in EXECUTABLE_EXTENSIONS or (suffix == ".xml" and is_sbml_xml(path)):
            strong_found = True
            break
        if suffix in SUPPORT_DATA_EXTENSIONS:
            support_count += 1

    if strong_found:
        return True

    if extra_links:
        repo_patterns = ["github.com", "gitlab.com", "bitbucket.org", "doi.org"]
        if any(any(p in x.lower() for p in repo_patterns) for x in extra_links):
            return True

    return support_count >= 2


# =============================================================================
# FILE I/O
# =============================================================================

def save_csv(rows: List[Dict], path: Path) -> None:
    """Save list of dicts as CSV."""
    if not rows:
        return

    all_keys = set()
    for row in rows:
        all_keys.update(row.keys())

    preferred = ["pmcid", "pmid", "filename", "category"]
    fieldnames = [k for k in preferred if k in all_keys] + sorted(all_keys - set(preferred))

    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def save_json(data: dict, path: Path) -> None:
    """Save dict as pretty JSON."""
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)


def summarize_by_category(rows: List[Dict]) -> Dict[str, int]:
    """Count files per category."""
    counts: Dict[str, int] = {}
    for r in rows:
        cat = r.get("category", "other")
        counts[cat] = counts.get(cat, 0) + 1
    return counts


# =============================================================================
# REPOSITORY LINK EXTRACTION
# =============================================================================

def extract_repo_links_from_text(text: str) -> List[str]:
    """Extract GitHub/GitLab/Bitbucket/DOI links."""
    urls = URL_REGEX.findall(text or "")
    kept = []
    for u in urls:
        low = u.lower()
        if any(p in low for p in ["github", "gitlab", "bitbucket", "doi.org"]):
            kept.append(u.rstrip(".,);]}>"))
    return sorted(set(kept))


def extract_repo_links_from_pubmed_records(records: List[Dict]) -> List[str]:
    """Extract repo links from PubMed metadata."""
    links = set()
    for rec in records:
        for value in rec.values():
            if isinstance(value, str):
                links.update(extract_repo_links_from_text(value))
    return sorted(links)


def extract_repo_links_from_directory(root: Path) -> List[str]:
    """Scan extracted files for repo links."""
    links = set()
    text_extensions = {".txt", ".md", ".rst", ".json", ".yaml", ".yml", ".py", ".r", ".m", ".xml", ".html"}

    for path in root.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in text_extensions:
            continue
        links.update(extract_repo_links_from_text(safe_read_text(path)))

    return sorted(links)


# =============================================================================
# ARCHIVE EXTRACTION
# =============================================================================

def recursively_extract_nested_archives(root_dir: Path, max_passes: int = 3) -> None:
    """Extract any nested ZIP/TAR files."""
    seen: Set[Path] = set()

    for _ in range(max_passes):
        found_new = False

        for path in root_dir.rglob("*"):
            if not path.is_file() or path in seen or not is_archive_file(path):
                continue

            lower = path.name.lower()
            target = path.parent / f"{path.stem}_extracted"
            if lower.endswith(".tar.gz"):
                target = path.parent / f"{path.name[:-7]}_extracted"
            elif lower.endswith(".tgz"):
                target = path.parent / f"{path.name[:-4]}_extracted"

            try:
                if lower.endswith(".zip"):
                    with zipfile.ZipFile(path, "r") as z:
                        z.extractall(target)
                else:
                    with tarfile.open(path, "r:gz") as t:
                        t.extractall(target, filter="data")
                seen.add(path)
                found_new = True
            except Exception as e:
                logger.warning("Nested extraction failed for %s: %s", path.name, e)

        if not found_new:
            break


def collect_and_copy_files(
    source_root: Path,
    pmcid: str,
    data_dir: Path,
    metadata_dir: Path
) -> Tuple[List[Dict], List[Dict]]:
    """Copy files to organized folders."""
    d_rows, m_rows = [], []
    copied_data, copied_meta = set(), set()

    for path in source_root.rglob("*"):
        if not path.is_file() or is_archive_file(path):
            continue

        kind = classify_file(path)
        if kind == "ignore":
            continue

        rel_path = path.relative_to(source_root)
        row = {
            "pmcid": pmcid,
            "filename": path.name,
            "relative_path": str(rel_path),
            "category": file_category_label(path),
        }

        if kind == "data":
            dest = data_dir / rel_path
            if dest not in copied_data:
                ensure_dir(dest.parent)
                shutil.copy2(path, dest)
                copied_data.add(dest)
            d_rows.append(row)
        else:
            dest = metadata_dir / rel_path
            if dest not in copied_meta:
                ensure_dir(dest.parent)
                shutil.copy2(path, dest)
                copied_meta.add(dest)
            m_rows.append(row)

    return d_rows, m_rows


# =============================================================================
# NCBI ENTREZ QUERIES
# =============================================================================

def build_pubmed_query(disease_term: str, keywords: Optional[List[str]] = None) -> str:
    """Build PubMed query with date filter."""
    if keywords is None:
        keywords = DEFAULT_KEYWORDS
    date_filter = f'("{MIN_PUBLICATION_YEAR}/01/01"[Date - Publication] : "3000"[Date - Publication])'
    return f"({disease_term}) AND ({' OR '.join(keywords)}) AND {date_filter}"


def search_all_pubmed_ids(query: str, batch_size: int = 200) -> List[str]:
    """Search PubMed and return all PMIDs."""
    pmids: List[str] = []
    retstart = 0
    total_count = None

    while True:
        current_retmax = batch_size
        if MAX_PMIDS_PER_DISEASE:
            remaining = MAX_PMIDS_PER_DISEASE - len(pmids)
            if remaining <= 0:
                break
            current_retmax = min(current_retmax, remaining)

        with Entrez.esearch(db="pubmed", term=query, retstart=retstart, retmax=current_retmax) as handle:
            result = Entrez.read(handle)

        if total_count is None:
            total_count = int(result.get("Count", 0))
            logger.info("PubMed total: %s", total_count)

        batch = result.get("IdList", [])
        if not batch or len(pmids) >= total_count:
            break

        pmids.extend(batch)
        retstart += len(batch)
        logger.info("PMIDs: %s / %s", len(pmids), total_count)
        time.sleep(0.34)

    return pmids


def fetch_pubmed_summaries(pmids: List[str]) -> List[Dict]:
    """Fetch metadata for PMIDs."""
    summaries: List[Dict] = []

    for i in range(0, len(pmids), ENTREZ_SUMMARY_BATCH):
        chunk = pmids[i:i + ENTREZ_SUMMARY_BATCH]

        with Entrez.esummary(db="pubmed", id=",".join(chunk), retmode="xml") as handle:
            result = Entrez.read(handle)

        docs = result if isinstance(result, list) else result.get("DocumentSummarySet", {}).get("DocumentSummary", [])
        for doc in docs:
            try:
                summaries.append(dict(doc))
            except Exception:
                summaries.append({"uid": str(doc)})

        logger.info("Summaries: %s / %s", min(i + ENTREZ_SUMMARY_BATCH, len(pmids)), len(pmids))
        time.sleep(0.34)

    return summaries


def pubmed_to_pmc(pmids: List[str]) -> Tuple[List[str], Dict[str, str]]:
    """Map PubMed IDs to PMC IDs."""
    if not pmids:
        return [], {}

    pmid_to_pmcid: Dict[str, str] = {}
    pmcids = set()

    for i in range(0, len(pmids), ENTREZ_SUMMARY_BATCH):
        chunk = pmids[i:i + ENTREZ_SUMMARY_BATCH]

        with Entrez.elink(dbfrom="pubmed", db="pmc", id=",".join(chunk), linkname="pubmed_pmc") as handle:
            results = Entrez.read(handle)

        for linkset in results:
            pmid_list = linkset.get("IdList", [])
            source_pmid = str(pmid_list[0]) if pmid_list else None

            for db in linkset.get("LinkSetDb", []):
                for item in db.get("Link", []):
                    pmcid = f"PMC{item['Id']}"
                    pmcids.add(pmcid)
                    if source_pmid:
                        pmid_to_pmcid[source_pmid] = pmcid

        logger.info("PMID→PMCID: %s / %s", min(i + ENTREZ_SUMMARY_BATCH, len(pmids)), len(pmids))
        time.sleep(0.34)

    return sorted(pmcids), pmid_to_pmcid


def get_pmc_oa_links(pmcid: str) -> List[str]:
    """Get OA download links for an article."""
    try:
        resp = requests.get(OA_URL, params={"id": pmcid}, timeout=REQUEST_TIMEOUT)
        resp.raise_for_status()
        hrefs = re.findall(r'href="([^"]+)"', resp.text)
        return [normalize_oa_url(h) for h in hrefs if is_archive_file(Path(h))]
    except Exception as e:
        logger.warning("OA links failed for %s: %s", pmcid, e)
        return []


# =============================================================================
# ARTICLE PROCESSING
# =============================================================================

def process_pmc_article(pmcid: str, disease_dir: Path, records: Dict, mapping: Dict) -> Dict:
    """Download, extract, and organize a single article."""
    pmc_dir = ensure_dir(disease_dir / pmcid)
    tmp_dir = ensure_dir(pmc_dir / "_tmp_raw")
    data_dir = ensure_dir(pmc_dir / "data")
    metadata_dir = ensure_dir(pmc_dir / "metadata")

    oa_links = get_pmc_oa_links(pmcid)

    linked_records = [records[pmid] for pmid in records if mapping.get(pmid) == pmcid]
    repo_links_from_records = extract_repo_links_from_pubmed_records(linked_records)

    report = {
        "pmcid": pmcid,
        "status": "no_archive_found",
        "kept": False,
        "repo_links": repo_links_from_records,
    }

    if not oa_links:
        shutil.rmtree(pmc_dir, ignore_errors=True)
        return report

    # Download and extract OA packages
    any_extracted = False
    for i, link in enumerate(oa_links, start=1):
        try:
            archive_path = tmp_dir / Path(link).name
            extract_dir = tmp_dir / f"archive_{i}"

            with requests.get(link, stream=True, timeout=REQUEST_TIMEOUT) as r:
                r.raise_for_status()
                with open(archive_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=CHUNK_SIZE):
                        if chunk:
                            f.write(chunk)

            lower = archive_path.name.lower()
            if lower.endswith(".zip"):
                with zipfile.ZipFile(archive_path, "r") as z:
                    z.extractall(extract_dir)
            else:
                with tarfile.open(archive_path, "r:gz") as t:
                    t.extractall(extract_dir, filter="data")

            any_extracted = True
        except Exception as e:
            logger.warning("Extraction failed for %s: %s", pmcid, e)

    if not any_extracted:
        shutil.rmtree(pmc_dir, ignore_errors=True)
        return {"pmcid": pmcid, "status": "extract_failed", "kept": False}

    # Extract nested archives
    recursively_extract_nested_archives(tmp_dir)

    # Find repo links in extracted files
    repo_links_from_files = extract_repo_links_from_directory(tmp_dir)
    all_repo_links = sorted(set(repo_links_from_records + repo_links_from_files))

    # Keep or discard based on content
    if not has_minimum_useful_content(tmp_dir, extra_links=all_repo_links):
        shutil.rmtree(pmc_dir, ignore_errors=True)
        return {"pmcid": pmcid, "status": "no_useful_content", "kept": False}

    # Copy files
    data_rows, metadata_rows = collect_and_copy_files(tmp_dir, pmcid, data_dir, metadata_dir)

    # Save manifests
    save_csv(data_rows, data_dir / "data_manifest.csv")
    save_json({
        "pmcid": pmcid,
        "repo_links": all_repo_links,
        "counts": summarize_by_category(data_rows),
        "files": data_rows
    }, data_dir / "data_manifest.json")

    if all_repo_links:
        repo_rows = [{"pmcid": pmcid, "url": u} for u in all_repo_links]
        save_csv(repo_rows, metadata_dir / "repository_links.csv")

    # Cleanup
    shutil.rmtree(tmp_dir, ignore_errors=True)

    return {
        "pmcid": pmcid,
        "status": "ok",
        "kept": True,
        "repo_links": all_repo_links,
        "data_counts": summarize_by_category(data_rows),
    }


# =============================================================================
# DISEASE PIPELINE
# =============================================================================

def fetch_ncbi_associated_data(disease_term: str, out_dir: str, email: str) -> Dict:
    """Run the complete pipeline for a disease."""
    Entrez.email = email
    base_dir = ensure_dir(Path(out_dir))

    query = build_pubmed_query(disease_term)
    logger.info("Query: %s", query)

    pmids = search_all_pubmed_ids(query)
    logger.info("PMIDs retrieved: %d", len(pmids))
    if not pmids:
        return {}

    summaries = fetch_pubmed_summaries(pmids)
    records = {str(s.get("Id") or s.get("uid")): s for s in summaries}

    # Save PubMed data
    pubmed_rows = [{
        "pmid": str(s.get("uid", s.get("Id", ""))),
        "title": str(s.get("Title", "")),
        "pubdate": str(s.get("PubDate", "")),
    } for s in summaries]
    save_csv(pubmed_rows, base_dir / "pubmed_records.csv")

    # Map to PMC
    pmcids, pmid_to_pmcid = pubmed_to_pmc(pmids)
    logger.info("PMCIDs mapped: %d", len(pmcids))

    # Process each article
    reports = [process_pmc_article(pmcid, base_dir, records, pmid_to_pmcid) for pmcid in pmcids]

    # Save summary
    save_json(reports, base_dir / "disease_summary.json")

    summary_rows = []
    for r in reports:
        row = {
            "pmcid": r.get("pmcid", ""),
            "status": r.get("status", ""),
            "kept": r.get("kept", False),
        }
        for k, v in r.get("data_counts", {}).items():
            row[f"data_{k}"] = v
        summary_rows.append(row)

    save_csv(summary_rows, base_dir / "disease_summary.csv")

    return {
        "disease_term": disease_term,
        "pmids": pmids,
        "pmcids": pmcids,
        "reports": reports,
    }


def run_all_diseases(email: str):
    """Download all diseases."""
    output_root = get_desktop_output_root()
    logger.info("Output root: %s\n", output_root)

    for name, query in DISEASE_QUERIES.items():
        logger.info("\n" + "=" * 70)
        logger.info("DISEASE: %s", name)
        logger.info("=" * 70)
        disease_out_dir = output_root / name
        fetch_ncbi_associated_data(query, str(disease_out_dir), email)


if __name__ == "__main__":
    EMAIL = "sebbahaya03@gmail.com"
    run_all_diseases(email=EMAIL)
