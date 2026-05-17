import os
import shutil
import json
import time
import datetime
import requests
import csv
import subprocess
from pathlib import Path
from bioservices import BioModels
import biomodels

# =============================================================================
# CONFIGURATION
# =============================================================================

RACINE_DATA      = Path("./BioModels_Database_Final")
STATS_OUTPUT_DIR = Path("./BioModels_Stats")   # all CSV/PNG outputs go here
R_SCRIPT_PATH    = Path("./stats_biomodels.R") # path to the R analysis script

EXTENSIONS_METADATA = {
    '.png', '.jpg', '.jpeg', '.pdf', '.txt',
    '.docx', '.doc', '.xlsx', '.xls', '.csv'
}
EXTENSIONS_MODELE = {
    '.xml', '.sbml', '.omex', '.sedml', '.cps', '.m', '.ode',
    '.py', '.f', '.java', '.vcml', '.zip', '.tsv', '.yaml', '.graphml'
}

DISEASES = {
    "dengue_files":       ('dengue', 'DENV'),
    "chikungunya_files":  ("chikungunya", "CHIKV"),
    "mpox_files":         ('mpox', 'monkeypox'),
    "west_nile_files":    ('west nile', 'WNV'),
    "influenza_files":    ("influenza", "influenza virus", "avian influenza", "H5N1"),
    "tuberculosis_files": ("tuberculosis", "TB", "mycobacterium"),
    "hiv_files":          ("HIV", "human immunodeficiency virus"),
    "covid_files":        ("covid", "SARS-CoV-2")
}

HEADERS = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 '
                  '(KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
}

# =============================================================================
# HELPER — DOI RESOLUTION (4-level cascade)
# =============================================================================

def fetch_doi(model_id: str, web_metadata: dict) -> str:
    """
    Extract the DOI of the linked publication using 4 strategies:

    1. web_metadata["publication"]["doi"]         — curated models (BIOMD)
    2. web_metadata["publication"]["link"]        — if it contains doi.org
    3. web_metadata["modelLevelAnnotations"]      — bqmodel:isDescribedBy
       (non-curated models with an article but no publication block)
    4. EBI Publications REST API                  — last resort

    Returns the raw DOI string (e.g. "10.1371/journal.pcbi.1000282") or "".
    """
    # Strategy 1 & 2 — publication block
    pub  = web_metadata.get("publication", {})
    doi  = pub.get("doi", "").strip()
    if doi:
        return doi

    link = pub.get("link", "").strip()
    if "doi.org" in link:
        return link.split("doi.org/")[-1].strip()

    # Strategy 3 — modelLevelAnnotations
    for annotation in web_metadata.get("modelLevelAnnotations", []):
        if annotation.get("qualifier") == "bqmodel:isDescribedBy":
            uri = annotation.get("uri", "")
            if "doi.org" in uri:
                return uri.split("doi.org/")[-1].strip()
            accession = annotation.get("accession", "")
            if accession.startswith("10."):
                return accession.strip()

    # Strategy 4 — EBI Publications REST API
    try:
        r = requests.get(
            f"https://www.ebi.ac.uk/biomodels/{model_id}/publication",
            headers=HEADERS, timeout=15
        )
        if r.status_code == 200:
            data = r.json()
            doi  = data.get("doi", "").strip()
            if doi:
                return doi
            url = data.get("url", "").strip()
            if "doi.org" in url:
                return url.split("doi.org/")[-1].strip()
    except Exception:
        pass

    return ""


# =============================================================================
# HELPER — YEAR RESOLUTION (3-level cascade)
# =============================================================================

def fetch_year(web_metadata: dict, doi: str) -> int | None:
    """
    Extract the publication year using 3 strategies in order of preference:

    1. web_metadata["publication"]["year"]   — explicit year in BioModels JSON
       (always present for curated BIOMD models)
    2. DOI → CrossRef API                   — article year for non-curated
       models that have a DOI but no publication block
    3. web_metadata["firstPublished"]        — Unix timestamp of BioModels
       deposit; last resort, this is NOT the article year but the deposit date

    Returns an int (year) or None if nothing could be found.
    The "year_source" key in the returned dict tells you which strategy worked.
    """
    # Strategy 1 — publication block
    year = web_metadata.get("publication", {}).get("year")
    if year:
        try:
            return int(year), "publication_block"
        except (ValueError, TypeError):
            pass

    # Strategy 2 — CrossRef via DOI
    if doi:
        try:
            r = requests.get(
                f"https://api.crossref.org/works/{doi}",
                headers={"User-Agent": "BioModelsPipeline/1.0 (mailto:your@email.com)"},
                timeout=15
            )
            if r.status_code == 200:
                parts = (
                    r.json()
                    .get("message", {})
                    .get("published", {})
                    .get("date-parts", [[]])[0]
                )
                if parts:
                    return int(parts[0]), "crossref"
        except Exception:
            pass

    # Strategy 3 — firstPublished timestamp (deposit date on BioModels)
    ts = web_metadata.get("firstPublished")
    if ts:
        try:
            return datetime.datetime.fromtimestamp(int(ts)).year, "firstPublished"
        except Exception:
            pass

    return None, "not_found"


# =============================================================================
# STEP 1 — DOWNLOAD & SORT
# =============================================================================

def step_download_and_sort():
    """
    Download models from BioModels and sort files into model/ and metadata/.

    - Writes a single consolidated {model_id}_metadata.json per model:
        model_id, doi, year, year_source, web_metadata, files_list
    - Deletes any pre-existing file whose name contains 'metadata' (except
      the newly generated JSON).
    - DOI via 4-level cascade, year via 3-level cascade.
    """
    RACINE_DATA.mkdir(parents=True, exist_ok=True)
    s = BioModels()

    for folder_name, queries in DISEASES.items():
        print(f"\n--- Category: {folder_name} ---")
        query_string = " OR ".join(f'"{q}"' for q in queries)
        results = s.search(query_string, numResults=100)

        if not isinstance(results, dict):
            time.sleep(5)
            continue

        ids = [res['id'] for res in results.get('models', [])]

        for model_id in ids:
            try:
                dossier_base = RACINE_DATA / folder_name / model_id
                d_metadata   = dossier_base / "metadata"
                d_model      = dossier_base / "model"
                d_metadata.mkdir(parents=True, exist_ok=True)
                d_model.mkdir(parents=True, exist_ok=True)

                # Fetch full web metadata JSON
                web_metadata = {}
                try:
                    web_res = requests.get(
                        f"https://www.ebi.ac.uk/biomodels/{model_id}?format=json",
                        headers=HEADERS, timeout=20
                    )
                    if web_res.status_code == 200:
                        web_metadata = web_res.json()
                except Exception:
                    pass

                # Resolve DOI and year
                doi              = fetch_doi(model_id, web_metadata)
                year, year_src   = fetch_year(web_metadata, doi)

                # Fetch file list
                files_objects = biomodels.get_metadata(model_id)
                if not files_objects:
                    continue

                # Write consolidated metadata JSON
                consolidated = {
                    "model_id":    model_id,
                    "doi":         doi,
                    "year":        year,
                    "year_source": year_src,
                    "web_metadata": web_metadata,
                    "files_list":  json.loads(json.dumps(files_objects, default=str))
                }
                consolidated_path = d_metadata / f"{model_id}_metadata.json"
                with open(consolidated_path, "w", encoding="utf-8") as f_out:
                    json.dump(consolidated, f_out, indent=4)

                # Remove stale metadata files
                for existing in list(d_metadata.iterdir()):
                    if existing == consolidated_path:
                        continue
                    if "metadata" in existing.name.lower():
                        existing.unlink()
                        print(f"    [-] Removed: {existing.name}")

                # Download and sort model files
                for target in files_objects:
                    nom_reel = getattr(target, 'name', str(target))
                    if not nom_reel or nom_reel == "None":
                        continue
                    ext     = Path(nom_reel).suffix.lower() or "no_ext"
                    is_meta = (
                        any(x in nom_reel.lower() for x in ["metadata", ".json", ".rdf", ".owl"])
                        or ext in EXTENSIONS_METADATA
                    )
                    dest   = d_metadata / nom_reel if is_meta else d_model / nom_reel
                    result = biomodels.get_file(target)
                    if isinstance(result, (str, Path)) and Path(result).exists():
                        shutil.copy(result, dest)

                print(f"  > {model_id}  DOI: {doi or 'N/A'}  "
                      f"Year: {year or 'N/A'} [{year_src}]")
                time.sleep(1.0)

            except Exception as e:
                print(f"Error {model_id}: {e}")


# =============================================================================
# STEP 2 — CLEAN EMPTY MODELS & GENERATE EXTENSION STATS CSV
# =============================================================================

def step_clean_and_stats():
    """
    Remove model folders with no recognised modelling files, then write
    statistiques_modeles_nettoyes.csv to STATS_OUTPUT_DIR.
    """
    STATS_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    stats_globales = []

    for chemin_model in list(RACINE_DATA.glob("**/model")):
        dossier_model = chemin_model.parent
        a_un_modele   = False
        comptage      = {}

        for f in chemin_model.iterdir():
            if f.is_file():
                ext = f.suffix.lower() if f.suffix else "no_ext"
                comptage[ext] = comptage.get(ext, 0) + 1
                if ext in EXTENSIONS_MODELE:
                    a_un_modele = True

        if not a_un_modele:
            print(f"[-] Deleting {dossier_model.name}: no modelling files.")
            shutil.rmtree(dossier_model)
        else:
            row = {"category": dossier_model.parts[1], "model_id": dossier_model.name}
            row.update(comptage)
            stats_globales.append(row)

    if stats_globales:
        all_keys = set().union(*(d.keys() for d in stats_globales))
        cols     = ['category', 'model_id'] + sorted(list(all_keys - {'category', 'model_id'}))
        out_path = STATS_OUTPUT_DIR / "statistiques_modeles_nettoyes.csv"
        with open(out_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=cols)
            writer.writeheader()
            for r in stats_globales:
                writer.writerow({c: r.get(c, 0) for c in cols})
        print(f"Extension stats saved to: {out_path}")

    print("Cleaning and statistics completed.")


# =============================================================================
# STEP 3 — SEPARATE BY CURATION STATUS
# =============================================================================

def step_separate_curation_status():
    """Reorganise model folders into disease/curated/ and disease/non_curated/."""
    print("Reorganising by curation status...")
    for category_dir in [d for d in RACINE_DATA.iterdir() if d.is_dir()]:
        for json_path in list(RACINE_DATA.glob(f"{category_dir.name}/**/*_metadata.json")):
            model_dir = json_path.parents[1]
            if any(x in model_dir.parts for x in ["curated", "non_curated"]):
                continue
            try:
                with open(json_path, 'r', encoding='utf-8') as f:
                    data     = json.load(f)
                    web_meta = data.get("web_metadata", data)
                    status   = web_meta.get("curationStatus", "UNKNOWN").lower()
                target_dir = category_dir / status
                target_dir.mkdir(exist_ok=True)
                shutil.move(str(model_dir), str(target_dir / model_dir.name))
            except Exception as e:
                print(f"  ! Error moving {model_dir.name}: {e}")


# =============================================================================
# STEP 4 — CLASSIFY BY MODELLING APPROACH
# =============================================================================

def step_classify_by_approach():
    """
    Reorganise folders by modelling approach and write
    modelling_approaches_summary.csv to STATS_OUTPUT_DIR.
    """
    STATS_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    data_export = []

    for json_path in list(RACINE_DATA.glob("**/*_metadata.json")):
        model_dir = json_path.parents[1]
        try:
            with open(json_path, 'r', encoding='utf-8') as f:
                data      = json.load(f)
                web_meta  = data.get("web_metadata", data)
                approach  = web_meta.get("modellingApproach", {}).get("name", "Not_specified")
                safe_name = approach.replace(" ", "_").replace("/", "-")

                data_export.append({
                    "disease_category":   model_dir.parts[1],
                    "model_id":           model_dir.name,
                    "modelling_approach": approach,
                    "doi":                data.get("doi", ""),
                    "year":               data.get("year", ""),
                    "year_source":        data.get("year_source", "")
                })

            new_parent = model_dir.parent / safe_name
            new_parent.mkdir(exist_ok=True)
            if model_dir.parent.name != safe_name:
                shutil.move(str(model_dir), str(new_parent / model_dir.name))
        except Exception:
            continue

    if data_export:
        out_path = STATS_OUTPUT_DIR / "modelling_approaches_summary.csv"
        with open(out_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=data_export[0].keys())
            writer.writeheader()
            writer.writerows(data_export)
        print(f"Approach classification saved to: {out_path}")


# =============================================================================
# STEP 5 — DELETE MODELS PUBLISHED BEFORE 2015
# =============================================================================

def step_delete_pre_2015():
    """
    Export stats_BEFORE_2015_filter.csv, ask for confirmation, delete old
    models, then export stats_AFTER_2015_filter.csv.
    Uses the pre-resolved year (and year_source) stored in the JSON.
    """
    STATS_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    all_models_data = []

    for json_path in RACINE_DATA.glob("**/*_metadata.json"):
        model_dir = json_path.parents[1]
        try:
            with open(json_path, 'r', encoding='utf-8') as f:
                data     = json.load(f)
                web_meta = data.get("web_metadata", data)
                status   = web_meta.get("curationStatus", "UNKNOWN")

                # Use pre-resolved year; fall back to live resolution
                year     = data.get("year")
                year_src = data.get("year_source", "")
                if not year:
                    year, year_src = fetch_year(web_meta, data.get("doi", ""))

                all_models_data.append({
                    "model_id":    model_dir.name,
                    "year":        year or 0,
                    "year_source": year_src,
                    "status":      status,
                    "category":    model_dir.parts[1],
                    "doi":         data.get("doi", "")
                })
        except Exception:
            continue

    if not all_models_data:
        print("No models found in the database.")
        return

    fieldnames  = ["model_id", "year", "year_source", "status", "category", "doi"]
    before_path = STATS_OUTPUT_DIR / "stats_BEFORE_2015_filter.csv"
    with open(before_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_models_data)
    print(f"Generated {before_path} ({len(all_models_data)} models)")

    # Show breakdown by year_source so the user knows where data came from
    from collections import Counter
    src_counts = Counter(m["year_source"] for m in all_models_data)
    print("  Year sources:")
    for src, count in src_counts.most_common():
        print(f"    {src}: {count}")

    confirm = input("\nConfirm deletion of models with year < 2015? (yes/no): ")
    if confirm.lower() != 'yes':
        return

    models_kept   = []
    deleted_count = 0
    for model in all_models_data:
        matching_dirs = list(RACINE_DATA.glob(f"**/{model['model_id']}"))
        if model['year'] and model['year'] < 2015:
            for d in matching_dirs:
                if d.is_dir():
                    shutil.rmtree(d)
            deleted_count += 1
        else:
            models_kept.append(model)

    after_path = STATS_OUTPUT_DIR / "stats_AFTER_2015_filter.csv"
    with open(after_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(models_kept)

    print(f"Cleanup complete. Deleted: {deleted_count} | Remaining: {len(models_kept)}")


# =============================================================================
# STEP 6 — GENERATE PUBLICATION DATES CSV
# =============================================================================

def step_generate_publication_dates_csv():
    """
    Scan all *_metadata.json files and write publication_dates_summary.csv
    to STATS_OUTPUT_DIR. Required by stats_biomodels.R (Block 1).

    Year resolution order per model:
      1. pre-resolved year stored in the JSON  (from step 1)
      2. CrossRef via DOI                      (live call if missing)
      3. firstPublished timestamp              (deposit date fallback)
    The year_source column records which strategy was used.
    """
    STATS_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    rows        = []
    missing_log = []   # models where year could not be resolved

    for json_path in RACINE_DATA.glob("**/*_metadata.json"):
        model_dir = json_path.parents[1]
        try:
            with open(json_path, 'r', encoding='utf-8') as f:
                data     = json.load(f)
                web_meta = data.get("web_metadata", data)
                pub      = web_meta.get("publication", {})
                doi      = data.get("doi", "")

                # Year — use stored value or resolve live
                year     = data.get("year")
                year_src = data.get("year_source", "")
                if not year:
                    year, year_src = fetch_year(web_meta, doi)

                if not year:
                    missing_log.append(model_dir.name)

                rows.append({
                    "model_id":         model_dir.name,
                    "disease_category": model_dir.parts[1],
                    "publication_year": year or "N/A",
                    "year_source":      year_src,
                    "journal":          pub.get("journal", "N/A"),
                    "doi":              doi
                })
        except Exception:
            continue

    if rows:
        out_path = STATS_OUTPUT_DIR / "publication_dates_summary.csv"
        with open(out_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        print(f"Publication dates saved to: {out_path} ({len(rows)} entries)")

        if missing_log:
            print(f"\nWarning: {len(missing_log)} model(s) with no year resolved:")
            for m in missing_log:
                print(f"  - {m}")
    else:
        print("No data found.")


# =============================================================================
# STEP 7 — CLEAN METADATA DIRS (keep only *_metadata.json)
# =============================================================================

def step_clean_metadata_dir():
    """
    For every metadata/ folder, remove all files that are NOT the consolidated
    {model_id}_metadata.json. Asks for confirmation first.
    """
    print("\n--- Clean Metadata Directories ---")
    print("This will delete every file in each 'metadata/' folder that is NOT")
    print("the consolidated *_metadata.json file.")
    confirm = input("Proceed? (yes/no): ")
    if confirm.lower() != "yes":
        print("Cancelled.")
        return

    removed_total = 0
    for meta_dir in RACINE_DATA.glob("**/metadata"):
        if not meta_dir.is_dir():
            continue
        model_id  = meta_dir.parent.name
        keep_name = f"{model_id}_metadata.json"
        for f in list(meta_dir.iterdir()):
            if f.is_file() and f.name != keep_name:
                f.unlink()
                removed_total += 1
                print(f"  [-] Removed: {f.relative_to(RACINE_DATA)}")

    print(f"\nDone. {removed_total} file(s) removed.")


# =============================================================================
# STEP 8 — RUN R STATISTICAL ANALYSIS
# =============================================================================

def step_run_r_stats():
    """
    Run stats_biomodels.R with STATS_OUTPUT_DIR as the working directory.
    Required CSVs (all generated by this pipeline):
      - publication_dates_summary.csv       (step 6)
      - modelling_approaches_summary.csv    (step 4)
      - statistiques_modeles_nettoyes.csv   (step 2)
      - stats_BEFORE_2015_filter.csv        (step 5)
    PNG outputs are written to STATS_OUTPUT_DIR/outputs/ by the R script.
    """
    if not R_SCRIPT_PATH.exists():
        print(f"R script not found: {R_SCRIPT_PATH.resolve()}")
        print("Place stats_biomodels.R next to this script and retry.")
        return

    if not STATS_OUTPUT_DIR.exists():
        print(f"Stats directory not found: {STATS_OUTPUT_DIR.resolve()}")
        return

    required = [
        "publication_dates_summary.csv",
        "modelling_approaches_summary.csv",
        "statistiques_modeles_nettoyes.csv",
        "stats_BEFORE_2015_filter.csv"
    ]
    missing = [c for c in required if not (STATS_OUTPUT_DIR / c).exists()]
    if missing:
        print("Missing CSV file(s):")
        for m in missing:
            print(f"  - {m}")
        return

    print(f"\nRunning {R_SCRIPT_PATH.name}")
    print(f"Working directory: {STATS_OUTPUT_DIR.resolve()}")
    print("PNG outputs → outputs/ subfolder\n")

    try:
        result = subprocess.run(
            ["Rscript", str(R_SCRIPT_PATH.resolve())],
            cwd=str(STATS_OUTPUT_DIR.resolve())
        )
        if result.returncode == 0:
            print("\nR analysis completed successfully.")
        else:
            print(f"\nR script exited with code {result.returncode}.")
    except FileNotFoundError:
        print("Rscript not found. Make sure R is installed and in your PATH.")


# =============================================================================
# MAIN MENU
# =============================================================================

def main():
    while True:
        print("\n========================================")
        print("     BIOMODELS ANALYSIS PIPELINE        ")
        print("========================================")
        print(f"  Data dir    : {RACINE_DATA.resolve()}")
        print(f"  Stats dir   : {STATS_OUTPUT_DIR.resolve()}")
        print(f"  R script    : {R_SCRIPT_PATH.resolve()}")
        print("----------------------------------------")
        print("1.  Download & Sort  (metadata vs model)")
        print("2.  Clean Empty Models & Extension Stats CSV")
        print("3.  Separate by Curation Status")
        print("4.  Classify by Modelling Approach")
        print("5.  Generate Publication Dates CSV")
        print("6.  Delete Models Before 2015  (+ CSV export)")
        print("7.  Clean Metadata Dirs  (keep *_metadata.json only)")
        print("8.  Run R Statistical Analysis  (stats_biomodels.R)")
        print("----------------------------------------")
        print("9.  Exit")

        choice = input("\nSelect an option (1-9): ").strip()

        if   choice == '1': step_download_and_sort()
        elif choice == '2': step_clean_and_stats()
        elif choice == '3': step_separate_curation_status()
        elif choice == '4': step_classify_by_approach()
        elif choice == '5': step_generate_publication_dates_csv()
        elif choice == '6': step_delete_pre_2015()
        elif choice == '7': step_clean_metadata_dir()
        elif choice == '8': step_run_r_stats()
        elif choice == '9': break
        else: print("Invalid choice.")

if __name__ == "__main__":
    main()