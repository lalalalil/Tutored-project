import os
import shutil
import json
import time
import requests
import csv
import subprocess
from pathlib import Path
from bioservices import BioModels
import biomodels

# --- CONFIGURATION ---
RACINE_DATA = Path("./BioModels_Database_Final")
STATS_OUTPUT_DIR = Path("./BioModels_Stats")          # [4] Single output dir for all count files
R_SCRIPT_PATH    = Path("./stats_biomodels.R")         # [5] Path to the R analysis script

EXTENSIONS_METADATA = {'.png', '.jpg', '.jpeg', '.pdf', '.txt', '.docx', '.doc', '.xlsx', '.xls', '.csv'}
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

# --- HELPER: DOI RESOLUTION ---

def fetch_doi(model_id: str, web_metadata: dict) -> str:
    """
    [2] Try to extract the DOI of the linked publication.
    Tries:  1) web_metadata["publication"]["doi"]
            2) web_metadata["publication"]["link"] if it contains 'doi.org'
            3) EBI Publications API fallback
    Returns the DOI string or an empty string.
    """
    pub = web_metadata.get("publication", {})

    # Direct field
    doi = pub.get("doi", "")
    if doi:
        return doi

    # Link field containing a DOI URL
    link = pub.get("link", "")
    if "doi.org" in link:
        return link.split("doi.org/")[-1]

    # Fallback: EBI Publications REST API
    try:
        pub_url = f"https://www.ebi.ac.uk/biomodels/{model_id}/publication"
        r = requests.get(pub_url, headers=HEADERS, timeout=15)
        if r.status_code == 200:
            data = r.json()
            doi = data.get("doi", "")
            if not doi:
                link = data.get("url", "")
                if "doi.org" in link:
                    doi = link.split("doi.org/")[-1]
            return doi
    except Exception:
        pass

    return ""


# --- PIPELINE FUNCTIONS ---

def step_download_and_sort():
    """
    Step 1: Download models and perform initial sorting (Model vs Metadata).

    Changes vs original:
      [1] After writing the consolidated JSON, delete every file in the
          metadata folder whose name contains 'metadata' (case-insensitive),
          so only the generated JSON survives.
      [2] The consolidated JSON now includes a 'doi' field fetched from the
          EBI API / publication record.
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

                # ---- Fetch web metadata ----
                web_metadata = {}
                web_url = f"https://www.ebi.ac.uk/biomodels/{model_id}?format=json"
                try:
                    web_res = requests.get(web_url, headers=HEADERS, timeout=20)
                    if web_res.status_code == 200:
                        web_metadata = web_res.json()
                except Exception:
                    pass

                # [2] Enrich with DOI
                doi = fetch_doi(model_id, web_metadata)
                web_metadata["doi"] = doi

                # ---- Fetch file list ----
                files_objects = biomodels.get_metadata(model_id)
                if not files_objects:
                    continue

                # [1] Build consolidated JSON (replaces separate metadata files)
                consolidated = {
                    "model_id":      model_id,
                    "doi":           doi,                          # [2]
                    "web_metadata":  web_metadata,
                    "files_list":    json.loads(
                        json.dumps(files_objects, default=str)
                    )
                }
                consolidated_path = d_metadata / f"{model_id}_metadata.json"
                with open(consolidated_path, "w", encoding="utf-8") as f_out:
                    json.dump(consolidated, f_out, indent=4)

                # [1] Delete every pre-existing file whose name contains 'metadata'
                #     (except the file we just created)
                for existing in list(d_metadata.iterdir()):
                    if existing == consolidated_path:
                        continue
                    if "metadata" in existing.name.lower():
                        existing.unlink()
                        print(f"    [-] Removed old metadata file: {existing.name}")

                # ---- Download model files ----
                for target in files_objects:
                    nom_reel = getattr(target, 'name', str(target))
                    if not nom_reel or nom_reel == "None":
                        continue

                    ext = Path(nom_reel).suffix.lower() or "no_ext"
                    is_meta = (
                        any(x in nom_reel.lower() for x in ["metadata", ".json", ".rdf", ".owl"])
                        or ext in EXTENSIONS_METADATA
                    )
                    dest = d_metadata / nom_reel if is_meta else d_model / nom_reel

                    result = biomodels.get_file(target)
                    if isinstance(result, (str, Path)) and Path(result).exists():
                        shutil.copy(result, dest)

                print(f"  > {model_id} processed.  DOI: {doi or 'N/A'}")
                time.sleep(1.0)

            except Exception as e:
                print(f"Error {model_id}: {e}")


def step_clean_metadata_dir():
    """
    [3] New step: for every metadata folder in the database, remove all files
    that are NOT the consolidated *_metadata.json generated in step 1.
    Asks for confirmation before deleting.
    """
    print("\n--- Clean Metadata Directories ---")
    print("This will delete every file in each 'metadata' folder that is NOT")
    print("the consolidated *_metadata.json file.")
    confirm = input("Proceed? (yes/no): ")
    if confirm.lower() != "yes":
        print("Cancelled.")
        return

    removed_total = 0
    for meta_dir in RACINE_DATA.glob("**/metadata"):
        if not meta_dir.is_dir():
            continue
        model_id = meta_dir.parent.name
        keep_name = f"{model_id}_metadata.json"
        for f in list(meta_dir.iterdir()):
            if f.is_file() and f.name != keep_name:
                f.unlink()
                removed_total += 1
                print(f"  [-] Removed: {f.relative_to(RACINE_DATA)}")

    print(f"\nDone. {removed_total} file(s) removed.")


def step_clean_and_stats():
    """
    Step 2: Remove empty models and generate statistics CSV.

    [4] All count/stats CSV files are written to STATS_OUTPUT_DIR instead of
        being scattered across RACINE_DATA.
    """
    STATS_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)   # [4]
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
            print(f"[-] Deleting {dossier_model.name}: No modeling files.")
            shutil.rmtree(dossier_model)
        else:
            row = {"category": dossier_model.parts[1], "model_id": dossier_model.name}
            row.update(comptage)
            stats_globales.append(row)

    if stats_globales:
        all_keys = set().union(*(d.keys() for d in stats_globales))
        cols     = ['category', 'model_id'] + sorted(list(all_keys - {'category', 'model_id'}))
        out_path = STATS_OUTPUT_DIR / "statistiques_modeles_nettoyes.csv"   # [4]
        with open(out_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=cols)
            writer.writeheader()
            for r in stats_globales:
                writer.writerow({c: r.get(c, 0) for c in cols})
        print(f"Extension stats saved to: {out_path}")

    print("Cleaning and statistics completed.")


def step_separate_curation_status():
    """Step 3: Separate models into 'curated' and 'non_curated' folders."""
    print("Reorganizing by Curation Status...")
    for category_dir in [d for d in RACINE_DATA.iterdir() if d.is_dir()]:
        for json_path in list(RACINE_DATA.glob(f"{category_dir.name}/**/*_metadata.json")):
            model_dir = json_path.parents[1]
            if any(x in model_dir.parts for x in ["curated", "non_curated"]):
                continue
            try:
                with open(json_path, 'r', encoding='utf-8') as f:
                    data       = json.load(f)
                    web_meta   = data.get("web_metadata", data)        # support both formats
                    status     = web_meta.get("curationStatus", "UNKNOWN").lower()
                target_dir = category_dir / status
                target_dir.mkdir(exist_ok=True)
                shutil.move(str(model_dir), str(target_dir / model_dir.name))
            except Exception as e:
                print(f"  ! Error moving {model_dir.name}: {e}")


def step_classify_by_approach():
    """
    Step 4: Reorganize folders by modeling approach.

    [4] modelling_approaches_summary.csv written to STATS_OUTPUT_DIR.
    """
    STATS_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)   # [4]
    data_export = []

    for json_path in list(RACINE_DATA.glob("**/*_metadata.json")):
        model_dir = json_path.parents[1]
        try:
            with open(json_path, 'r', encoding='utf-8') as f:
                data     = json.load(f)
                web_meta = data.get("web_metadata", data)
                approach = web_meta.get("modellingApproach", {}).get("name", "Not_specified")
                safe_name = approach.replace(" ", "_").replace("/", "-")

                data_export.append({
                    "disease_category":   model_dir.parts[1],
                    "model_id":           model_dir.name,
                    "modelling_approach": approach,
                    "doi":                data.get("doi", "")    # [2] carry DOI into CSV
                })

            new_parent = model_dir.parent / safe_name
            new_parent.mkdir(exist_ok=True)
            if model_dir.parent.name != safe_name:
                shutil.move(str(model_dir), str(new_parent / model_dir.name))
        except Exception:
            continue

    if data_export:
        out_path = STATS_OUTPUT_DIR / "modelling_approaches_summary.csv"   # [4]
        with open(out_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=data_export[0].keys())
            writer.writeheader()
            writer.writerows(data_export)
        print(f"Approach classification saved to: {out_path}")


def step_delete_pre_2015():
    """
    Step 5: Delete models published before 2015 and export comparative CSVs.

    [4] CSV files written to STATS_OUTPUT_DIR.
    """
    STATS_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)   # [4]
    all_models_data = []

    for json_path in RACINE_DATA.glob("**/*_metadata.json"):
        model_dir = json_path.parents[1]
        try:
            with open(json_path, 'r', encoding='utf-8') as f:
                data     = json.load(f)
                web_meta = data.get("web_metadata", data)
                year     = web_meta.get("publication", {}).get("year")
                status   = web_meta.get("curationStatus", "UNKNOWN")

                all_models_data.append({
                    "model_id": model_dir.name,
                    "year":     int(year) if year else 0,
                    "status":   status,
                    "category": model_dir.parts[1],
                    "doi":      data.get("doi", "")    # [2]
                })
        except Exception:
            continue

    if not all_models_data:
        print("No models found in the database.")
        return

    fieldnames = ["model_id", "year", "status", "category", "doi"]
    before_path = STATS_OUTPUT_DIR / "stats_BEFORE_2015_filter.csv"   # [4]
    with open(before_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_models_data)
    print(f"Generated {before_path} ({len(all_models_data)} models)")

    confirm = input("Confirm deletion of models published before 2015? (yes/no): ")
    if confirm.lower() != 'yes':
        return

    models_kept   = []
    deleted_count = 0

    for model in all_models_data:
        matching_dirs = list(RACINE_DATA.glob(f"**/{model['model_id']}"))
        if model['year'] < 2015 and model['year'] != 0:
            for d in matching_dirs:
                if d.is_dir():
                    shutil.rmtree(d)
            deleted_count += 1
        else:
            models_kept.append(model)

    after_path = STATS_OUTPUT_DIR / "stats_AFTER_2015_filter.csv"   # [4]
    with open(after_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(models_kept)

    print(f"Cleanup complete. Deleted: {deleted_count} | Remaining: {len(models_kept)}")


def step_generate_publication_dates_csv():
    """
    Step 6 : Generate publication_dates_summary.csv required by stats_biomodels.R.

    [4] Written to STATS_OUTPUT_DIR.
    [2] DOI column included.
    """
    STATS_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []

    for json_path in RACINE_DATA.glob("**/*_metadata.json"):
        model_dir = json_path.parents[1]
        try:
            with open(json_path, 'r', encoding='utf-8') as f:
                data     = json.load(f)
                web_meta = data.get("web_metadata", data)
                pub      = web_meta.get("publication", {})
                rows.append({
                    "model_id":         model_dir.name,
                    "disease_category": model_dir.parts[1],
                    "publication_year": pub.get("year", "N/A"),
                    "journal":          pub.get("journal", "N/A"),
                    "doi":              data.get("doi", "")    # [2]
                })
        except Exception:
            continue

    if rows:
        out_path = STATS_OUTPUT_DIR / "publication_dates_summary.csv"
        with open(out_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        print(f"Publication dates saved to: {out_path}")
    else:
        print("No data found.")


def step_run_r_stats():
    """
    [5] Run stats_biomodels.R using the CSV files in STATS_OUTPUT_DIR.
    Sets the working directory to STATS_OUTPUT_DIR so the R script finds
    all three expected CSVs without modification:
      - publication_dates_summary.csv
      - modelling_approaches_summary.csv
      - statistiques_modeles_nettoyes.csv
    PNG outputs are written to STATS_OUTPUT_DIR/outputs/ by the R script.
    """
    if not R_SCRIPT_PATH.exists():
        print(f"R script not found: {R_SCRIPT_PATH.resolve()}")
        print("Place stats_biomodels.R next to this script and retry.")
        return

    if not STATS_OUTPUT_DIR.exists():
        print(f"Stats directory not found: {STATS_OUTPUT_DIR.resolve()}")
        print("Run steps 2 and 4 first to generate the CSV files.")
        return

    missing = []
    for csv_name in ["publication_dates_summary.csv",
                     "modelling_approaches_summary.csv",
                     "statistiques_modeles_nettoyes.csv"]:
        if not (STATS_OUTPUT_DIR / csv_name).exists():
            missing.append(csv_name)
    if missing:
        print("Missing CSV file(s) in stats directory:")
        for m in missing:
            print(f"  - {m}")
        print("Run the relevant pipeline steps first.")
        return

    print(f"\nRunning {R_SCRIPT_PATH.name} on: {STATS_OUTPUT_DIR.resolve()}")
    print("(PNG outputs will appear in the 'outputs' subfolder)\n")

    try:
        result = subprocess.run(
            ["Rscript", str(R_SCRIPT_PATH.resolve())],
            cwd=str(STATS_OUTPUT_DIR.resolve()),   # R's working dir = stats folder
            capture_output=False
        )
        if result.returncode == 0:
            print("\nR analysis completed successfully.")
        else:
            print(f"\nR script exited with code {result.returncode}.")
    except FileNotFoundError:
        print("Rscript not found. Make sure R is installed and in your PATH.")


# --- MAIN MENU ---

def main():
    while True:
        print("\n--- BIOMODELS ANALYSIS PIPELINE ---")
        print(f"  Stats output dir : {STATS_OUTPUT_DIR.resolve()}")
        print(f"  R script         : {R_SCRIPT_PATH.resolve()}")
        print()
        print("1.  Download & Sort (Metadata vs Model)")
        print("2.  Clean Empty Models & Generate Extension Stats CSV")
        print("3.  Separate by Curation Status")
        print("4.  Classify Folders by Modelling Approach")
        print("5.  Delete Models Published Before 2015 (With CSV Stats)")
        print("6.  Generate Publication Dates CSV")
        print("--- Other Options ---")
        print("7.  Clean Metadata Dirs (keep only *_metadata.json)")
        print("8.  Run R Statistical Analysis (stats_biomodels.R)")
        print("---")
        print("9.  Exit")

        choice = input("\nSelect an option (1-9): ").strip()

        if   choice == '1': step_download_and_sort()
        elif choice == '2': step_clean_and_stats()
        elif choice == '3': step_separate_curation_status()
        elif choice == '4': step_classify_by_approach()
        elif choice == '5': step_delete_pre_2015()
        elif choice == '6': step_generate_publication_dates_csv()
        elif choice == '7': step_clean_metadata_dir()
        elif choice == '8': step_run_r_stats()
        elif choice == '9': break
        else: print("Invalid choice.")

if __name__ == "__main__":
    main()
