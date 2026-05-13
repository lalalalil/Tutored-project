import os
import shutil
import json
import time
import requests
import csv
from pathlib import Path
from bioservices import BioModels
import biomodels

# --- CONFIGURATION ---
RACINE_DATA = Path("./BioModels_Database_Final")
EXTENSIONS_METADATA = {'.png', '.jpg', '.jpeg', '.pdf', '.txt', '.docx', '.doc', '.xlsx', '.xls', '.csv'}
EXTENSIONS_MODELE = {
    '.xml', '.sbml', '.omex', '.sedml', '.cps', '.m', '.ode', 
    '.py', '.f', '.java', '.vcml', '.zip', '.tsv', '.yaml', '.graphml'
}

DISEASES = {
    "dengue_files": ('dengue','DENV'),
    "chikungunya_files": ("chikungunya", "CHIKV"),
    "lyme_files": ('lyme', 'borrelia', 'borreliosis'),
    "mpox_files": ('mpox', 'monkeypox'),
    "west_nile_files": ('west nile', 'WNV'),
    "influenza_files": ("influenza", "influenza virus","avian influenza", "H5N1"),
    "tuberculosis_files": ("tuberculosis", "TB", "mycobacterium"),
    "hiv_files": ("HIV", "human immunodeficiency virus"),
    "covid_files": ("covid", "SARS-CoV-2")
}

HEADERS = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
}

# --- PIPELINE FUNCTIONS ---

def step_download_and_sort():
    """Step 1: Download models and perform initial sorting (Model vs Metadata)"""
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
                d_metadata = dossier_base / "metadata"
                d_model = dossier_base / "model"
                d_metadata.mkdir(parents=True, exist_ok=True)
                d_model.mkdir(parents=True, exist_ok=True)

                files_to_preserve = [f"{model_id}_web_metadata.json", "files_list.json"]

                # 1. Fetch Web JSON Metadata
                web_url = f"https://www.ebi.ac.uk/biomodels/{model_id}?format=json"
                try:
                    web_res = requests.get(web_url, headers=HEADERS, timeout=20)
                    if web_res.status_code == 200:
                        with open(d_metadata / files_to_preserve[0], "w", encoding="utf-8") as f_web:
                            json.dump(web_res.json(), f_web, indent=4)
                except: pass

                # 2. Retrieve internal files list
                files_objects = biomodels.get_metadata(model_id)
                if not files_objects: continue
                with open(d_metadata / files_to_preserve[1], "w", encoding="utf-8") as f_meta:
                    json.dump(files_objects, f_meta, indent=4, default=str)

                # 3. Download & Sort
                for target in files_objects:
                    nom_reel = getattr(target, 'name', str(target))
                    if not nom_reel or nom_reel == "None": continue
                    
                    if nom_reel.lower().startswith("metadata."):
                        old_file = d_model / nom_reel
                        if old_file.exists(): old_file.unlink()

                    ext = Path(nom_reel).suffix.lower() or "no_ext"
                    is_meta = any(x in nom_reel.lower() for x in ["metadata", ".json", ".rdf", ".owl"]) or ext in EXTENSIONS_METADATA
                    dest = d_metadata / nom_reel if is_meta else d_model / nom_reel
                    if is_meta: files_to_preserve.append(nom_reel)

                    result = biomodels.get_file(target)
                    if isinstance(result, (str, Path)) and Path(result).exists():
                        shutil.copy(result, dest)
                
                print(f"  > {model_id} processed.")
                time.sleep(1.0)
            except Exception as e:
                print(f"Error {model_id}: {e}")

def step_clean_and_stats():
    """Step 2: Remove empty models and generate statistics CSV"""
    stats_globales = []
    # Using glob to find 'model' folders at any depth (in case already moved)
    for chemin_model in list(RACINE_DATA.glob("**/model")):
        dossier_model = chemin_model.parent
        a_un_modele = False
        comptage = {}
        
        for f in chemin_model.iterdir():
            if f.is_file():
                ext = f.suffix.lower() if f.suffix else "no_ext"
                comptage[ext] = comptage.get(ext, 0) + 1
                if ext in EXTENSIONS_MODELE: a_un_modele = True

        if not a_un_modele:
            print(f"[-] Deleting {dossier_model.name}: No modeling files.")
            shutil.rmtree(dossier_model)
        else:
            # Detect category from path
            row = {"category": dossier_model.parts[1], "model_id": dossier_model.name}
            row.update(comptage)
            stats_globales.append(row)

    if stats_globales:
        all_keys = set().union(*(d.keys() for d in stats_globales))
        extension_keys = all_keys - {'category', 'model_id'}
        cols = ['category', 'model_id'] + sorted(list(extension_keys))
        
        with open(RACINE_DATA / "extension_stats_summary.csv", "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=cols)
            writer.writeheader()
            for r in stats_globales:
                writer.writerow({c: r.get(c, 0) for c in cols})
        print("Cleaning and statistics completed.")

def step_separate_curation_status():
    """Step 5: Separate models into 'curated' and 'non_curated' subfolders under disease category"""
    print("Reorganizing: Disease/Status/Model_ID...")
    for category_dir in [d for d in RACINE_DATA.iterdir() if d.is_dir()]:
        # We find all metadata files in this disease folder
        for json_path in list(RACINE_DATA.glob(f"{category_dir.name}/**/*_web_metadata.json")):
            model_dir = json_path.parents[1]
            
            # Avoid moving if it's already in curated/non_curated
            if "curated" in model_dir.parts or "non_curated" in model_dir.parts:
                continue

            try:
                with open(json_path, 'r', encoding='utf-8') as f:
                    status = json.load(f).get("curationStatus", "UNKNOWN").lower()
                
                target_dir = category_dir / status
                target_dir.mkdir(exist_ok=True)
                
                # Physical move
                shutil.move(str(model_dir), str(target_dir / model_dir.name))
            except Exception as e:
                print(f"  ! Error moving {model_dir.name}: {e}")
    print("Curation status separation completed.")

def step_classify_by_approach():
    """Step 3: Reorganize folders by modeling approach (now supports status subfolders)"""
    data_export = []
    for category_dir in [d for d in RACINE_DATA.iterdir() if d.is_dir()]:
        for json_path in list(RACINE_DATA.glob(f"{category_dir.name}/**/*_web_metadata.json")):
            model_dir = json_path.parents[1]
            try:
                with open(json_path, 'r', encoding='utf-8') as f:
                    content = json.load(f)
                    approach = content.get("modellingApproach", {}).get("name", "Not_specified")
                    safe_name = approach.replace(" ", "_").replace("/", "-")
                    
                    data_export.append({
                        "disease_category": category_dir.name,
                        "model_id": model_dir.name,
                        "modelling_approach": approach
                    })

                # Move into a subfolder at the same level as the model folder
                new_parent = model_dir.parent / safe_name
                new_parent.mkdir(exist_ok=True)
                
                if model_dir.parent.name != safe_name:
                    shutil.move(str(model_dir), str(new_parent / model_dir.name))
            except: continue

    if data_export:
        with open(RACINE_DATA / "modelling_approaches_summary.csv", "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=data_export[0].keys())
            writer.writeheader()
            writer.writerows(data_export)
        print("Approach classification completed.")

def step_delete_pre_2015():
    """Step 4: Delete models published before 2015"""
    confirm = input("Confirm deletion of models before 2015? (yes/no): ")
    if confirm.lower() != 'yes': return

    for json_path in RACINE_DATA.glob("**/*_web_metadata.json"):
        model_dir = json_path.parents[1]
        try:
            with open(json_path, 'r', encoding='utf-8') as f:
                year = json.load(f).get("publication", {}).get("year")
                if year and int(year) < 2015:
                    print(f"Deleting {model_dir.name} ({year})")
                    shutil.rmtree(model_dir)
        except: continue

# --- MAIN MENU ---

def main():
    while True:
        print("\n--- BIOMODELS ANALYSIS PIPELINE ---")
        print("1. Download & Sort (Metadata vs Model)")
        print("2. Clean Empty Models & Generate Extension Stats CSV")
        print("3. Separate by Status (Disease/Status/Model_ID)")
        print("4. Classify Folders by Approach (Deep Classification)")
        print("5. Delete Models Published Before 2015")
        print("6. Exit")
        
        choice = input("\nSelect an option (1-6): ")
        if choice == '1': step_download_and_sort()
        elif choice == '2': step_clean_and_stats()
        elif choice == '3': step_separate_curation_status()
        elif choice == '4': step_classify_by_approach()
        elif choice == '5': step_delete_pre_2015()
        elif choice == '6': break
        else: print("Invalid choice.")

if __name__ == "__main__":
    main()