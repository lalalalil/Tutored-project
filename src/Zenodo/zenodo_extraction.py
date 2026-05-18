import requests
import json
import time
import os
from pathlib import Path
import shutil

# File type definitions 
EXECUTABLE_EXTENSIONS = (
    ".py", ".r", ".m", ".ipynb", ".sbml", ".sedml", 
    ".omex", ".cps", ".cellml", ".jl", ".cpp", ".c"
)

# Technical keywords 
KEYWORDS_QUERY = ' OR '.join([
    "SBML", "OMEX", "notebook", "model", '"source code"', "python", '"R script"'
])

def recherche_zenodo(query, size=20): 
    url = "https://zenodo.org/api/records" # Querying Zenodo API
    filtered_hits = [] # Store retrieved hits
    current_page = 1 # Start at page 1

    # Combine virus name with technical keywords to filter results
    full_q = f'title:\"{query.strip("\"")}\" AND ({KEYWORDS_QUERY}) -title:epidemic* -title:transmission -title:climate -title:bed'

    # Search parameters
    params = {
        "q": full_q,
        "size": size,
        "all_versions": "false", # Filter to avoid duplicate versions
        "type": "software", # Filter by type ("publication", "dataset", "software")
        "status": "published", 
        "sort": "mostrecent", # Sort by most recent first
    }

    # Loop to handle pagination
    while len(filtered_hits) < size and current_page <= 5:
        params["page"] = current_page
        response = requests.get(url, params=params) # Execute the request

        if response.status_code == 200: # Success code
          data = response.json()
          raw_hits = data["hits"]["hits"]

          if not raw_hits: # If the page is empty, stop the loop
              break

          # Filter results by file extension (case insensitive)
          new_hits = [
              h for h in raw_hits
              if any(f["key"].lower().endswith(EXECUTABLE_EXTENSIONS) for f in h.get("files", []))
          ]

          filtered_hits.extend(new_hits)
          current_page += 1 # Move to the next page for the next iteration
          time.sleep(1) # Pause between pages to avoid overloading Zenodo
        else :
          print(f"Erreur API: {response.status_code}")
          break

    # Return only the requested number of results
    final_results = filtered_hits[:size]
    print(f"{len(final_results)} valid models retrieved for {query}")
    return final_results

def file_extraction(doi):
    url = f"https://zenodo.org/api/records/{doi}" # Connection to the URL
    response = requests.get(url) # Request

    if response.status_code == 200: # Success code
        data = response.json() # Retrieving metadata
        
        # Identify valid files based on the synchronized extension list
        files_to_download = [
            f for f in data.get("files", []) 
            if f["key"].lower().endswith(EXECUTABLE_EXTENSIONS)
        ]

        # Only proceed if there are valid files to download
        if files_to_download:
            # Directory path creation
            dossier_article = Path(f"DOI{doi}") # Creating the various paths
            dossier_metadata = dossier_article / "metadata"
            dossier_data = dossier_article / "data"

            # If the directory for this DOI doesn't exist, extract the data
            if not dossier_article.exists():
                dossier_metadata.mkdir(parents=True, exist_ok=True) # Creating directories
                dossier_data.mkdir(parents=True, exist_ok=True)

                # Create the metadata file
                metadata_file = dossier_metadata / f"metadata_{doi}.json"
                with open(metadata_file, "w", encoding="utf-8") as f:
                    json.dump(data, f, indent=4, ensure_ascii=False) # Converts Python dict to formatted JSON text

                # Downloading model files
                for f in files_to_download:
                    download_url = f["links"]["self"]
                    filename = dossier_data / f["key"] # Writing directly into the data sub-directory
                    
                    time.sleep(1)
                    r = requests.get(download_url)
                    with open(filename, "wb") as local_file:
                        local_file.write(r.content)


def extraction_finale(liste):
    racine = Path("/home/bravais/Documents/P_tut/Models")
    racine.mkdir(parents=True, exist_ok=True)
    os.chdir(racine)

    # Iterate through each virus in our list
    for virus in liste: 
        nom_dossier = virus[0] # Create a directory for each virus
        dossier_virus = Path(nom_dossier)
        dossier_virus.mkdir(parents=True, exist_ok=True)
        os.chdir(dossier_virus) # Move into the virus directory

        dossier_zenodo = racine / dossier_virus / "Zenodo"
        dossier_zenodo.mkdir(parents=True, exist_ok=True)
        os.chdir(dossier_zenodo)
        
        for autre_nom in virus: 
            hits = recherche_zenodo(autre_nom, 20)
            for results in hits:
                doi = results["links"]["doi"].split(".")[-1] # Extract DOI from the hit
                file_extraction(doi) # Retrieve files using the DOI
        os.chdir("..") 
        time.sleep(1)
        
        os.chdir(racine)
        
        # Remove empty directories
        if dossier_zenodo.exists() and not any(dossier_zenodo.iterdir()):
            shutil.rmtree(dossier_virus)


extraction_finale([
    ["Dengue", "DENV"],
    ["Chikungunya", "CHIKV"],
    ["Mpox", "monkeypox"],
    ['West nile', "WNV"],
    ["Influenza", '"influenza virus"', '"avian influenza"', "H5N1"],
    ["Tuberculosis", "mycobacterium"],
    ["HIV", '"human immunodeficiency virus"'],
    ["Covid", "SARS-CoV-2"]
])


