import requests
import json
import time
import os
from pathlib import Path

def recherche_zenodo(query, size=20): 
    url = "https://zenodo.org/api/records" # Querying Zenodo API
    filtered_hits = [] # Store retrieved hits
    current_page = 1 # Start at page 1

    # Search parameters
    params = {
        "q": 'title :'+query,
        "size": size,
        "all_versions": "false", # Filter to avoid duplicate versions1
        "type": "software", # Filter by type ("publication", "dataset", "software")
        "status": "published", 
        "sort": "mostrecent", # Sort by most recent first ("mostviewed", "-publicationdate")
    }

    # Loop to handle pagination
    while len(filtered_hits) < size:
        params["page"] = current_page
        response = requests.get(url, params=params) # Execute the request

        if response.status_code == 200: # Success code
          data = response.json()
          raw_hits = data["hits"]["hits"]

          if not raw_hits: # If the page is empty, stop the loop
              break

          # Filter results by file extension
          new_hits = [
              h for h in raw_hits
              if any(f["key"].endswith((".xml", ".sbml", ".omex", ".R")) for f in h.get("files", []))
          ]

          filtered_hits.extend(new_hits)
          current_page += 1 # Move to the next page for the next iteration
          time.sleep(1) # Pause between pages to avoid overloading Zenodo
        else :
          print(f"Erreur API: {response.status_code}")

    # Return only the requested number of results
    final_results = filtered_hits[:size]
    print(f"{len(final_results)} modèles valides récupérés pour {query}")
    return final_results

def file_extraction(doi):
    url = f"https://zenodo.org/api/records/{doi}" # Connection to the URL
    response = requests.get(url) # Request

    if response.status_code == 200: # Success code

        # Directory path creation
        dossier_article = Path(f"DOI{doi}") # Creating the various paths
        dossier_metadata = dossier_article / "metadata"
        dossier_data = dossier_article / "data"

        # If the directory for this DOI doesn't exist, extract the data
        if not dossier_article.exists():

            dossier_metadata.mkdir(parents=True, exist_ok=True) # Creating directories
            dossier_data.mkdir(parents=True, exist_ok=True)

            data = response.json() # Retrieving metadata

            metadata_file = dossier_metadata / f"metadata_{doi}.json" # Create the metadata file
            with open(metadata_file, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=4, ensure_ascii=False) # json.dump converts Python dict to formatted JSON text

            files = data.get("files", []) # Downloading model files
            for f in files:
                download_url = f["links"]["self"]
                filename = dossier_data / f["key"] # Writing directly into the data sub-directory
                time.sleep(1)

                r = requests.get(download_url)
                with open(filename, "wb") as local_file:
                    local_file.write(r.content)

def extraction_finale(liste):
    racine = Path("/home/bravais/Documents/Ptut/Models") # Create the root directory
    racine.mkdir(parents=True, exist_ok=True)
    os.chdir(racine) # Move into the root directory

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
    ["\"West nile\"", "WNV"],
    ["Influenza", "\"influenza virus\"", "\"avian influenza\"", "H5N1"],
    ["Tuberculosis", "TB", "mycobacterium"],
    ["HIV", "\"human immunodeficiency virus\""]
])