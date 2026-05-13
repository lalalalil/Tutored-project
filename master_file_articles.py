import os
import requests
import time
from Bio import Entrez, Medline

# --- CONFIGURATION ---
Entrez.email = "ton.nom@etudiant.univ.fr" 

def get_citations_from_doi(doi):
    """Récupère le nombre de citations officiel via l'API Crossref."""
    if not doi or "no_doi" in doi:
        return "0"
    try:
        url = f"https://api.crossref.org/works/{doi}"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            return str(data['message'].get('is-referenced-by-count', 0))
    except:
        pass
    return "0"

def fetch_pubmed_with_counter(disease, limit=500):
    output_dir = f"data_{disease.lower()}_pubmed"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print(f"--- Recherche PubMed : {disease} ---")
    
    query = f"({disease}[Title] AND model*[Title]) AND English[Language] AND 2015:2026[DP]"
    
    # 1. Recherche des IDs
    search_handle = Entrez.esearch(db="pubmed", term=query, retmax=limit)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    
    ids = search_results["IdList"]
    if not ids:
        print("Aucun article trouvé.")
        return

    # 2. Récupération des données
    fetch_handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
    records = list(Medline.parse(fetch_handle))
    fetch_handle.close()

    # --- INITIALISATION DU COMPTEUR ---
    files_downloaded = 0
    total_to_fetch = len(records)
    print(f"Préparation de l'extraction pour {total_to_fetch} articles...\n")

    # 3. Traitement
    for rec in records:
        doi = rec.get("LID", "no_doi").split(" ")[0]
        pmid = rec.get("PMID")
        
        citation_count = get_citations_from_doi(doi)
        
        full_data = {
            "SOURCE": "PubMed",
            "PMID": pmid,
            "DOI": doi,
            "TITLE": rec.get("TI"),
            "DATE": rec.get("DP"),
            "JOURNAL": rec.get("JT"),
            "AUTHORS": ", ".join(rec.get("AU", [])),
            "KEYWORDS_EN": ", ".join(rec.get("OT", []) + rec.get("MH", [])),
            "ABSTRACT": rec.get("AB", "No abstract available"),
            "CITATIONS_COUNT": citation_count
        }

        filename = f"master_pmid_{pmid}.txt"
        file_path = os.path.join(output_dir, filename)
        
        with open(file_path, "w", encoding="utf-8") as f:
            for key, val in full_data.items():
                f.write(f"[{key}]\n{val}\n\n")
        
        # --- MISE À JOUR DU COMPTEUR ---
        files_downloaded += 1
        print(f"[{files_downloaded}/{total_to_fetch}] Article {pmid} téléchargé avec succès.")
        time.sleep(0.2)

    # --- AFFICHAGE FINAL ---
    print("-" * 30)
    print(f"RÉSUMÉ : {files_downloaded} fichiers ont été téléchargés et sauvegardés dans '{output_dir}'.")
    print("-" * 30)

# Lancement
fetch_pubmed_with_counter("Covid", limit=500)