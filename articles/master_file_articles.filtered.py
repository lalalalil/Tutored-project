import os
import requests
import time
from Bio import Entrez, Medline

# --- CONFIGURATION ---
Entrez.email = "ton.nom@etudiant.univ.fr" 

def get_citations_from_doi(doi):
    """Récupère le nombre de citations via Crossref."""
    if not doi or "no_doi" in doi:
        return 0
    try:
        url = f"https://api.crossref.org/works/{doi}"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            # On retourne un entier (int) pour pouvoir comparer
            return int(data['message'].get('is-referenced-by-count', 0))
    except:
        pass
    return 0

def fetch_filtered_pubmed_data(disease, limit=500, min_citations=5):
    output_dir = f"data_{disease.lower()}_filtered"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # REQUÊTE STRICTE : La maladie ET le mot "model" doivent être dans le TITRE
    query = f"({disease}[Title] AND model*[Title]) AND English[Language] AND 2015:2026[DP]"
    
    print(f"--- Recherche PubMed (Filtre: >{min_citations} citations) : {disease} ---")
    
    search_handle = Entrez.esearch(db="pubmed", term=query, retmax=limit)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    
    ids = search_results["IdList"]
    if not ids:
        print("Aucun article trouvé avec ces mots-clés dans le titre.")
        return

    fetch_handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
    records = list(Medline.parse(fetch_handle))
    fetch_handle.close()

    files_saved = 0
    ignored_low_impact = 0
    total_found = len(records)

    for i, rec in enumerate(records, 1):
        doi = rec.get("LID", "no_doi").split(" ")[0]
        pmid = rec.get("PMID")
        
        # 1. RÉCUPÉRATION DU NOMBRE DE CITATIONS
        citation_count = get_citations_from_doi(doi)
        
        # 2. LE FILTRE (Condition : Citations > 5)
        if citation_count > min_citations:
            full_data = {
                "SOURCE": "PubMed",
                "PMID": pmid,
                "DOI": doi,
                "TITLE": rec.get("TI"),
                "DATE": rec.get("DP"),
                "JOURNAL": rec.get("JT"),
                "AUTHORS": ", ".join(rec.get("AU", [])),
                "KEYWORDS": ", ".join(rec.get("OT", []) + rec.get("MH", [])),
                "ABSTRACT": rec.get("AB", "No abstract available"),
                "CITATIONS_COUNT": str(citation_count)
            }

            filename = f"master_pmid_{pmid}.txt"
            with open(os.path.join(output_dir, filename), "w", encoding="utf-8") as f:
                for key, val in full_data.items():
                    f.write(f"[{key}]\n{val}\n\n")
            
            files_saved += 1
            print(f"[{i}/{total_found}] ✅ Article {pmid} sauvegardé ({citation_count} citations).")
        else:
            ignored_low_impact += 1
            print(f"[{i}/{total_found}] ❌ Article {pmid} ignoré (Seulement {citation_count} citations).")
        
        time.sleep(0.2) # Courtoisie API

    print("\n" + "="*40)
    print(f"RÉSUMÉ POUR : {disease.upper()}")
    print(f"- Articles analysés : {total_found}")
    print(f"- Articles sauvegardés (> {min_citations} citations) : {files_saved}")
    print(f"- Articles rejetés (<= {min_citations} citations) : {ignored_low_impact}")
    print("="*40)

# Lancement
fetch_filtered_pubmed_data("Covid", limit=500, min_citations=5)
