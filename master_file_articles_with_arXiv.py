import os
import arxiv
import requests
from Bio import Entrez, Medline
import time

Entrez.email = "ton.nom@etudiant.univ.fr" 

def get_citations_from_doi(doi):
    """Interroge Crossref pour obtenir le nombre de citations via le DOI."""
    if not doi or "no_doi" in doi or "arxiv" in doi:
        return "0"
    
    try:
        # API Crossref gratuite
        url = f"https://api.crossref.org/works/{doi}"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            return str(data['message'].get('is-referenced-by-count', 0))
    except Exception:
        pass
    return "0"

def fetch_exhaustive_data_with_citations(disease, limit=1000):
    output_dir = f"data_{disease.lower()}_final"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    articles_db = {}

    # --- 1. PUBMED ---
    print(f"--- PubMed : {disease} ---")
    query = f"({disease}[Title/Abstract] AND model*[Title]) AND English[Language] AND 2015:2026[DP]"
    search_handle = Entrez.esearch(db="pubmed", term=query, retmax=limit)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    
    ids = search_results["IdList"]
    if ids:
        fetch_handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
        records = list(Medline.parse(fetch_handle))
        for rec in records:
            doi = rec.get("LID", "no_doi").split(" ")[0]
            articles_db[doi] = {
                "SOURCE": "PubMed",
                "PMID": rec.get("PMID"),
                "DOI": doi,
                "TITLE": rec.get("TI"),
                "DATE": rec.get("DP"),
                "ABSTRACT": rec.get("AB", "No abstract"),
                "KEYWORDS": ", ".join(rec.get("OT", []) + rec.get("MH", [])),
                "AUTHORS": ", ".join(rec.get("AU", []))
            }

    # --- 2. arXiv ---
    print(f"--- arXiv : {disease} ---")
    client = arxiv.Client()
    search = arxiv.Search(query = f'ti:"{disease}" AND ti:"model"', max_results = limit)
    for result in client.results(search):
        doi = result.doi if result.doi else f"arxiv_{result.get_short_id()}"
        if doi not in articles_db:
            articles_db[doi] = {
                "SOURCE": "arXiv",
                "DOI": doi,
                "TITLE": result.title,
                "DATE": result.published.strftime("%Y-%m-%d"),
                "ABSTRACT": result.summary,
                "KEYWORDS": ", ".join(result.categories),
                "AUTHORS": ", ".join([a.name for a in result.authors])
            }

    # --- 3. RÉCUPÉRATION DES CITATIONS (Crossref) ---
    print(f"--- Récupération des citations via Crossref ---")
    for doi in articles_db:
        # On ne cherche les citations que si on a un vrai DOI
        if "arxiv" not in doi:
            citations = get_citations_from_doi(doi)
            articles_db[doi]["CITATIONS_COUNT"] = citations
            time.sleep(0.2) # Courtoisie pour l'API
        else:
            articles_db[doi]["CITATIONS_COUNT"] = "0 (Preprint)"

    # --- 4. SAUVEGARDE ---
    for doi, data in articles_db.items():
        filename = "".join([c for c in doi if c.isalnum() or c in ("_", "-")])
        with open(os.path.join(output_dir, f"master_{filename}.txt"), "w", encoding="utf-8") as f:
            for key, val in data.items():
                f.write(f"[{key}]\n{val}\n\n")
            
    print(f"Terminé : {len(articles_db)} articles sauvegardés avec compte de citations.")

fetch_exhaustive_data_with_citations("Covid", limit=500)