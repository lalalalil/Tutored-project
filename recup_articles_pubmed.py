import os
import time
from Bio import Entrez

# Configuration
Entrez.email = "ton.nom@etudiant.univ.fr"  # À remplacer par ton mail
OUTPUT_DIR = "articles_pubmed"

def fetch_infectio_data(disease_name, limit=20):
    # 1. Création du dossier de stockage
    path = os.path.join(OUTPUT_DIR, disease_name.lower())
    if not os.path.exists(path):
        os.makedirs(path)

    # 2. Requête optimisée : Titre contient Maladie ET (Model ou variantes) + Date
    # La recherche est insensible à la casse (Maj/Min)
    query = f'({disease_name}[Title] AND model*[Title]) AND 2015:2026[DP]'
    
    print(f"--- Recherche PubMed pour : {disease_name} ---")
    
    try:
        # Recherche des IDs
        search_handle = Entrez.esearch(db="pubmed", term=query, retmax=limit)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        id_list = search_results["IdList"]
        print(f"Nombre d'articles trouvés : {len(id_list)}")

        # 3. Récupération individuelle et sauvegarde
        for pmid in id_list:
            # Petit délai pour respecter les serveurs du NCBI (3 requêtes/sec max)
            time.sleep(0.3) 
            
            fetch_handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
            article_metadata = fetch_handle.read()
            fetch_handle.close()

            # Nom du fichier : article_PMID.txt
            file_name = f"article_{pmid}.txt"
            with open(os.path.join(path, file_name), "w", encoding="utf-8") as f:
                f.write(article_metadata)
            
            print(f"  > Article {pmid} sauvegardé dans {path}")

    except Exception as e:
        print(f"Erreur lors de la recherche pour {disease_name}: {e}")

# Exemple d'utilisation pour ton projet
fetch_infectio_data("Covid", limit=10)
# fetch_infectio_data("Ebola", limit=10)