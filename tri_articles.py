from Bio import Entrez
import os

Entrez.email = "ton.nom@etudiant.univ.fr"

def search_covid_models(limit=10):
    # Création du répertoire de sortie
    output_dir = "covid_models_data"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # REQUÊTE : (Covid dans le titre) ET (Model dans le titre) ET (Après 2015)
    query = '(Covid[Title] AND Model[Title]) AND 2015:2026[DP]'
    
    print(f"Lancement de la requête : {query}")

    # 1. Recherche des IDs
    handle = Entrez.esearch(db="pubmed", term=query, retmax=limit)
    record = Entrez.read(handle)
    handle.close()

    id_list = record["IdList"]
    print(f"Articles trouvés : {len(id_list)}")

    # 2. Récupération et sauvegarde
    for pmid in id_list:
        fetch_handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
        data = fetch_handle.read()
        fetch_handle.close()

        with open(f"{output_dir}/article_{pmid}.txt", "w", encoding="utf-8") as f:
            f.write(data)
        print(f"PMID {pmid} sauvegardé.")

search_covid_models(limit=10)