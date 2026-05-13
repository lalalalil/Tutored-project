import os
import time
from scholarly import scholarly

OUTPUT_DIR = "articles_scholar"

def fetch_scholar_data(disease_name, limit=5):
    print(f"--- Recherche Google Scholar pour : {disease_name} ---")
    
    path = os.path.join(OUTPUT_DIR, disease_name.lower())
    if not os.path.exists(path):
        os.makedirs(path)

    # Requête : allintitle force la recherche des mots dans le titre uniquement
    search_query = scholarly.search_pubs(f'allintitle: "{disease_name}" "model"')

    for i in range(limit):
        try:
            # Récupération du résultat suivant
            result = next(search_query)
            bib = result['bib']
            pmid_fallback = result.get('pub_url', 'no_url').split('/')[-1] # ID de secours
            
            # On prépare le texte à sauvegarder
            content = f"TITRE: {bib.get('title')}\n"
            content += f"ANNEE: {bib.get('pub_year', 'N/A')}\n"
            content += f"CITATIONS: {result.get('num_citations', 0)}\n"
            content += f"URL: {result.get('pub_url', 'N/A')}\n"
            content += f"RESUME: {bib.get('abstract', 'N/A')}\n"

            file_name = f"scholar_{i}.txt"
            with open(os.path.join(path, file_name), "w", encoding="utf-8") as f:
                f.write(content)
            
            print(f"  > Résultat Scholar {i+1} sauvegardé. (Citations: {result.get('num_citations', 0)})")
            
            # IMPORTANT : Délai long pour éviter le bannissement IP par Google
            time.sleep(5) 

        except StopIteration:
            print("Fin des résultats disponibles.")
            break
        except Exception as e:
            print(f"Erreur Scholar : {e}")
            break

# Exemple d'utilisation
fetch_scholar_data("Covid", limit=3)