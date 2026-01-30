import os

# Configuration
source_dir = "AllBiomodels.xml"
output_file = "liste_modeles_infectieux.txt"

# Mots-clés basés sur votre document (Hôte-pathogène, immunité, etc.)
keywords = ["virus", "bacteria", "infection",]

print(f"Analyse des fichiers dans {source_dir}...")

infectious_count = 0

with open(output_file, "w") as out:
    for filename in os.listdir(source_dir):
        if filename.endswith(".xml"):
            path = os.path.join(source_dir, filename)
            try:
                with open(path, 'r', encoding='utf-8') as f:
                    content = f.read().lower()
                    # Vérification si un mot-clé est présent
                    if any(word in content for word in keywords):
                        out.write(f"{filename}\n")
                        infectious_count += 1
            except:
                continue

print(f"Terminé ! {infectious_count} modèles infectieux détectés.")
print(f"La liste a été sauvegardée dans : {output_file}")
