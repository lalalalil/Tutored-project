#!/bin/env python

import sys

if len(sys.argv) < 3 : 
	print("Bienvenue sur convertiseur_GFF. Pour utiliser le programme, veuillez saisir votre fichier en format LST et le nom de votre fichier de sortie GFF.")
	sys.exit(1)

fichier=sys.argv[1]
f_sortie=sys.argv[2]
methode=0 # initialisation des variables
numero_cds=1
numero_version=0
dernier_fin=0
dernier_debut=1

with open(fichier, "r") as fichier_entree, open(f_sortie, "w") as fich_sortie : # ouverture des fichiers 
    for ligne in fichier_entree : # découpage des lignes et suppression des espaces
        ligne=ligne.strip()
        colonnes=ligne.split()

        if methode==0 : # récupération de la méthode de prédiction
            methode=colonnes[0]
            
        if not colonnes : # si colonne pas vide
            continue
        
        if methode=="GENEMARK" :
            if colonnes[0]=="Sequence:" : # récupération du nom de la séquence
                nom_seq=colonnes[1]

            elif len(colonnes) == 7 and colonnes[0].isdigit(): # si on arrive à une colonne de prédiction de CDS
                debut=int(colonnes[0]) # récupération de la position de début
                fin=int(colonnes[1]) # récupération de la position de fin
                score=float(colonnes[5]) # récupération du score
                frame=colonnes[4]
                
                if colonnes[2]=="direct" : # si brin direct
                    brin="+"
                    if fin!=dernier_fin : # si direct alors on regarde la fin de la cds précèdente
                        numero_cds+=1
                        numero_version=1
                    else : 
                        numero_version+=1
                else : # brin complementaire
                    brin="-"
                    if debut!=dernier_debut :  # si complémentaire alors on regarde le début de la cds précédente
                        numero_cds+=1
                        numero_version=1
                    else : 
                        numero_version+=1
                dernier_fin=fin # mise à jour debut et fin
                dernier_debut=debut
                note=f"gene GM_CDS_{numero_cds}.{numero_version}" # formatage des résultats
                sortie=f"{nom_seq}\t{methode}\tCDS\t{debut}\t{fin}\t{score}\t{brin}\t{frame}\t{note}\n"
                fich_sortie.write(sortie)
                if brin=='+' and debut!=1 : # formatage de l'ATG
                    sortie=f"{nom_seq}\t{methode}\tATG\t{debut}\t{debut+2}\t{colonnes[6]}\t{brin}\t.\n"
                    fich_sortie.write(sortie)
                else :
                    sortie=f"{nom_seq}\t{methode}\tATG\t{fin-2}\t{fin}\t{colonnes[6]}\t{brin}\t.\n"
                    fich_sortie.write(sortie)
            
        elif methode=='GeneMark.hmm':
            if colonnes[0]=="FASTA" : # récupération du nom de la séquence
                nom_seq=colonnes[3]

            elif len(colonnes) == 6 and colonnes[0].isdigit(): # si on arrive à une colonne de prédiction de CDS
                debut=int(colonnes[2]) # récupération de la position de début
                fin=int(colonnes[3]) # récupération de la position de fin
                
                sortie=f"{nom_seq}\t{methode}\tCDS\t{debut}\t{fin}\t.\t{colonnes[1]}\t.\tgene GMhmm_CDS_{colonnes[0]}\n"
                fich_sortie.write(sortie)
                if colonnes[1]=='+' and debut!=1 : # formatage de l'ATG
                    sortie=f"{nom_seq}\t{methode}\tATG\t{debut}\t{debut+2}\t.\t{colonnes[1]}\t.\n"
                    fich_sortie.write(sortie)
                else :
                    sortie=f"{nom_seq}\t{methode}\tATG\t{fin-2}\t{fin}\t.\t{colonnes[1]}\t{colonnes[1]}\t.\n"                
                    fich_sortie.write(sortie)
        else :
            print("ERREUR : type de fichier non accepté")

