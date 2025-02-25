# !/usr/bin/ python3
#-*- coding : utf-8 -*-
import requests
import os
import sys

# Extraire les uniprot_ID
def uniprot_ID (dico_espece) : 

    gene = dico_espece["gene_symbol"]
    species = dico_espece["species"]
      
    url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene}+AND+{species}&format=json&fields=accession"
    response = requests.get(url)

    if response.ok : 
        data= response.json()
        results = data.get("results", [])
        
        if results : 
            dico_espece["UniprotID"] = results[0].get("primaryAccession") 
        else : 
            dico_espece["UniprotID"] = None
        
    else : 
        print(f"Aucune information trouver pour {species} : {gene} ")
 
    return dico_espece



def QuickGO (dico_espece) :

    dico_aspect={}
    uniprot_id = dico_espece["UniprotID"]

    # requete pour avoir le goID
    serveur = "https://www.ebi.ac.uk/QuickGO/services"
    ext1 = f"/annotation/search?geneProductId={uniprot_id}&limit=200&page=1"
    response = requests.get(serveur+ext1, headers={"Accept": "application/json"})
    if response.ok:
        data = response.json()
        results = data.get("results", [])
        for res in results :
            aspect = res.get("goAspect")                
            goID = res.get("goId")

            #requete API pour avoir les GO term 
            ext2 = f"/ontology/go/terms/{goID}"
            response = requests.get(serveur+ext2, headers={"Accept": "application/json"})
            if response.ok:
                data = response.json()
                results = data.get("results", [])
                if results : 
                    go_name=results[0].get("name")
                if aspect not in dico_aspect:
                    dico_aspect[aspect]={}   # Créer un nouveau sous-dictionnaire si l'aspect n'existe pas encore
                dico_aspect[aspect][goID] = go_name

                dico_espece["GO"]= dico_aspect
                                     
    return  dico_espece

def main_GO (dico_espece) : 

    dico_espece = uniprot_ID(dico_espece)
    dico_espece = QuickGO(dico_espece)
            
    return dico_espece

# En cas de lancement par ligne de commande : lit le fichier en argument
if __name__ == '__main__':
    gene_symbols = sys.argv[1]
    if os.path.isfile(gene_symbols): # Vérifie que l'argument soit un fichier
        with open(gene_symbols, "r") as infos:
            for line in infos:
                symbol = line.split(",")[0]
                current_species = line[:-1].split(",")[1]

                # Dictionnaire initial à passer en argument aux sous-scripts
                species_info = {}
                species_info["species"] = current_species
                species_info["gene_symbol"] = symbol

                print("\t",main_GO(species_info))
    else :
        print("Erreur : Veuillez donnez un nom de fichier accessible comme seul argument \n"\
              "Format de chaque ligne : [Symbole de gène],[Espèce] \n" \
              "Ex: RAD51,homo_sapiens")
