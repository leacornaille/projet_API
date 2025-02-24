import requests, json

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



def QuickGO (uniprot_id) :

    dico_aspect={}
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
                    #lien GO_term
                    go_link = f"https://amigo.geneontology.org/amigo/term/{goID}"
                    go= f"name: {go_name}, link: {go_link}"
                if aspect not in dico_aspect:
                    dico_aspect[aspect]={}   # Créer un nouveau sous-dictionnaire si l'aspect n'existe pas encore
                dico_aspect[aspect][goID] = set()
                dico_aspect[aspect][goID].add(go)
                                     
    return  dico_aspect

def main_GO (dico_species) : 

    dico ={}
    dico_uniprot_id = uniprot_ID(dico_species)

    for species, uniprot in dico_uniprot_id.items() :
        if uniprot : 
            res = QuickGO(uniprot)
            dico[species]=res
            
    return dico

# Pour tester
#fichier = input("Entrez un fichier avec pour chaque ligne gene,espece : ")
#print(main_GO(fichier))
