# !/usr/bin/ python3
#-*- coding : utf-8 -*-
import requests
import sys
import os

def extraire_info_uniprot(dico_espece_gene):

    uniprot =[]
    pdb_entries = []
    for espece, gene in dico_espece_gene.items():
        print(espece)
        query = f"gene:{gene} AND reviewed:true AND {espece}"
        params = {
        "query": query,
        "fields": "accession,protein_name,xref_pdb",
        "sort": "accession desc",
        "size": 1  # Limite à 3 résultats
    }

        headers = {
        "accept": "application/json"
    }
        base_url = "https://rest.uniprot.org/uniprotkb/search"

        response = requests.get(base_url, headers=headers, params=params)
    
        if not response.ok:
            continue  # Passe au gène suivant au lieu d'arrêter

    # Extraction des données JSON
        data = response.json()
        results = data.get("results", [])

        if not results:
            query = f"gene:{gene} AND {espece}"  # Suppression de "reviewed:true"
            params = {
        "query": query,
        "fields": "accession,protein_name,xref_pdb",
        "sort": "accession desc",
        "size": 1  # Limite à 3 résultats
    }
            response = requests.get(base_url, headers=headers, params=params)
        
            if not response.ok:
                continue
        
            data = response.json()
            results = data.get("results", [])


    # Affichage des résultats pour ce gène

        if not results:
            continue
        else:
            for entry in results:
                accession = entry.get('primaryAccession', 'N/A')
                protein_desc = entry.get('proteinDescription', {}).get('recommendedName', {})
                protein_name = protein_desc.get('fullName', {}).get('value', 'N/A')
                for xref in entry.get("uniProtKBCrossReferences"):
                    if xref.get("database") == "PDB":
                        pdb_entries.append(xref.get("id"))

            # Stockage dans les dictionnaires
                dico_uniprot_entry = {
                        'gene_symbol': gene,
                        'uniprot_id': accession,
                        'uniprot_links': f'<br>\n\t\t\t\t\t<a href="https://www.uniprot.org/uniprot/{accession}" target="_blank">{accession}</a>',
                        'protein_name': protein_name,
                        'pdb_id': pdb_entries,
                        'pdb_links': []
                    }
                for pdb_entry in pdb_entries :
                    dico_uniprot_entry["pdb_links"].append(f'<br>\n\t\t\t\t\t<a href="https://www.rcsb.org/structure/{pdb_entry}" target="_blank">{pdb_entry}</a>' )
                

                print(f"Ajouté: {dico_uniprot_entry}") 
                uniprot.append(dico_uniprot_entry)
    return uniprot

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

                print("\t",extraire_info_uniprot(species_info))
    else :
        print("Erreur : Veuillez donnez un nom de fichier accessible comme seul argument \n"\
              "Format de chaque ligne : [Symbole de gène],[Espèce] \n" \
              "Ex: RAD51,homo_sapiens")