# !/usr/bin/ python3
#-*- coding : utf-8 -*-
import requests
import sys
import os

def extraire_info_uniprot(dico_espece_gene):
    dico_uniprot_entry = {  # Initialisation vide pour éviter l'erreur
        'species': dico_espece_gene["species"],
        'gene_symbol': dico_espece_gene["gene_symbol"],
        'uniprot_id': 'N/A',
        'uniprot_links': "<br>\n\t\t\t\t\t<p>Data not found</p>",
        'protein_name': 'N/A',
        'pdb_id': [],
        'pdb_links': []
    }

    pdb_entries = []
    espece = dico_espece_gene["species"]
    gene = dico_espece_gene["gene_symbol"]
    query = f"gene:{gene} AND reviewed:true AND {espece}"
    params = {
        "query": query,
        "fields": "accession,protein_name,xref_pdb",
        "sort": "accession desc",
        "size": 1
    }
    headers = {
        "accept": "application/json"
    }
    base_url = "https://rest.uniprot.org/uniprotkb/search"

    response = requests.get(base_url, headers=headers, params=params)
    data = response.json()
    results = data.get("results", [])

    if not results:
        query = f"gene:{gene} AND {espece}"  # Suppression de "reviewed:true"
        params["query"] = query
        response = requests.get(base_url, headers=headers, params=params)
        data = response.json()
        results = data.get("results", [])

    if results:
        entry = results[0]  # Prend le premier résultat
        accession = entry.get('primaryAccession', 'N/A')
        protein_desc = entry.get('proteinDescription', {}).get('recommendedName', {})
        protein_name = protein_desc.get('fullName', {}).get('value', 'N/A')

        for xref in entry.get("uniProtKBCrossReferences", []):
            if xref.get("database") == "PDB":
                pdb_entries.append(xref.get("id"))

        # Mise à jour du dictionnaire
        dico_uniprot_entry.update({
            'uniprot_id': accession,
            'uniprot_links': f'<br>\n\t\t\t\t\t<a href="https://www.uniprot.org/uniprot/{accession}" target="_blank">{accession}</a>',
            'protein_name': protein_name,
            'pdb_id': pdb_entries,
            'pdb_links': [f'<br>\n\t\t\t\t\t<a href="https://www.rcsb.org/structure/{pdb}" target="_blank">{pdb}</a>' for pdb in pdb_entries]
        })

    return dico_uniprot_entry


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