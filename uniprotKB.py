import requests
from NCBI import GeneSymbol 

def extraire_info_uniprot(filename = "GeneSymbols_45.txt"):

    dico_espece_gene = GeneSymbol(filename)
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
                        'species' : espece,
                        'uniprot_id': accession,
                        'uniprot_links': f'<a href="https://www.uniprot.org/uniprot/{accession}" target="_blank">{accession}</a>',
                        'protein_name': protein_name,
                        'pdb_id': pdb_entries,
                        'pdb_links': []
                    }
                for pdb_entry in pdb_entries :
                    dico_uniprot_entry["pdb_links"].append(f'<a href="https://www.rcsb.org/structure/{pdb_entry}" target="_blank">{pdb_entry}</a>' )
                

                print(f"Ajouté: {dico_uniprot_entry}") 
                uniprot.append(dico_uniprot_entry)
    return uniprot

# Pour tester
fichier = input("Entrez un fichier avec pour chaque ligne gene,espece : ")
print(extraire_info_uniprot(fichier))
