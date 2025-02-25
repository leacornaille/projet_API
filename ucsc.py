import requests

# Dictionnaire pour stocker les génomes
def ucsc_link(gene,species_name):
    genome_browser = []

# URL de l'API UCSC pour récupérer la liste des génomes de référence
    url = 'https://api.genome.ucsc.edu/list/ucscGenomes'

# Envoi de la requête GET pour récupérer les données
    response = requests.get(url)

# Vérification si la requête a réussi
    if response.status_code == 200:
    # Récupérer les données JSON
        data = response.json()

    # Accéder au dictionnaire des génomes
        genomes = data.get('ucscGenomes', {})

        # Filtrer les génomes correspondant à l'espèce cible
        filtered_genomes = {genome: info for genome, info in genomes.items() if info.get('scientificName') == species_name}

        if filtered_genomes:
        # Ne conserver que le dernier génome trouvé pour l'espèce
            last_genome = list(filtered_genomes.keys())[-1]  # Dernier génome dans la liste triée
            genome_browser = {"gene_symbol": gene, "genome": last_genome, "species": species_name, 
                                "lien_ucsc":f'</br>\n\t\t\t\t\t<a href="https://genome.ucsc.edu/cgi-bin/hgTracks?db={last_genome}&position={gene}">{gene}</a>' }
        else:
            genome_browser = {"gene_symbol": gene, "genome": "Aucun genome trouvé", "species": species_name,"lien_ucsc": "Non présent sur l'UCSC"}
    return genome_browser

'''
    with open(species_file, 'r') as file:
        species_lines = file.readlines()

    for line in species_lines:
        # Extraction du nom de l'espèce et du symbole du gène
        parts = line.strip().split(',')
        if len(parts) < 2:
            continue  # Ignorer les lignes mal formatées

        gene, species_name = parts[0], parts[1]

        # Correction du format des espèces
        species_name = species_name.replace('_', ' ')
        species_parts = species_name.split(' ', 1)
        species_name = ' '.join([species_parts[0].capitalize()] + [species_parts[1].lower()]) if len(species_parts) > 1 else species_parts[0].capitalize()
'''