# !/usr/bin/ python3
#-*- coding : utf-8 -*- 

from Bio import Entrez 
import json 


def get_gene_id (gene_name, organism, mail) : 
    Entrez.email = mail

    # extraction des ID pour chaque gene ID 
    handle=Entrez.esearch(db="gene", term = f'{organism} [Orgn] AND {gene_name} [gene]')
    record = Entrez.read(handle)
    ids=record["IdList"]
    gene_id = ids[0]
    handle.close()

    # extraction du nom officiel 
    handle = Entrez.esummary(db="gene", id= gene_id, retmode="json")
    result = json.loads(handle.read())
    handle.close()

    official_name = result["result"].get(gene_id, {}).get("description")

    return ids, official_name


def get_linked_ids (gene_ids, target_db, mail) : 
    Entrez.email = mail

    # Requete pour avoir les IDs dans les autres db à partir du gene ID 
    handle = Entrez.elink(dbfrom = "gene", db= target_db, id = gene_ids, retmode= 'json')
    result = json.loads(handle.read())
    handle.close()

    # extraction des IDs pour les différents db 
    linked_ids = []
    for linkset in result.get("linksets", []) : 
        for linksetdb in linkset.get("linksetdbs", []) :
            linked_ids.extend(linksetdb.get("links", []))
    
    return linked_ids


def get_seq_info (seq_ids, db, mail) : 
    Entrez.email = mail

    sequences = []
    
    # Requete pour avoir les RNA number acces ou transcrit number access selon la db 
    handle = Entrez.esummary(db = db, id = ",".join(seq_ids), retmode = 'json')
    result = json.loads(handle.read())
    handle.close()

    # extraction des informatiosn RNA number access ou transcrit number access
    for seq_id in seq_ids : 
        seq_info = result.get("result", {}).get(str(seq_id), {})
        if seq_info : 
            acc_version = seq_info.get("accessionversion","")  
            if acc_version.startswith(("NP_", "NM_", "XP_", "XM_")) : 
                if acc_version not in sequences : 
                    sequences.append(acc_version)
    
    return sequences


def extract_info (gene_symbol, organism, mail) : 

    gene_id, official_name = get_gene_id(gene_symbol, organism, mail)

    if gene_id :
        # extraction des ids pour les db prot et nucléotide
        protein_ids = get_linked_ids(gene_id, "protein", mail)  # db protein
        nucleotide_ids = get_linked_ids(gene_id, 'nuccore', mail) # db nucleotide

        # extraction des numéros d'access prot et transcrit
        protein_info = get_seq_info(protein_ids, "protein", mail) or ["Data not found"]
        nucleotide_info = get_seq_info(nucleotide_ids, 'nuccore', mail) or ["Data not found"]

        #Creaction des lien
        lien_RNA = [ f'<br>\n\t\t\t\t\t<a href="https://www.ncbi.nlm.nih.gov/nuccore/{rna}" target="_blank">{rna}</a>'
                         for rna in nucleotide_info ] if nucleotide_info[0]!= "Data not found" else ["Data not found"]

        lien_prot = [f'<br>\n\t\t\t\t\t<a href="https://www.ncbi.nlm.nih.gov/protein/{prot}" target="_blank">{prot}</a>' 
                          for prot in protein_info ] if protein_info[0]!= "Data not found" else ["Data not found"]

        lien_gene =  [f'<br>\n\t\t\t\t\t<a href="https://www.ncbi.nlm.nih.gov/gene/{gene_id[0]}" target="_blank">{gene_id[0]}</a>']


    return   {"Species" : organism,
              "Gene Symbol " : gene_symbol,
              "Official name" : official_name,
              "Gene ID" : gene_id, 
              "Protein" : protein_info, 
              "Transcript" : nucleotide_info,
              "Links RNA" : lien_RNA,
              "Links prot" : lien_prot,
              "Links gene" : lien_gene,
              }

