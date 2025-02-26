# !/usr/bin/ python3
#-*- coding : utf-8 -*-
import sys
import ensembl
import GO_term
import uniprotKB
import NCBI
import ucsc

def html_table(file, mail) :
    ## HTML Body
    body_html ="" # Initialisation du corps du HTML

    # Boucle à travers le fichier d'origine pour créer les lignes du tableau
    with open(file, "r") as infos:
        for line in infos:        
            current_species = line[:-1].split(",")[1]
            current_species_ss_gca = current_species.split("_gca")[0]
            current_species_ss = current_species_ss_gca.replace("_"," ").capitalize()
            symbol = line.split(",")[0]

            # Dictionnaire initial à passer en argument aux sous-scripts
            species_info = {}
            species_info_ncbi_uniprot = {}
            species_info_ucsc = {}
            species_info["species"] = current_species
            species_info["gene_symbol"] = symbol
            species_info_ncbi_uniprot["species"] = current_species_ss_gca
            species_info_ncbi_uniprot["gene_symbol"] = symbol
            species_info_ucsc["species"] = current_species_ss
            species_info_ucsc["gene_symbol"] = symbol
            print(species_info_ucsc)

        # Appel des fonctions principales des sous-scripts
            embl_info = ensembl.InfoGene(species_info)
            go_info = GO_term.main_GO(species_info)
            uniprot_info = uniprotKB.extraire_info_uniprot(species_info_ncbi_uniprot)
            ncbi_info = NCBI.extract_info(species_info_ncbi_uniprot["gene_symbol"],species_info_ncbi_uniprot["species"], mail)
            ucsc_link = ucsc.ucsc_link(species_info_ucsc["gene_symbol"],species_info_ucsc["species"])

            link = f"https://{embl_info['division']}.ensembl.org/{embl_info['species'].capitalize()}" # Raccourci
            
        ## Liens Ensembl pour protéines et transcrits
            # Initialisation : Premiers lien Ensembl
            transcript = embl_info["transcript_id"][0]
            embl_transcript = f"<a href = {link}/Transcript/Summary?db=core;g={embl_info["species"]};t={transcript}>{transcript}</a>"
            if embl_info["prot_id"][0] == "Not translated": # En cas de non traduction
                embl_prot = "<a>Not translated</a>"
            else :
                if embl_info['division'] == 'bacteria':
                    embl_prot = f"<a href = {link}/Transcript/ProteinSummary_{transcript}?db=core;g={embl_info["species"]};t={transcript}>{embl_info["prot_id"][0]}</a>"
                else:
                    embl_prot = f"<a href = {link}/Transcript/ProteinSummary?db=core;g={embl_info["species"]};t={transcript}>{embl_info["prot_id"][0]}</a>"


            # Suite des liens Ensembl
            for i in range(1, len(embl_info["transcript_id"])):
                transcript = embl_info["transcript_id"][i]
                embl_transcript += f"<br>\n\t\t\t\t\t<a href = {link}/Transcript/Summary?db=core;g={embl_info['species']};t={transcript}>{transcript}</a>"

                # Vérifie si l'index 'i' est valide pour 'prot_id'
                if i < len(embl_info["prot_id"]):
                    if embl_info["prot_id"][i] == "Not translated":  # En cas de non traduction
                        embl_prot += "<br>\n\t\t\t\t\t<a>Not translated</a>"
                    else:
                        if embl_info['division'] == 'bacteria':
                            embl_prot += f"<br>\n\t\t\t\t\t<a href = {link}/Transcript/ProteinSummary_{transcript}?db=core;g={embl_info['species']};t={transcript}>{embl_info['prot_id'][i]}</a>"
                        else:
                            embl_prot += f"<br>\n\t\t\t\t\t<a href = {link}/Transcript/ProteinSummary?db=core;g={embl_info['species']};t={transcript}>{embl_info['prot_id'][i]}</a>"
                else:
                    # Si 'prot_id' n'a pas d'élément pour cet 'i', afficher un message indiquant que ce n'est pas traduit
                    embl_prot += "<br>\n\t\t\t\t\t<a>Not translated</a>"
        ## Liens GO term
            # Faire les liens si un UniprotID a été trouvé
            if go_info["UniprotID"] != None :
                mf = go_info["GO"]["molecular_function"]
                bp = go_info["GO"]["biological_process"]
                cc = go_info["GO"]["cellular_component"]

                for key, value in mf.items() :
                    mf_link = f"<br>\n\t\t\t\t\t<a href = https://amigo.geneontology.org/amigo/term/{key}>{key}</a><p> : {value}</p>"

                for key, value in bp.items() :
                    bp_link = f"<br>\n\t\t\t\t\t<a href = https://amigo.geneontology.org/amigo/term/{key}>{key}</a><p> : {value}</p>"

                for key, value in cc.items() :
                    cc_link = f"<br>\n\t\t\t\t\t<a href = https://amigo.geneontology.org/amigo/term/{key}>{key}</a><p> : {value}</p>"

            # En cas d'absence d'UniprotID
            else :
                mf_link = "<br>\n\t\t\t\t\t<p>Data not found</p>"
                bp_link = "<br>\n\t\t\t\t\t<p>Data not found</p>"
                cc_link = "<br>\n\t\t\t\t\t<p>Data not found</p>"

            body_html += f"""
                    <tr>                    
                        <td><div class=header_1>
                            {embl_info['species']}
                        </div></td>
                        <td>
                            {embl_info["gene_symbol"]}
                        </td>
                        <td><a href = {embl_info["gene_browser"]}>
                            View {embl_info["gene_symbol"]} in gene browser
                        </a></td>
                        <td><a href = {link}/Gene/Summary?db=core;g={embl_info['gene_id']}>
                            {embl_info["gene_id"]}
                        </a></td>
                        <td>
                            {embl_transcript}
                        </td>
                        <td>
                            {embl_prot}
                        </td>
                        <td>{"".join(ncbi_info['Links gene'])}
                            
                        </td>
                        <td>{"".join(ncbi_info['Links RNA'])}
                            
                        </td>
                        <td>{"".join(ncbi_info['Links prot'])}
                        </td>
                        <td>
                            {"".join(uniprot_info["uniprot_links"])}
                        </td>
                        <td>
                            {uniprot_info['protein_name']}
                        </td>
                        <td>
                            {"".join(uniprot_info["pdb_links"])}
                        </td>
                        <td>
                            {bp_link}
                        </td>
                        <td>
                            {mf_link}
                        </td>
                        <td>
                            {cc_link}
                        </td>
                        <td>{ucsc_link["lien_ucsc"]}
                        </td>
                    </tr>"""

    ## HTML Head
    head_html = """
    <!DOCTYPE HTML>
    <html lang="fr-FR">
        <head>
            <title> Projet API </title>
            <meta http-equiv="Content-type" content="text/html; charset=utf-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">

            <!-- CSS -->
            <link type="text/css" href="https://cdn.datatables.net/2.2.1/css/dataTables.dataTables.min.css" rel="stylesheet">
            <link type="text/css" href="https://cdn.datatables.net/buttons/3.2.0/css/buttons.dataTables.min.css" rel="stylesheet">
            <link type="text/css" href="https://cdn.datatables.net/colreorder/2.0.4/css/colReorder.dataTables.min.css" rel="stylesheet">
            <link type="text/css" href="https://cdn.datatables.net/fixedheader/4.0.1/css/fixedHeader.dataTables.min.css" rel="stylesheet">
            <link type="text/css" href="https://cdn.datatables.net/fixedcolumns/5.0.4/css/fixedColumns.dataTables.min.css" rel="stylesheet">
            <link type="text/css" href="https://cdn.datatables.net/responsive/3.0.3/css/responsive.dataTables.min.css" rel="stylesheet">
            <link type="text/css" href="https://cdn.datatables.net/rowreorder/1.5.0/css/rowReorder.dataTables.min.css" rel="stylesheet">

            <!-- JavaScript -->
            <script type="text/javascript" language="javascript" src="https://code.jquery.com/jquery-3.7.0.js"></script> 
            <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/2.2.1/js/dataTables.min.js"></script>
            <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/buttons/3.2.0/js/dataTables.buttons.min.js"></script>
            <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/buttons/3.2.0/js/buttons.colVis.min.js"></script>
            <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/colreorder/2.0.4/js/dataTables.colReorder.min.js"></script>
            <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/fixedcolumns/5.0.4/js/dataTables.fixedColumns.min.js"></script>
            <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/fixedheader/4.0.1/js/dataTables.fixedHeader.min.js"></script>
            <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/responsive/3.0.3/js/dataTables.responsive.min.js"></script>
            <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/rowreorder/1.5.0/js/dataTables.rowReorder.min.js"></script>

            <!-- Script Supplémentaire-->
            <script type="text/javascript" class="init">
                $(document).ready(function() {
                    const table = $('#gene').DataTable({
                        <!--Configure le placement des éléments autour du tableau-->
                        layout: { 
                            topStart : 'pageLength',
                            top1Start: {
                                buttons: [{
                                        extend: 'colvis',
                                        columns: 'th:nth-child(n+2)'
                                    }]
                            },
                            bottomStart: 'pageLength',
                            bottom1Start: {
                                buttons: [{
                                        extend: 'colvis',
                                        columns: 'th:nth-child(n+2)'
                                    }]
                            },
                            topEnd : "paging",
                            top1End : "search",
                            bottomEnd : "paging",
                            bottom1End : "search",
                        },
                        fixedHeader: true,
                        "scrollY": true,
                        "scrollX": true,
                        lengthMenu: [5, 10, 25, 50, 75, 100],
                        colReorder: {
                            columns: ':gt(1)'
                        },                                       
                    });
                });
                </script>
        </head>
        <body>
            <h1>Resultat du scripting d'aggrégation automatique des annotations</h1>

            <table id="gene" class="display nowrap cell-border" style="width:100%">
                <thead>
                    <tr>
                        <th><div class=header_1>Species</div></th>
                        <th>Gene Symbol</th>
                        <th>Gene Browser</th>
                        <th>Gene - EMBL ID</th>
                        <th>Transcript - EMBL ID</th>
                        <th>Protein - EMBL ID</th>
                        <th>Gene - NCBI ID</th>
                        <th>Transcript - NCBI ID</th>
                        <th>Protein - NCBI ID</th>
                        <th>Uniprot ID</th>
                        <th>Protein Name</th>
                        <th>3D structure - PDB</th>
                        <th>Biological process - GO term</th>
                        <th>Molecular function - GO term</th>
                        <th>Cellular component - GO term</th>
                        <th> Lien UCSC Genome Browser</th>
                    </tr>
                </thead>
                <tbody>
    """

    ## HTML Tail
    tail_html = f""" 
                </tbody>
            </table>
        </body>
    </html>
    """

    html = open("Results.html","w")
    html.write(head_html[1:]+body_html+tail_html)
    html.close()

if __name__ == '__main__':
    filename = sys.argv[1]
    mail = sys.argv[2]
    html_table(filename, mail)