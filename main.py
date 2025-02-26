# !/usr/bin/ python3
#-*- coding : utf-8 -*-
import sys
import ensembl
import GO_term
import uniprotKB
import NCBI
import ucsc
import time

start = time.time()

def html_table(file, mail, filename_output) :
    """
Point d'entrée principal du script.

Premier argument : Nom du fichier contenant les gènes et les espèces (exemple : GeneSymbols_45.txt).
Deuxième argument : Email pour faire les requêtes à l'API NCBI (exemple : exemple@mail.com).
Troisième argument : Nom du fichier de sortie (exemple : Results.html).
    """
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
            species_info["species"] = current_species
            species_info["gene_symbol"] = symbol

            species_info_ncbi_uniprot = {}
            species_info_ncbi_uniprot["species"] = current_species_ss_gca
            species_info_ncbi_uniprot["gene_symbol"] = symbol

            species_info_ucsc = {}
            species_info_ucsc["species"] = current_species_ss
            species_info_ucsc["gene_symbol"] = symbol
            print(species_info_ucsc)

        # Appel des fonctions principales des sous-scripts
            embl_info = ensembl.InfoGene(species_info)
            go_info = GO_term.main_GO(species_info_ncbi_uniprot)
            uniprot_info = uniprotKB.extraire_info_uniprot(species_info_ncbi_uniprot)
            ncbi_info = NCBI.extract_info(species_info_ncbi_uniprot["gene_symbol"],species_info_ncbi_uniprot["species"], mail)
            ucsc_link = ucsc.ucsc_link(species_info_ucsc["gene_symbol"],species_info_ucsc["species"])

            link = f"https://{embl_info['division']}.ensembl.org/{embl_info['species'].capitalize()}"  # Raccourci
            
            if embl_info['division'] == 'vertebrates':
                link = f"https://www.ensembl.org/{embl_info['species'].capitalize()}"

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
                try :
                    mf = go_info["GO"]["molecular_function"]
                    bp = go_info["GO"]["biological_process"]
                    cc = go_info["GO"]["cellular_component"]

                    # Initialisation des variables de liens
                    mf_link = ""
                    bp_link = ""
                    cc_link =""

                    # Création et concaténation des liens
                    for key, value in mf.items() :
                        mf_link += f"<br>\n\t\t\t\t\t<a href = https://amigo.geneontology.org/amigo/term/{key}>{key}</a> : {value}"

                    for key, value in bp.items() :
                        bp_link += f"<br>\n\t\t\t\t\t<a href = https://amigo.geneontology.org/amigo/term/{key}>{key}</a> : {value}"

                    for key, value in cc.items() :
                        cc_link += f"<br>\n\t\t\t\t\t<a href = https://amigo.geneontology.org/amigo/term/{key}>{key}</a> : {value}"

            # En cas d'absence d'UniprotID
                except :
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

                        <td>
                            {embl_info['division']}
                        </td>

                        <td>
                            {ncbi_info['Official name']}
                        </td>

                        <td><a href = {link}/Gene/Summary?db=core;g={embl_info['gene_id']}>
                            {embl_info["gene_id"]}
                        </a></td>

                        <td>
                            {"".join(ncbi_info['Links gene'])}                            
                        </td>

                        <td><a href = {embl_info["gene_browser"]}>
                            View {embl_info["gene_symbol"]} in Ensembl genome browser
                        </a><br>
                            {ucsc_link["lien_ucsc"]}
                        </td>

                        <td><div class='scroll'>
                            {embl_transcript}
                        </div></td>

                        <td><div class='scroll'>
                            {"".join(ncbi_info['Links RNA'])}                            
                        </div></td>
                        
                        <td>
                            {uniprot_info['protein_name']}
                        </td>

                        <td><div class='scroll'>
                            {embl_prot}
                        </div></td>

                        <td><div class='scroll'>
                            {"".join(ncbi_info['Links prot'])}
                        </div></td>

                        <td><div class='scroll'>
                            {"".join(uniprot_info["uniprot_links"])}
                        </div></td>

                        <td><div class='scroll'>
                            {"".join(uniprot_info["pdb_links"])}
                        </div></td>

                        <td><div class='scroll'>
                            {bp_link}
                        </div></td>

                        <td><div class='scroll'>
                            {mf_link}
                        </div></td>

                        <td><div class='scroll'>
                            {cc_link}
                        </div></td>
                        
                    </tr>"""

    ## HTML Head
    head_html = """
    <!DOCTYPE HTML>
    <html lang="fr-FR">
        <head>
            <title>Results : Gene Annotation Automatic Table</title>
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
            <script type="text/javascript" language="javascript" src="https://cdnjs.cloudflare.com/ajax/libs/jszip/3.10.1/jszip.min.js"></script>
            <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/buttons/3.2.0/js/buttons.html5.min.js"></script>
            <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/colreorder/2.0.4/js/dataTables.colReorder.min.js"></script>
            <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/fixedcolumns/5.0.4/js/dataTables.fixedColumns.min.js"></script>
            <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/fixedheader/4.0.1/js/dataTables.fixedHeader.min.js"></script>
            <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/responsive/3.0.3/js/dataTables.responsive.min.js"></script>
            <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/rowreorder/1.5.0/js/dataTables.rowReorder.min.js"></script>

            <!-- Script Supplémentaire-->
            <script type="text/javascript" class="init">
                $(document).ready(function() {
                    $('#gene').DataTable({
                        <!--Configure le placement des éléments autour du tableau-->
                        layout: { 
                            topStart : 'pageLength',
                            top1Start: {
                                buttons: [
                                    {
                                        extend: 'colvis',
                                        columns: 'th:nth-child(n+2)'
                                    },
                                    {
                                        extend: 'collection',
                                        text: 'Export',
                                        buttons: ['csv','excel']
                                    }
                                ]
                            },
                            bottomStart: 'pageLength',
                            bottom1Start: {
                                buttons: [
                                    {
                                        extend: 'colvis',
                                        columns: 'th:nth-child(n+2)'
                                    },
                                    {
                                        extend: 'collection',
                                        text: 'Export',
                                        buttons: ['csv','excel']
                                    }
                                ]
                            },
                            topEnd : "paging",
                            top1End : "search",
                            bottomEnd : "paging",
                            bottom1End : "search",
                        },
                        fixedHeader: true,
                        scrollY: true,
                        scrollX: true,
                        lengthMenu: [5, 10, 25, 50, 75, 100],
                        fixedColumns: {
                          leftColumns: 2
                        },
                        colReorder: {
                            columns: ':gt(1)'
                        },                                       
                    });
                });
            </script>

            <style>

            h1{
                padding-top: 40px;
                padding-bottom: 40px;
                text-align: center;
            }

            .header{
                background-color: #E8E8E8;
                margin: 0 ;
                height: 120px;
                margin-bottom: 10px;
                border-bottom: 1px solid #ddd;
            }

            a{
                color: black;
            }
            
            .scroll {white-space:nowrap;
                max-height: 120px;
                max-width: 230px;
                overflow: auto;
	        }
            
            td{
                text-align: justify;
            }

            </style>
        </head>

        <body>
            <div class="header">
            <h1>Automatic Gene Annotation Table</h1>
            </div>
            <table id="gene" class="display nowrap cell-border" style="width:100%">
                <thead>
                    <tr>
                        <th><div class=header_1>Organism</div></th>
                        <th>Gene Symbol</th>
                        <th>Division</th>
                        <th>Official Full Name</th>
                        <th>Gene Access Number - ENSEMBL</th>
                        <th>Gene Access Number - NCBI</th>
                        <th>Genome Browser</th>
                        <th>RNA Access Numbers - ENSEMBL</th>
                        <th>RNA Access Numbers - NCBI</th>
                        <th>Protein Name</th>
                        <th>Protein Access Numbers - ENSEMBL</th>
                        <th>Protein Access Numbers - NCBI</th>
                        <th>Protein Access Numbers - UniProt</th>
                        <th>3D structure - PDB</th>
                        <th>Biological process - GO term</th>
                        <th>Molecular function - GO term</th>
                        <th>Cellular component - GO term</th>
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

    html = open(filename_output,"w")
    html.write(head_html[1:]+body_html+tail_html)
    html.close()

if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) < 3:
        print("Usage: python3 main.py <fichier_genes> <email> <nom_fichier_sortie>")
        sys.exit(1)
    filename = args[0]
    mail = args[1]
    filename_output = args[2]
    html_table(filename, mail, filename_output)
    print(f"Fichier HTML {filename_output} créée")
    print(time.ctime(time.time() - start)[11:19])