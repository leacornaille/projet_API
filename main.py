# !/usr/bin/ python3
#-*- coding : utf-8 -*-
import sys
import NCBI
import GO_term
import uniprotKB
import ensembl


## HTML Body
body_html ="" # Initialisation du corps du HTML


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
                {html_transcript}
            </td>
            <td>
                {html_prot}
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