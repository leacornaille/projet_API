import Ensembl as ens
import NCBI
import GO_term
import uniprotKB
import sys
from collections import defaultdict

### Fonctions
def read(fileName_input):
  """
Lit un fichier d'entrée et extrait un dictionnaire d'espèces et de gènes.

Args: 
  fileName_input (str) : Nom du fichier d'entrée. Chaque ligne doit être au format "gene,espece".
    
Returns:
  dict: Un dictionnaire avec les espèces comme clés et leurs gènes associés comme valeurs.    
  """
  dico = {}
  try:
    with open(fileName_input, "r") as fichier :
      for line in fichier :
        li = line[:-1].split(",")
        espece = li[1]
        gene = li[0]
        if espece in dico:
          dico[espece] += gene
        else :
          dico[espece] = gene
  except FileNotFoundError:
    print(f"Fichier introuvable : {fileName_input}")
  return dico


def main():
  """
Point d'entrée principal du script.

Premier argument : Nom du fichier contenant les gènes et les espèces (exemple : GeneSymbols_45.txt).
Deuxième argument : Email pour faire les requêtes à l'API NCBI (exemple : exemple@mail.com).
  """
  args = sys.argv[1:]

  if len(args) < 3:
    print("Usage: python3 main.py <fichier_genes> <email> <nom_fichier_sortie>")
    sys.exit(1)

  if not dico:
    print("Aucun gène ou espèce valide trouvé dans le fichier d'entrée.")
    sys.exit(1)
  else:
    print("Données extraites")

  ### Récupération des résultats
  
  #### Ensembl
  dico = read(args[0])
  try:
    ensembl_results = ens.ensembl_extract(dico)
    print("Résultats Ensembl récupérés")
  except Exception as e:
    print(f"Erreur lors de l'extraction Ensembl: {e}")

  #### NCBI
  try:
    ncbi_results = NCBI.lien_NCBI(args[0])
    print("Résultats NCBI récupérés")
  except Exception as e:
    print(f"Erreur lors de l'extraction NCBI: {e}")

  #### GO Term
  try:
    GO_results = GO_term.main_GO(args[0])
    print("Résultats NCBI récupérés")
  except Exception as e:
    print(f"Erreur lors de l'extraction NCBI: {e}")

  #### uniprot
  try:
    uniprot_results = uniprotKB.extraire_info_uniprot(args[0])
    print("Résultats NCBI récupérés")
  except Exception as e:
    print(f"Erreur lors de l'extraction NCBI: {e}")


  ## Création du corps et de la fin du fichier HTML
  html_content =""

  # Ajout des données dans le tableau HTML
  print("Création du fichier html en cours...")
  html_content += f"""
                  <tr>
                      <td>{ensembl_results['gene_symbol']}</td>
                      <td>{ensembl_results['species']}</td>
                      <td>{ensembl_results['gene_name']}</td>
                      <td>
                          <a href="{result['gene_access_number_link']}">{result['gene_access_number']}</a>
                      </td>
                      <td>
                          <a href="{result['gene_link']}">{result['gene_id']}</a>
                      </td>
                      <td>
                          <a href="{result['genome_browser_link']}">Ensembl Genome Browser</a>
                      </td>
                      <td>
      """
      
      # Ajout des RNA Access Numbers
      if result['rna_access_number_link'] is not None :
        if len(result['rna_access_number_link']) > 1:
          for rna, link in zip(result['rna_access_number'], result['rna_access_number_link']):
            html_content += f'<a href="{link}">{rna}</a><br>'
        else:        
          html_content += f'<a href="{result["rna_access_number_link"][0]}">{result["rna_access_number"][0]}</a>'
      else:
        html_content += f'{result["rna_access_number"]}'

      html_content += "</td><td>"

      if result['rna_link'] is not None :
        if len(result['rna_link']) > 1:
            for rna, link in zip(result['rna_id'], result['rna_link']):
                html_content += f'<a href="{link}">{rna}</a><br>'
        else:        
            html_content += f'<a href="{result["rna_link"][0]}">{result["rna_id"][0]}</a>'
      else:
        html_content += f'{result["rna_id"]}'

      html_content += "</td><td>"

      if result['rna_pred_link'] is not None :
        if len(result['rna_pred_link']) > 1:
            for rna, link in zip(result['rna_pred_id'], result['rna_pred_link']):
                html_content += f'<a href="{link}">{rna}</a><br>'
        else:        
            html_content += f'<a href="{result["rna_pred_link"][0]}">{result["rna_pred_id"][0]}</a>'
      else:
         html_content += f'{result["rna_pred_id"]}'

      html_content += "</td><td>"

      # Ajout des Protein Access Numbers
      if result['protein_access_number_link'] is not None :
        if len(result['protein_access_number_link']) > 1:
            for prot, link in zip(result['protein_access_number'], result['protein_access_number_link']):
                html_content += f'<a href="{link}">{prot}</a><br>'
        else:        
            html_content += f'<a href="{result["protein_access_number_link"][0]}">{result["protein_access_number"][0]}</a>'
      else:
        html_content += f'{result["protein_access_number"]}'

      html_content += "</td><td>"

      if result['prot_link'] is not None :
        if len(result['prot_link']) > 1:
            for prot, link in zip(result['prot_id'], result['prot_link']):
                html_content += f'<a href="{link}">{prot}</a><br>'
        else:        
            html_content += f'<a href="{result["prot_link"][0]}">{result["prot_id"][0]}</a>'
      else:
        html_content += f'{result["prot_id"]}'

      html_content += "</td><td>"

      if result['prot_pred_link'] is not None:
        if len(result['prot_pred_link']) > 1:
            for prot, link in zip(result['prot_pred_id'], result['prot_pred_link']):
                html_content += f'<a href="{link}">{prot}</a><br>'
        else:        
            html_content += f'<a href="{result["prot_pred_link"][0]}">{result["prot_pred_id"][0]}</a>'
      else:
        html_content += f'{result["prot_pred_id"]}'

      html_content += "</td></tr>"
  print("Création du fichier html en cours.......")

  # Fermeture du tableau et du fichier HTML
  html_content += """
              </tbody>
          </table>

          <script>
              $(document).ready(function () {
                  $('#geneTable').DataTable({
                      responsive: true,
                      scrollX: true,
                      scrollY: '800px',
                      fixedHeader: true,
                      colReorder: true,
                      fixedColumns: {
                          leftColumns: 2
                      }
                  });
              });
          </script>
      </body>
  </html>
  """

  ## Création du corps et de la fin du fichier HTML
  html_head = """
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
          
      </style>
      
      </head>
      <body>
          <div class="header">
          <h1>Gene Annotation Automatic Table</h1>
          </div>
          <table id="geneTable" class="display" style="width:100%">
              <thead>
                  <tr>
                      <th>Gene Symbol</th>
                      <th>Organism</th>
                      <th>Official full name</th>
                      <th>Ensembl Gene Access Number</th>
                      <th>NCBI Gene Access Number</th>
                      <th>Genome Browser</th>
                      <th>Ensembl RNA Access Numbers</th>
                      <th>NCBI RNA Access Numbers</th>
                      <th>Predicted RNA Access Number(s)</th>
                      <th>Ensembl Protein Access Numbers</th>
                      <th>NCBI Protein Access Numbers</th>
                      <th>Predicted Protein Access Number(s)</th>
                  </tr>
              </thead>
              <tbody>
  """

  print("Création du fichier html en cours..........")

  # Écriture dans le fichier HTML
  fileName_output = args[2]
  with open(fileName_output, "w") as html_file:
      html_file.write(html_head + html_content)

  print(f"Fichier HTML généré : {fileName_output}")

if __name__ == "__main__":
  main()