# 🧬 Script d'Annotation Automatique de Gènes

## 📌 Description

Ce script génère un tableau HTML contenant des annotations pour des gènes spécifiés. Il récupère des informations depuis plusieurs bases de données :

- **Ensembl** (gènes, transcrits, protéines)
- **NCBI** (noms officiels, liens)
- **UniProtKB** (protéines, structures 3D, liens UniProt)
- **GO terms** (processus biologiques, fonctions moléculaires, composants cellulaires)
- **UCSC** (liens vers le navigateur UCSC)

## 🚀 Installation

### 1️⃣ Prérequis

- **Python 3** installé sur votre machine
- Connexion Internet pour accéder aux bases de données
- Bibliothèques Python nécessaires (voir ci-dessous)

### 2️⃣ Installer les dépendances

Utilisez `pip` pour installer les bibliothèques nécessaires :

```bash
pip install requests biopython
```

## 🔧 Utilisation

Le script s'exécute en ligne de commande avec les arguments suivants :

```bash
python3 main.py <fichier_genes> <email> <fichier_sortie>
```

### 📌 Arguments

- `fichier_genes` : Nom du fichier contenant les gènes et les espèces (ex: `GeneSymbols.txt`)
- `email` : Adresse email pour les requêtes NCBI
- `fichier_sortie` : Nom du fichier de sortie (ex: `Results.html`)

### 🎯 Exemple d'exécution

```bash
python3 main.py GeneSymbols.txt exemple@mail.com Results.html
```

