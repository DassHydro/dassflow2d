import os
import sys

def process_results_for_hbanks(res_directory="simus/bin_reference_copy/figures_total/res_without_poro", output_filename="porosity_hbanks.csv"):
    """
    Parcourt tous les fichiers de résultats .dat dans un dossier,
    trouve la hauteur d'eau maximale (h) pour chaque cellule sur toute la simulation,
    et écrit le résultat dans un fichier de sortie.
    """
    # Dictionnaire pour stocker la hauteur max pour chaque ID de cellule
    # Format: {cell_id: max_h}
    hbanks_per_cell = {}

    print(f"Analyse des résultats dans le dossier '{res_directory}'...")

    # 1. Vérifier que le dossier de résultats existe
    if not os.path.isdir(res_directory):
        print(f"ERREUR : Le dossier '{res_directory}' est introuvable.")
        sys.exit(1)

    # 2. Lister tous les fichiers de résultats
    result_files = [f for f in os.listdir(res_directory) if f.startswith("result_") and f.endswith(".dat")]
    
    if not result_files:
        print("Aucun fichier de résultat trouvé. Assurez-vous que la simulation de référence a bien été exécutée.")
        return

    print(f"{len(result_files)} fichiers trouvés. Début du traitement...")

    # 3. Parcourir chaque fichier et chaque ligne pour trouver le max
    for filename in result_files:
        filepath = os.path.join(res_directory, filename)
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    if line.strip().startswith("#"):
                        continue
                    
                    parts = line.strip().split()
                    
                    if len(parts) >= 5: # On a besoin au minimum des 5 premières colonnes
                        cell_id = int(parts[0])
                        h_value = float(parts[5]) # La hauteur 'h' est la 5ème colonne (indice 4)
                        
                        # Mettre à jour la valeur maximale pour cette cellule
                        current_max = hbanks_per_cell.get(cell_id, 0.0)
                        if h_value > current_max:
                            hbanks_per_cell[cell_id] = h_value
                            
        except (ValueError, IndexError) as e:
            print(f"  Avertissement : Ligne ignorée dans {filename} à cause d'une erreur de formatage ({e})")
            continue
    
    print("Traitement terminé.")

    # 4. Écrire les résultats dans le fichier de sortie
    if not hbanks_per_cell:
        print("Aucune donnée valide n'a pu être lue. Le fichier de sortie n'a pas été créé.")
        return

    with open(output_filename, 'w') as f:
        f.write("# i Hbanks \n")
        # Trier par numéro de cellule pour un fichier propre
        for cell_id in sorted(hbanks_per_cell.keys()):
            max_h = hbanks_per_cell[cell_id]
            f.write(f"{cell_id};{max_h} \n")
            
    print(f"Fichier de sortie '{output_filename}' créé avec les données de {len(hbanks_per_cell)} cellules.")
    
    with open("porosity_beta.csv", 'w') as f:
        f.write("# i beta \n")
        for cell_id in sorted(hbanks_per_cell.keys()):
            f.write(f"{cell_id};2 \n")

# --- Point d'entrée du script ---
if __name__ == "__main__":
    process_results_for_hbanks()
