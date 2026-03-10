#
# Librairie de fonctions pour la mise à jour des fichiers DASSFLOW2D
#
# Fonctions principales :
#   - apply_bathymetry : applique une nouvelle bathymétrie à un fichier
#   - update_boundary_files : met à jour les fichiers de conditions aux limites (hydrographes, hpresc, etc.)
#

import numpy as np
import os


def parse_geo_file(file_path):
    """
    Fonction de récupération des noeuds et cellules depuis un fichier .geo DASSFLOW2D.
    Retourne deux dictionnaires : nodes et cells.
    nodes : {node_id: (x, y)}
    cells : {cell_id: [node_id1, node_id2, node_id3, (node_id4)]}  # node_id4 optionnel pour triangles
    
    :param file_path: Chemin vers le fichier .geo
    """
    nodes = {}
    cells = {}
    current_section = None
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('$') or line.startswith('['):
                    if "$ Nodes" in line: current_section = "nodes"
                    if "$ Cells" in line: current_section = "cells"
                    continue
                
                parts = line.split()
                
                if current_section == "nodes":
                    try:
                        n_id = int(parts[0])
                        x, y = float(parts[1]), float(parts[2])
                        nodes[n_id] = (x, y)
                    except (ValueError, IndexError): continue
                        
                elif current_section == "cells":
                    try:
                        c_id = int(parts[0])
                        n1, n2, n3, n4 = int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4])                      
                        # DETECTION TRIANGLE : Si le dernier noeud est égal au premier
                        if n4 == n1:
                            # C'est un triangle, on ne garde que les 3 premiers
                            cells[c_id] = [n1, n2, n3]
                        else:
                            # C'est un vrai quad, on garde les 4
                            cells[c_id] = [n1, n2, n3, n4]                           
                    except (ValueError, IndexError): continue                        
        return nodes, cells
    except Exception as e:
        print(f"Erreur : {e}")
        return {}, {}

def calcul_center_cell(nodes, cells, cell_id):
    """
    Calcule le centre (barycentre) d'une cellule donnée.
    :param nodes: Dictionnaire des noeuds {node_id: (x, y)}
    :param cells: Dictionnaire des cellules {cell_id: [node_id1, node_id2, node_id3, (node_id4)]}
    :param cell_id: ID de la cellule à traiter
    :return: (x_center, y_center) du centre de la cellule
    """
    node_indices = cells[cell_id]
    cell_coords = np.array([nodes[nid] for nid in node_indices])
    center = np.mean(cell_coords, axis=0) 
    return center

def calcul_ghost_center(nodes, cells, cell_id, edge_id):
    """
    Calcule le centre de la cellule fantôme associée à l'arête edge_id de la cellule cell_id.
    :param nodes: Dictionnaire des noeuds {node_id: (x, y)}
    :param cells: Dictionnaire des cellules {cell_id: [node_id1, node_id2, node_id3, (node_id4)]}
    :param cell_id: ID de la cellule à traiter
    :param edge_id: ID de l'arête (1 à 4) de la cellule
    :return: (x_ghost, y_ghost) du centre fantôme de la cellule
    """
    if cell_id not in cells:
        print(f"Cellule {cell_id} introuvable.")
        return None
        
    c_nodes_ids = cells[cell_id]
    num_nodes = len(c_nodes_ids) # Sera 3 pour triangle, 4 pour quad
    
    # --- VERIFICATION ---
    # On ne peut pas demander l'arête 4 d'un triangle !
    if edge_id > num_nodes:
        print(f"Erreur : L'arête {edge_id} n'existe pas pour une cellule à {num_nodes} noeuds.")
        return None

    # --- 1. BARYCENTRE PHYSIQUE ---
    # Grâce au parser nettoyé, la moyenne est maintenant correcte (divisé par 3 ou 4 selon le cas)
    sum_x = sum(nodes[nid][0] for nid in c_nodes_ids if nid in nodes)
    sum_y = sum(nodes[nid][1] for nid in c_nodes_ids if nid in nodes)
    valid_pts = sum(1 for nid in c_nodes_ids if nid in nodes)
    
    if valid_pts != num_nodes:
        print("Erreur : Certains noeuds sont manquants dans la liste des noeuds.")
        return None
        
    c_phy_x = sum_x / valid_pts
    c_phy_y = sum_y / valid_pts
    
    # --- 2. SELECTION DES NOEUDS DE L'ARETE ---
    # Logique universelle : L'arête i connecte le noeud (i-1) au noeud (i)
    # Si on arrive au bout de la liste, on boucle sur le premier (modulo)
    
    idx_a = edge_id - 1
    idx_b = (edge_id) % num_nodes # Si edge=4 et num=4 -> index 0. Si edge=3 et num=3 -> index 0.
    
    n_a = c_nodes_ids[idx_a]
    n_b = c_nodes_ids[idx_b]
    
    # --- 3. CALCUL GEOMETRIQUE ---
    xa, ya = nodes[n_a]
    xb, yb = nodes[n_b]
    
    m_edge_x = (xa + xb) / 2.0
    m_edge_y = (ya + yb) / 2.0
    
    c_ghost_x = 2 * m_edge_x - c_phy_x
    c_ghost_y = 2 * m_edge_y - c_phy_y
    
    return (c_ghost_x, c_ghost_y)

# ------------------------------------------------------------------------------------------------ #
#                   FONCTIONS PRINCIPALES DE MODIFICATION DE FICHIERS DASSFLOW2D
# ------------------------------------------------------------------------------------------------ #

def apply_bathymetry(f, geo_file_path):
    """
    Parcourt le fichier .geo, calcule la bathymétrie pour chaque cellule (interne et ghost)
    en utilisant la fonction f(x, y) et modifie le fichier sur place.
    
    Args:
        f (function): Une fonction qui prend (x, y) et retourne une valeur float (la profondeur).
        geo_file_path (str): Le chemin vers le fichier .geo.
    """
    
    # --- ETAPE 1 : Chargement de la géométrie ---
    print(f"Chargement du maillage {geo_file_path}...")
    nodes, cells = parse_geo_file(geo_file_path)
    
    if not nodes or not cells:
        print("Erreur critique : Impossible de parser le fichier ou fichier vide.")
        return

    # --- ETAPE 2 : Traitement et Réécriture ---
    new_lines = []
    current_section = None
    
    print("Application de la bathymétrie...")
    
    try:
        with open(geo_file_path, 'r', encoding='utf-8') as file:
            for line in file:
                original_line = line.strip()
                
                # Si la ligne est vide, on garde telle quelle
                if not original_line:
                    new_lines.append(line)
                    continue
                
                # --- DETECTION DES SECTIONS ---
                if "$ Nodes" in line:
                    current_section = "nodes"
                    new_lines.append(line)
                    continue
                elif "$ Cells" in line:
                    current_section = "cells"
                    new_lines.append(line)
                    continue
                elif "INLET" in line or "OUTLET" in line:
                    # Ces lignes marquent souvent le début d'un bloc de conditions aux limites
                    current_section = "boundary"
                    new_lines.append(line)
                    continue
                elif line.startswith('$'):
                    # Autres sections avec $ (ex: $ EndNodes)
                    current_section = None
                    new_lines.append(line)
                    continue

                parts = original_line.split()
                
                # Si la ligne ne commence pas par un chiffre, c'est probablement du texte/header, on ignore
                if not parts[0].isdigit():
                    new_lines.append(line)
                    continue

                # --- TRAITEMENT SELON LA SECTION ---
                
                # CAS 1 : C'est une cellule physique ($ Cells)
                if current_section == "cells":
                    try:
                        c_id = int(parts[0])
                        
                        # 1. Calcul du centre physique
                        center = calcul_center_cell(nodes, cells, c_id)
                        
                        # 2. Application de la fonction f(x,y)
                        new_bathy = f(center[0], center[1])
                        
                        # 3. Modification de la ligne (Dernière colonne = bathy)
                        parts[-1] = f"{new_bathy:.6f}"
                        
                        # On reconstruit la ligne avec les espacements originaux ou tabulations
                        new_lines.append("\t".join(parts) + "\n")
                        
                    except Exception as e:
                        print(f"Erreur calcul cell {parts[0]}: {e}")
                        new_lines.append(line) # En cas d'erreur, on garde l'originale

                # CAS 2 : C'est une condition limite (INLET/OUTLET)
                elif current_section == "boundary":
                    try:
                        # Format attendu : CellID EdgeID Flag Bathy Flag
                        # Ex: 472 4 0 0.095058 1
                        
                        c_id = int(parts[0])
                        edge_id = int(parts[1])
                        
                        # 1. Calcul du centre fantôme
                        ghost_center = calcul_ghost_center(nodes, cells, c_id, edge_id)
                        
                        if ghost_center is None:
                            # Si le calcul échoue (ex: arête invalide), on garde la ligne
                            new_lines.append(line)
                            continue

                        # 2. Application de la fonction f(x,y) sur le fantôme
                        new_ghost_bathy = f(ghost_center[0], ghost_center[1])
                        
                        # 3. Modification de la ligne (4ème colonne, index 3)
                        parts[3] = f"{new_ghost_bathy:.6f}"
                        
                        new_lines.append("\t".join(parts) + "\n")
                        
                    except Exception as e:
                        print(f"Erreur calcul boundary cell {parts[0]}: {e}")
                        new_lines.append(line)

                # CAS 3 : Autres (Nodes, etc.), on ne touche pas
                else:
                    new_lines.append(line)

        # --- ETAPE 3 : Ecriture finale ---
        with open(geo_file_path, 'w', encoding='utf-8') as f_out:
            f_out.writelines(new_lines)
            
        print("Fichier mis à jour avec succès.")

    except Exception as e:
        print(f"Erreur globale lors du traitement : {e}")




def update_boundary_files(file_path, target_group_id=1, times=None, values=None):
    """
    Met à jour les fichiers de conditions aux limites (hydrograph.txt, hpresc.txt, etc.).
    Fonctionne pour hydrographe, hpresc, zpresc. Ratcurve à tester.
    
    Args:
        file_path (str): Chemin complet du fichier.
        target_group_id (int): L'index du bloc à modifier (1 = 1er bloc, 2 = 2ème bloc, etc.). (si plusieurs hydrographes par exemple)
        times (list): Liste des temps.
        values (list): Liste des valeurs.
    """
    if times is None or values is None or len(times) != len(values):
        print(f"Erreur : Les listes de temps ({len(times) if times is not None else 'None'}) et de valeurs ({len(values) if values is not None else 'None'}) doivent avoir la même taille.")
        return

    if not os.path.exists(file_path):
        print(f"Erreur : Le fichier cible n'existe pas : {file_path}")
        return

    target_idx = target_group_id - 1
    if target_idx < 0:
        print("Erreur : L'index de position doit être supérieur ou égal à 1.")
        return

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            
        iterator = iter(lines)
        blocks = []
        global_header = []

        try:
            global_header.append(next(iterator))
            global_header.append(next(iterator))
            global_header.append(next(iterator))
            num_blocks_line = next(iterator)
            total_blocks = int(num_blocks_line.strip())
        except StopIteration:
            print(f"Erreur : Le fichier {file_path} est vide ou mal formaté.")
            return

        # --- 3. Lecture de tous les blocs ---
        for i in range(total_blocks):
            block = {}
            try:
                # 3 lignes d'en-tête de bloc (on les garde telles quelles)
                block['header'] = [next(iterator), next(iterator), next(iterator)]
                
                # Nombre de points du bloc
                count_line = next(iterator)
                num_points = int(count_line.strip())
                
                # Données du bloc
                block['data_lines'] = []
                for _ in range(num_points):
                    block['data_lines'].append(next(iterator))
                
                blocks.append(block)
            except StopIteration:
                print(f"Attention : Fin de fichier inattendue au bloc n°{i+1}.")
                break

        if target_idx >= len(blocks):
            print(f"Erreur : Vous demandez le bloc n°{target_group_id}, mais le fichier n'en contient que {len(blocks)}.")
            return

        print(f" -> Mise à jour du bloc n°{target_group_id} (sur {len(blocks)} blocs présents)...")

        new_data_lines = []
        for t, v in zip(times, values):
            new_data_lines.append(f"{float(t)} \t {float(v)}\n")
            
        # Remplacement
        blocks[target_idx]['data_lines'] = new_data_lines
        
        # --- 5. Réécriture du fichier ---
        with open(file_path, 'w', encoding='utf-8') as f:
            for line in global_header:
                f.write(line)
            f.write(f"{len(blocks)}\n") # Nombre total
            
            for block in blocks:
                for h_line in block['header']:
                    f.write(h_line)
                f.write(f"{len(block['data_lines'])}\n")
                for d_line in block['data_lines']:
                    f.write(d_line)
                    
        print(f"Succès : Le fichier {file_path} a été mis à jour.")

    except Exception as e:
        print(f"Erreur critique lors du traitement : {e}")

if __name__ == "__main__":
    print("=== Update file Tools librairy ===")