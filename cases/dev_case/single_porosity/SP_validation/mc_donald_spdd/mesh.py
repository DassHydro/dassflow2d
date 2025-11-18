import os

def gen_channel_3cells(filename, L, dx, H_top=0.0):
    """
    Génère un canal rectiligne avec 3 cellules par section (3 rangs en largeur).
    
    filename : nom du fichier .geo à créer
    L : longueur totale en mètres
    dx : résolution en x
    H_top : altitude de base (bathymétrie uniforme)
    """

    nx = int(L // dx)        # nombre de colonnes
    ny = 3                   # 3 cellules verticales → 4 noeuds

    dy = 3.0                 # épaisseur verticale arbitraire (1 m)

    nodes = []
    cells = []

    # ==========================
    # 1. GÉNÉRATION DES NŒUDS
    # ==========================

    node_id = 1
    for i in range(nx + 1):
        x = i * dx
        for j in range(ny + 1):
            y = j * dy
            z = H_top
            nodes.append((node_id, x, y, z))
            node_id += 1

    # ==========================
    # 2. GÉNÉRATION DES CELLULES
    # ==========================

    cell_id = 1
    for i in range(nx):
        for j in range(ny):
            # indices des 4 noeuds
            n1 = i*(ny+1) + j + 1
            n2 = i*(ny+1) + j + 2
            n3 = (i+1)*(ny+1) + j + 2
            n4 = (i+1)*(ny+1) + j + 1

            cells.append((cell_id, n1, n2, n3, n4, 1, H_top))
            cell_id += 1

    # ==========================
    # 3. ÉCRITURE DU FICHIER GEO
    # ==========================

    with open(filename, "w") as f:
        f.write("# Nodes\n")
        for n in nodes:
            f.write(f"{n[0]:7d} {n[1]:14.7E} {n[2]:14.7E} {n[3]:14.7E}\n")

        f.write("\n# Cells\n")
        for c in cells:
            f.write(f"{c[0]:7d} {c[1]:7d} {c[2]:7d} {c[3]:7d} {c[4]:7d} {c[5]:4d} {c[6]:14.7E}\n")

    return filename

code_dir =  os.getcwd() #os.path.join(dassflow_dir,"code")
bin_dir = os.path.join(code_dir,"bin_A")
mesh_name = "channel3.geo"
gen_channel_3cells(os.path.join(bin_dir, mesh_name), L=1000, dx=10, H_top=0.0)