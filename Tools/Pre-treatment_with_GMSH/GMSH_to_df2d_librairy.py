#
# Librairie principale de conversion GMSH -> DASSFLOW2D
#
# convert_gmsh_to_df2d : fonction principale de conversion/génaation des fichiers
#
# /!\ ATTENTION : les fichiers maillages .msh doivent respecter certaines conventions de nommage
#              pour les groupes physiques (INLET_*, OUTLET_*, WALL_*) notamment.
#       --> Voir tuto (Tuto_mesh_conventions.py) associé pour plus de détails.
#

import gmsh, os
from convertion_tools import process_bc, generate_hydrograph_file, generate_hpresc_file, generate_zpresc_file, generate_ratcurve_file, generate_land_uses_file, gen_input


def convert_gmsh_to_df2d(Input_MSH, Sim_duration = 3600, Default_Q=10.0, Default_H=1.0, Default_Z=1.0, Output_BC="bc.txt", Output_LAND="land_uses.txt", Land_uses_config={1: {"friction": 0.033, "beta": 0.0} }, folder_output="./"):
    if folder_output != "./":
        os.makedirs(folder_output, exist_ok=True)
    Output_BC = os.path.join(folder_output, Output_BC)
    Output_LAND =  os.path.join(folder_output, Output_LAND)
    
    Output_GEO = Input_MSH.replace(".msh", ".geo")
    Output_GEO_NAME = Output_GEO.split("/")[-1]

    Output_GEO = os.path.join(folder_output, Output_GEO_NAME)
    
    gmsh.initialize()
    try:
        print(f"--- 1. Chargement du maillage {Input_MSH} ---")
        gmsh.open(Input_MSH)
        # --- A. NOEUDS ---
        nodeTags, nodeCoords, _ = gmsh.model.mesh.getNodes()
        node_map = {tag: i+1 for i, tag in enumerate(nodeTags)}
        nodes_list = [(nodeCoords[i*3], nodeCoords[i*3+1]) for i in range(len(nodeTags))]

        node_z_coords = [nodeCoords[i*3+2] for i in range(len(nodeTags))]
        # --- B. CELLULES & TOPOLOGIE ---
        all_cells = [] 
        edge_map = {} 

        cells_bathy = [] 

        for dim in [2, 3]:
            tags, nodes = gmsh.model.mesh.getElementsByType(dim)
            n_p = 3 if dim == 2 else 4
            for i in range(len(tags)):
                cid = len(all_cells) + 1 
                el = [node_map[nodes[i*n_p + j]] for j in range(n_p)]

                # --- CALCUL BATHYMETRIE (Moyenne Z) ---
                z_sum = 0.0
                for nid in el:
                    z_sum += node_z_coords[nid-1]
                cells_bathy.append(z_sum / n_p)

                if dim == 2: el.append(el[0]) # Tri -> Quad
                all_cells.append(el)
                for k in range(n_p):
                    edge_map[tuple(sorted((el[k], el[(k+1)%n_p])))] = (cid, k+1)

        # --- C. CONDITIONS LIMITES (MULTI-FILES) ---
        inlet_lines, outlet_lines = [], []
        bc_entries = [] # (ID, Type, File)
        files_data = {} # Map: 'hydrograph.txt' -> [{'id': 1, 'variant': 1}, ...]

        phys_groups = gmsh.model.getPhysicalGroups(1)
        print(f"--- 2. Analyse des conditions aux limites ---")

        # 1. Séparation
        inlet_tags, outlet_tags = [], []
        for p_dim, p_tag in phys_groups:
            p_name = gmsh.model.getPhysicalName(p_dim, p_tag).upper()
            if "WALL" in p_name: continue
            if "INLET" in p_name: inlet_tags.append((p_tag, p_name))
            elif "OUTLET" in p_name: outlet_tags.append((p_tag, p_name))


        # 2. Traitement
        inlet_lines = process_bc(inlet_tags, True, 0, bc_entries, files_data, node_map, edge_map)
        outlet_lines = process_bc(outlet_tags, False, len(inlet_tags), bc_entries, files_data, node_map, edge_map)
        

        # --- D. ECRITURE GEO ---
        print(f"--- 3. Ecriture de {Output_GEO} ---")
        with open(Output_GEO, 'w') as f:
            f.write(f"$ Generated\n{len(nodes_list)} {len(all_cells)} 0\n$ Nodes\n")
            for i, c in enumerate(nodes_list): f.write(f"{i+1} {c[0]:.6f} {c[1]:.6f} 0.0\n")
            f.write("$ Cells\n")
            for i, n in enumerate(all_cells): 
                z_val = cells_bathy[i]
                f.write(f"{i+1} {n[0]} {n[1]} {n[2]} {n[3]} 1 {z_val:.6f}\n")
            
            f.write(f"$ INLET\nINLET {len(inlet_lines)} {len(inlet_tags)}\n")
            for l in inlet_lines: 
                ghost_z = cells_bathy[l[0]-1]
                f.write(f"{l[0]} {l[1]} 0 {ghost_z:.6f} {l[2]}\n")
            f.write(f"OUTLET {len(outlet_lines)} {len(outlet_tags)}\n")
            for l in outlet_lines: 
                ghost_z = cells_bathy[l[0]-1]
                f.write(f"{l[0]} {l[1]} 0 {ghost_z:.6f} {l[2]}\n")
                
        # --- E. ECRITURE BC.TXT ---
        print(f"--- 4. Ecriture de {Output_BC} ---")
        with open(Output_BC, 'w') as f:
            f.write("!\n!\n!\n" + f"{len(bc_entries)}\n" + "!\n!\n!\n")
            for bid, btype, fname in bc_entries:
                f.write(f"{bid}\t{btype}\t{fname}\n")

        # --- F. GENERATION DATA (INTELLIGENTE) ---
        print(f"--- 5. Génération data (avec variantes) ---")
        generate_land_uses_file(Output_LAND, Land_uses_config)

        for fname, entries in files_data.items():
            if fname == "hydrograph.txt":
                generate_hydrograph_file(fname, entries, Sim_duration=Sim_duration, Default_Q=Default_Q, folder_output=folder_output)
            elif fname == "hpresc.txt":
                #generate_hpresc_file(fname, Sim_duration=Sim_duration, Default_H=Default_H, folder_output=folder_output)
                generate_hpresc_file(fname, entries, Sim_duration=Sim_duration, Default_H=Default_H, folder_output=folder_output)
            elif fname == "zpresc.txt":
                #generate_zpresc_file(fname, Sim_duration=Sim_duration, Default_Z=Default_Z, folder_output=folder_output)
                generate_zpresc_file(fname, entries, Sim_duration=Sim_duration, Default_Z=Default_Z, folder_output=folder_output)
            elif fname == "rating_curve.txt":
                generate_ratcurve_file(fname, entries, folder_output=folder_output)
            else:
                with open(fname, 'w') as f: f.write("$ Placeholder\n")

        print("--- TERMINE AVEC SUCCES ---")

    except Exception as e:
        print(f"ERREUR : {e}")
    finally:
        gmsh.finalize()

    input_params = { "mesh_name":Output_GEO_NAME,
                                            "ts":Sim_duration, 
                                            "dta":100, 
                                            "dtw":10, 
                                            "dtp":10,
                                            "dt":1, 
                                            "temp_scheme":'euler', 
                                            "spatial_scheme":'first_b1', 
                                            "adapt_dt":1, 
                                            "cfl":0.8, 
                                            "feedback_inflow":1,                                    
                                            "coef_feedback":0.8,
                                            "heps":0, 
                                            "friction":1, 
                                            "g":10, 
                                            "w_tecplot":1,
                                            "w_vtk":1, 
                                            "w_gnuplot":1, 
                                            "w_obs":0, 
                                            "use_obs":0, 
                                            "max_nt_for_adjoint":2500,
                                            "c_manning":0,
                                            "c_manning_beta":0,
                                            "c_bathy":0,
                                            "c_hydrograph":0,
                                            "c_ratcurve":0,
                                            "c_rain":0,
    #										"c_infil":0,
                                            "c_ic":0,
                                            "restart_min":0,
                                            "eps_min":0.0001}

    gen_input(input_params, folder_output=folder_output)


if __name__ == "__main__":
    Input_mesh = "YOUR_MESH.msh"
    print(f"=== CONVERSION Test for {Input_mesh} ===")
    convert_gmsh_to_df2d(Input_MSH=Input_mesh, Sim_duration=3600, Default_Q=10.0, Default_H=1.0, Default_Z=1.0)