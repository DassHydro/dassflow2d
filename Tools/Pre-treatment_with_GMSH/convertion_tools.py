# 
# Librairie d'outils pour la conversion GMSH -> DASSFLOW2D, utilisée par GMSH_to_df2d_librairy.py
#
# Pas vraiment utilisable en dehors de ce contexte hormis gen_input()
#


import gmsh,os

# =================================================

def parse_bc_name(name):
    """
    Découpe 'INLET:discharg1_2' en :
    - type: 'discharg1' (pour bc.txt)
    - variant: 2 (pour générer des données différentes)
    """
    # Nettoyage initial
    clean = name.upper().replace('"', '').strip()
    
    # Extraction de la partie après ':'
    raw_type = clean.split(":")[1].lower() if ":" in clean else "unknown"
    
    # Gestion du cas par défaut si pas de ':'
    if raw_type == "unknown":
        if "INLET" in clean: raw_type = "discharg1"
        elif "OUTLET" in clean: raw_type = "hpresc"
        
    # Extraction du suffixe _X
    if "_" in raw_type:
        base_type, suffix = raw_type.split("_")
        try:
            variant = int(suffix)
        except ValueError:
            variant = 1
    else:
        base_type = raw_type
        variant = 1
        
    return base_type, variant

# Fonction helper pour traiter une CL
def process_bc(tag_list, is_inlet,start_id_offset=0, bc_entries = [], files_data={}, node_map={}, edge_map={}):
    BC_FILE_MAP = {"discharg1": "hydrograph.txt","discharg2": "hydrograph.txt","hpresc":    "hpresc.txt","zpresc":    "zpresc.txt","ratcurve":  "rating_curve.txt"}
 
    processed_lines = []
    
    for idx, (p_tag, p_name) in enumerate(tag_list):
        # ID Logic
        local_id = idx + 1
        global_id = start_id_offset + local_id
        
        # Parsing Intelligent (Type + Variante)
        base_type, variant = parse_bc_name(p_name)
        fname = BC_FILE_MAP.get(base_type, "unknown.txt")
        
        # Ajout BC.TXT (On écrit le type de base ex: 'discharg1')
        bc_entries.append((global_id, base_type, fname))
        
        # Stockage Données (On garde la variante pour le générateur)
        if fname not in files_data: files_data[fname] = []
        files_data[fname].append({'id': global_id, 'variant': variant})
        
        # Géométrie
        for entity in gmsh.model.getEntitiesForPhysicalGroup(1, p_tag):
            l_tags, l_nodes = gmsh.model.mesh.getElementsByType(1, tag=entity)
            for k in range(len(l_tags)):
                key = tuple(sorted((node_map[l_nodes[2*k]], node_map[l_nodes[2*k+1]])))
                if key in edge_map:
                    # .GEO utilise le LOCAL ID
                    processed_lines.append([edge_map[key][0], edge_map[key][1], local_id])
    return processed_lines

# =================================================

def generate_hydrograph_file(filename, entries, Sim_duration=3600, Default_Q=10.0, folder_output="./"):
    """
    entries: liste de dictionnaires {'id': GlobalID, 'variant': int}
    """
    filename = os.path.join(folder_output, filename)
    print(f" -> Création de {filename} avec {len(entries)} courbe(s)...")
    with open(filename, 'w') as f:
        f.write("$ Hydrograph file generated automatically\n$\n$\n")
        f.write(f"{len(entries)}\n") # Nombre total d'hydrogrammes
        
        for i, entry in enumerate(entries):
            gid = entry['id']
            var = entry['variant']
            f.write(f"$ Hydrograph {i+1} for Group ID {gid} (Variant {var})\n$\n$\n")
            f.write("2\n") # Nombre de points
            f.write(f"0.0     {Default_Q}\n")
            f.write(f"{Sim_duration} {Default_Q}\n")

"""
def generate_hpresc_file(filename, Sim_duration=3600, Default_H=1.0, folder_output="./"):
    filename = os.path.join(folder_output, filename)
    print(f" -> Création de {filename}...")
    with open(filename, 'w') as f:
        f.write(f"$ {filename} file generated automatically\n$\n$\n")
        f.write("1\n$\n$\n$\n2\n")
        f.write(f"0.0     {Default_H}\n")
        f.write(f"{Sim_duration} {Default_H}\n")
"""
def generate_hpresc_file(filename, entries, Sim_duration=3600, Default_H=1.0, folder_output="./"):
    """
    Génère un fichier de hauteurs imposées (hpresc.txt) supportant plusieurs courbes.
    entries: liste de dictionnaires {'id': GlobalID, 'variant': int}
    """
    filename = os.path.join(folder_output, filename)
    print(f" -> Création de {filename} avec {len(entries)} courbe(s)...")
    with open(filename, 'w') as f:
        # En-tête global
        f.write("$ Hpresc file generated automatically\n$\n$\n")
        f.write(f"{len(entries)}\n") # Nombre total de courbes hpresc
        
        # Boucle sur chaque entrée
        for i, entry in enumerate(entries):
            gid = entry['id']
            var = entry['variant']
            f.write(f"$ Hpresc {i+1} for Group ID {gid} (Variant {var})\n$\n$\n")
            f.write("2\n") # Nombre de points (ici constant : t_start et t_end)
            f.write(f"0.0     {Default_H}\n")
            f.write(f"{Sim_duration} {Default_H}\n")
"""
def generate_zpresc_file(filename, Sim_duration=3600, Default_Z=1.0, folder_output="./"):
    filename = os.path.join(folder_output, filename)
    print(f" -> Création de {filename}...")
    with open(filename, 'w') as f:
        f.write(f"$ {filename} file generated automatically\n$\n$\n")
        f.write("1\n$\n$\n$\n2\n")
        f.write(f"0.0     {Default_Z}\n")
        f.write(f"{Sim_duration} {Default_Z}\n")
"""
def generate_zpresc_file(filename, entries, Sim_duration=3600, Default_Z=1.0, folder_output="./"):
    """
    Génère un fichier de cotes imposées (zpresc.txt) supportant plusieurs courbes.
    entries: liste de dictionnaires {'id': GlobalID, 'variant': int}
    """
    filename = os.path.join(folder_output, filename)
    print(f" -> Création de {filename} avec {len(entries)} courbe(s)...")
    with open(filename, 'w') as f:
        # En-tête global
        f.write("$ Zpresc file generated automatically\n$\n$\n")
        f.write(f"{len(entries)}\n") # Nombre total de courbes zpresc
        
        # Boucle sur chaque entrée
        for i, entry in enumerate(entries):
            gid = entry['id']
            var = entry['variant']
            f.write(f"$ Zpresc {i+1} for Group ID {gid} (Variant {var})\n$\n$\n")
            f.write("2\n") # Nombre de points
            f.write(f"0.0     {Default_Z}\n")
            f.write(f"{Sim_duration} {Default_Z}\n")

def generate_ratcurve_file(filename, entries, folder_output="./"):
    filename = os.path.join(folder_output, filename)
    print(f" -> Création de {filename}...")
    with open(filename, 'w') as f:
        f.write("$ Rating Curve\n$\n$\n")
        f.write(f"{len(entries)}\n")
        for entry in entries:
            f.write("$\n$\n$\n2 0\n")
            f.write(f"1  1\n")
            f.write(f"1 1\n")

def generate_land_uses_file(filename, Land_uses_config = {1: {"friction": 0.033, "beta": 0.0} }):
    print(f" -> Création de {filename} (Format Fortran corrigé)...")
    with open(filename, 'w') as f:
        # Structure stricte 7 lignes d'en-tête
        f.write("$ Comment 1\n$ Comment 2\n$ Comment 3\n")
        f.write(f"{len(Land_uses_config)}\n")
        f.write("$ Comment 5\n$ Comment 6\n$ Comment 7\n")
        for code, props in Land_uses_config.items():
            f.write(f"{code} {props['friction']} {props['beta']}\n")

# =================================================


def gen_input(input_param= { "mesh_name":'channel.geo',
                                        "ts":100, 
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
                                        "w_tecplot":0,
                                        "w_vtk":0, 
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
                                        #"c_infil":0,
                                        "c_ic":0,
                                        "restart_min":0,
                                        "eps_min":0.0001}, folder_output="./" ):
    filename = os.path.join(folder_output, "input.txt")
    with open(filename, 'w') as f:
        f.write('!======================================================================================================================!')
        f.write('\n!Input File for Shallow-Water Model')
        f.write('\n!======================================================================================================================!')
        f.write('\n\n    &list_input')
        
        f.write('\n\n!======================================================================================================================!')
        f.write('\n!Mesh Type')
        f.write('\n!======================================================================================================================!')
        f.write("\n\n	mesh_type	=	'dassflow',")
        tmp_name = f"""\n    mesh_name    =   '{input_param["mesh_name"]}', """
        f.write(tmp_name)
        
        f.write('\n\n!======================================================================================================================!')
        f.write('\n!Simulation parameters')
        f.write('\n!======================================================================================================================!')
        tmp_name = f"""\n\n  ts                 =   {input_param["ts"]},      ! Simulation Time """
        f.write(tmp_name)
        tmp_name = f"""\n    dtw                =   {input_param["dtw"]},       ! data assimilation"""
        f.write(tmp_name)
        tmp_name = f"""\n    dtp                =   {input_param["dtp"]}, """
        f.write(tmp_name)    
        tmp_name = f"""\n    dta                =   {input_param["dta"]}, """
        f.write(tmp_name)    
        tmp_name = f"""\n    temp_scheme        =   '{input_param["temp_scheme"]}', """
        f.write(tmp_name)    
        tmp_name = f"""\n    spatial_scheme     =   '{input_param["spatial_scheme"]}', """
        f.write(tmp_name)    
        tmp_name = f"""\n    friction           =   {input_param["friction"]}, """
        f.write(tmp_name)    
        tmp_name = f"""\n    adapt_dt           =   {input_param["adapt_dt"]}, """
        f.write(tmp_name)    
        tmp_name = f"""\n    dt                 =   {input_param["dt"]}, """
        f.write(tmp_name)      
        tmp_name = f"""\n    cfl                =   {input_param["cfl"]}, """
        f.write(tmp_name)     
        tmp_name = f"""\n    feedback_inflow    =   {input_param["feedback_inflow"]}, """
        f.write(tmp_name)     
        tmp_name = f"""\n    coef_feedback      =   {input_param["coef_feedback"]}, """
        f.write(tmp_name)     
        f.write('\n\n!======================================================================================================================!')
        f.write('\n!PHYSICAL PARAMETERS')
        f.write('\n!======================================================================================================================!')
        tmp_name = f"""\n\n    g                  =   {input_param["g"]}, """
        f.write(tmp_name) 
        f.write('\n\n!======================================================================================================================!')
        f.write('\n!OUTPUT RESULTS')
        f.write('\n!======================================================================================================================!')     
        tmp_name = f"""\n\n    w_tecplot          =   {input_param["w_tecplot"]}, """
        f.write(tmp_name)   
        tmp_name = f"""\n      w_gnuplot          =   {input_param["w_gnuplot"]}, """
        f.write(tmp_name)   
        tmp_name = f"""\n      w_vtk              =   {input_param["w_vtk"]}, """
        f.write(tmp_name)   
        f.write('\n\n!======================================================================================================================!')
        f.write('\n!ASSIMILATION PARAMETER')
        f.write('\n!======================================================================================================================!') 
        tmp_name = f"""\n\n    w_obs                 =   {input_param["w_obs"]}, """
        f.write(tmp_name)   
        tmp_name = f"""\n    use_obs               =   {input_param["use_obs"]}, """
        f.write(tmp_name)   
        tmp_name = f"""\n    max_nt_for_adjoint    =   {input_param["max_nt_for_adjoint"]}, """
        f.write(tmp_name) 
       
        tmp_name = f"""\n\n    c_hydrograph          =   {input_param["c_hydrograph"]}, """
        f.write(tmp_name)    
        tmp_name = f"""\n    c_ratcurve            =   {input_param["c_ratcurve"]}, """
        f.write(tmp_name) 
        tmp_name = f"""\n    c_manning             =   {input_param["c_manning"]}, """
        f.write(tmp_name) 
        tmp_name = f"""\n    c_manning_beta        =   {input_param["c_manning_beta"]}, """
        f.write(tmp_name)
        tmp_name = f"""\n    c_bathy               =   {input_param["c_bathy"]}, """
        f.write(tmp_name)
        tmp_name = f"""\n    c_rain                =   {input_param["c_rain"]}, """
        f.write(tmp_name)
       # tmp_name = f"""\n    c_infil               =   {input_param["c_infil"]}, """
       # f.write(tmp_name)
        tmp_name = f"""\n    c_ic                  =   {input_param["c_ic"]}, """
        f.write(tmp_name)
        
        tmp_name = f"""\n\n  restart_min         =   {input_param["restart_min"]}, """
        f.write(tmp_name)                
        tmp_name = f"""\n    eps_min               =   {input_param["eps_min"]}, """
        f.write(tmp_name)    
        f.write("\n //  ")         


if __name__ == "__main__":
    print("=== CONVERSION Tools librairy ===")