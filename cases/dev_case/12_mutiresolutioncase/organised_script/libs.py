import numpy as np
import os

def read_mesh_from_textfile(path, read_boundary = False ):
# >>> Source File
    with open(path, 'r') as f:
        a = f.readlines()    
    
# >>> extract "metadata" + build list index
    nb_node, nb_cell, scale = np.asanyarray(a[1].split(), dtype = "float")
    read_node_indexes = np.arange(nb_node, dtype = "int") +3
    read_cell_indexes = np.arange(nb_cell, dtype = "int") + read_node_indexes[-1] + 2

# same for boundary
    if read_boundary:
      # inlet
        inlet_header_index = read_cell_indexes[-1]+2
        nb_cells_inlet = int(a[inlet_header_index].split()[1]    )
        nb_inlet = int(a[inlet_header_index].split()[2])
        read_inlet_indexes = np.arange(nb_cells_inlet, dtype = "int") +inlet_header_index +1
      # outlet
        outlet_header_index = read_inlet_indexes[-1] +1
        nb_cells_outlet = int(a[outlet_header_index].split()[1])
        nb_outlet = int(a[outlet_header_index].split()[2])
        read_outlet_indexes = np.arange(nb_cells_outlet, dtype = "int") +outlet_header_index +1
    
# >>> Define variables to fill
    nodes_x = []
    nodes_y = []
    nodes_z = []
    nodes_index = []    
    
    all_cell_index = []
    node_1 = []
    node_2 = []
    node_3 = []
    node_4 = []
    all_land_use =[]
    all_bathy = []
    
    if read_boundary:     
        INLET_index_fortran_cell =  []
        INLET_index_relative_edge = []
        INLET_ghost_bathy = []
        INLET_group = []        
        # 
        OUTLET_index_fortran_cell =  []
        OUTLET_index_relative_edge = []
        OUTLET_ghost_bathy = []
        OUTLET_group = []
        
# Nodes section ---------------------------------------
    for node_index in read_node_indexes:
        index_node, x_coord_node, y_coord_node, z_coord_node =  np.asanyarray(a[node_index].split(), dtype = "float") 
        nodes_x.append(x_coord_node)
        nodes_y.append(y_coord_node)
        nodes_z.append(z_coord_node)        
        nodes_index.append(int(index_node))
        
# Cells section ---------------------------------------        
    for cell_index in read_cell_indexes: 
        index_cell, node1, node2, node3, node4, land_use, bathy =  np.asanyarray(a[cell_index].split(), dtype = "float")      
        all_cell_index.append(int(index_cell))
        node_1.append(int(node1))
        node_2.append(int(node2))
        node_3.append(int(node3))
        node_4.append(int(node4))
        all_land_use.append( int(land_use ))
        all_bathy.append(bathy)
        
        
# Boundary section ---------------------------------------
    if read_boundary:
        # read inlet
        for inlet_index in read_inlet_indexes :
            index_fortran_cell, index_relative_edge , useless, ghost_bathy, group =  np.asanyarray(a[inlet_index].split(), dtype = "float") 
            index_fortran_cell = np.array(index_fortran_cell, dtype = "int")
            index_relative_edge = np.array(index_relative_edge, dtype = "int")
            useless = np.array(useless, dtype = "int")
            group = np.array(group, dtype = "int")
               
            
            INLET_index_fortran_cell.append(index_fortran_cell)
            INLET_index_relative_edge.append(index_relative_edge)
            INLET_ghost_bathy.append(ghost_bathy)
            INLET_group.append(group)   
            
        for outlet_index in read_outlet_indexes :
            index_fortran_cell, index_relative_edge , useless, ghost_bathy, group =  np.asanyarray(a[outlet_index].split(), dtype = "float")  
            
            index_fortran_cell = np.array(index_fortran_cell, dtype = "int")
            index_relative_edge = np.array(index_relative_edge, dtype = "int")
            useless = np.array(useless, dtype = "int")
            group = np.array(group, dtype = "int")
             
            OUTLET_index_fortran_cell.append(index_fortran_cell)
            OUTLET_index_relative_edge.append(index_relative_edge)
            OUTLET_ghost_bathy.append(ghost_bathy)
            OUTLET_group.append(group)   
        
        
        
        
    res = dict()
    
    res["header"] = dict()
    res["header"]["nb_node"]=nb_node
    res["header"]["nb_cell"]=nb_cell
    
    res["node"] = dict()
    res["node"]["x"]=nodes_x
    res["node"]["y"]=nodes_y
    res["node"]["z"]=nodes_z    
    res["node"]["index"]=nodes_index
    
    res["cell"] = dict()
    res["cell"]["index"]    = all_cell_index
    res["cell"]["node_1"]   = node_1
    res["cell"]["node_2"]   = node_2
    res["cell"]["node_3"]   = node_3
    res["cell"]["node_4"]   = node_4
    res["cell"]["land_use"] = all_land_use
    res["cell"]["bathy"]    = all_bathy

    if read_boundary:
        
        res["boundaries"] = dict()
        res["boundaries"]["header"] = dict()
        
        res["boundaries"]["header"]["nb_inflow"] = nb_inlet
        res["boundaries"]["header"]["nb_cell_inflow"] = nb_cells_inlet
        res["boundaries"]["header"]["nb_outflow"] = nb_outlet
        res["boundaries"]["header"]["nb_cell_outflow"] = nb_cells_outlet
        
        res["boundaries"]["inlet"] = dict()
        for key in range(len(read_inlet_indexes)):
             res["boundaries"]["inlet"][key]  = dict()
             
             res["boundaries"]["inlet"][key]["index_fortran_cell"] = INLET_index_fortran_cell[key]
             res["boundaries"]["inlet"][key]["index_relative_edge"] = INLET_index_relative_edge[key]
             res["boundaries"]["inlet"][key]["ghost_bathy"] = INLET_ghost_bathy[key]
             res["boundaries"]["inlet"][key]["group"] = INLET_group[key]
        
        res["boundaries"]["outlet"] = dict()
        for key in range(len(read_outlet_indexes)):
            res["boundaries"]["outlet"][key]  = dict()
            
            res["boundaries"]["outlet"][key]["index_fortran_cell"]  = OUTLET_index_fortran_cell[key]
            res["boundaries"]["outlet"][key]["index_relative_edge"] = OUTLET_index_relative_edge[key]
            res["boundaries"]["outlet"][key]["ghost_bathy"]         = OUTLET_ghost_bathy[key]
            res["boundaries"]["outlet"][key]["group"]               = OUTLET_group[key]

        
    return(res)
    
    
# write mesh from  mesh_dict, wich has the same structure as the one gotten with read_mesh
def gen_mesh_dassflow(mesh_dict, path = 'dassflow_mesh.geo'):
       
    nb_node = int(mesh_dict["header"]["nb_node"])
    nb_cell = int(mesh_dict["header"]["nb_cell"])
    
    with open(path, 'w') as f:
# Write HEADER
        f.write('# mesh generated from gmsh using gen_empty_mesh_dassflow() \n')  # python will convert \n to os.linesep
        f.write(f'{nb_node} {nb_cell} 0. \n') 
# Write Nodes block
        f.write('#Nodes||| id node, x coord, y coord, bathymetry\n ') 
        for id_node in range(nb_node):
              towrite_id_node = mesh_dict["node"]["index"][id_node]
              x = mesh_dict["node"]["x"][id_node]
              y = mesh_dict["node"]["y"][id_node]
              z = mesh_dict["node"]["z"][id_node]
              f.write(f'{towrite_id_node} {x} {y} {z} \n') 
# Write cell block
        f.write('#cells||| id cell, id_node1, id_node2, id_node3, id_node4, patch_manning, bathymetry \n ') 
                
        for id_cell in range(nb_cell):              
              cell_index = mesh_dict["cell"]["index"][id_cell]
              id_node1 =  mesh_dict["cell"]["node_1"][id_cell]
              id_node2 =  mesh_dict["cell"]["node_2"][id_cell]
              id_node3 = mesh_dict["cell"]["node_3"][id_cell]
              id_node4 =  mesh_dict["cell"]["node_4"][id_cell]
              land_use =   mesh_dict["cell"]["land_use"][id_cell]
              bathy =   mesh_dict["cell"]["bathy"][id_cell]                  
              f.write(f'{cell_index} {id_node1} {id_node2} {id_node3} {id_node4} {land_use} {bathy}\n') 

        if "boundaries" in mesh_dict.keys():

              f.write("# boundaries \n")                      

              f.write(f'INLET  {mesh_dict["boundaries"]["header"]["nb_cell_inflow"]} {mesh_dict["boundaries"]["header"]["nb_inflow"]}\n')
              for key, values in mesh_dict["boundaries"]["inlet"].items():
                  f.write(f'{values["index_fortran_cell"]} {values["index_relative_edge"]} {values["index_relative_edge"]} {values["ghost_bathy"]}  {values["group"]}\n')                  

              f.write(f'OUTLET  {mesh_dict["boundaries"]["header"]["nb_cell_outflow"]} {mesh_dict["boundaries"]["header"]["nb_outflow"]} \n')
              for key, values in mesh_dict["boundaries"]["outlet"].items():
                  f.write(f'{values["index_fortran_cell"]} {values["index_relative_edge"]} {values["index_relative_edge"]} {values["ghost_bathy"]}  {values["group"]}\n')
        f.close()
    return()
    
   

def gen_land_use(path, patch_values):
    with open(path, 'w') as f:
        f.write('######################## \n')
        f.write('#  Number of Land Uses \n')
        f.write('######################## \n')                
        f.write(str(len(patch_values)) + "\n")
        
        f.write('########################\n')
        f.write('# List of land_use \n')
        f.write('########################\n')
        for i in range(len(patch_values)):
            f.write(f"{i+1} {patch_values[i]} 0 \n")
            
    return()
        
        


def rewrite_friction_dassflow(bin_dir, mesh_name, land_use_name, mesh_dict, mesh_corresp, land_use_value):        
    mesh_dict["cell"]["land_use"] = mesh_corresp    
    gen_mesh_dassflow(mesh_dict = mesh_dict, path = bin_dir +"/" + mesh_name)
    gen_land_use(path = bin_dir +"/" + land_use_name,patch_values =  land_use_value)
    
    print("new files created :")
    
    print(f">>> {bin_dir}/{mesh_name}")
    print(f">>> {bin_dir}/{land_use_name}")
    return()
     

def source_min_res(min_dir):
    min_files = os.listdir(min_dir)
    tmp = [x[:-4] for x in min_files]
    def find_indices(lst, condition):
                           return [int(i) for i, elem in enumerate(lst) if elem == condition]
    
    min_files = [x for i, x in enumerate(min_files) if i in find_indices(tmp, "manning")]
    indices = np.array([int(x[-3:]) for x in min_files])
    
    min_file = [x for i, x in enumerate(min_files) if i ==np.array(np.where(indices == max(indices)))[0]][0]
    min_path = min_dir +"/"+ min_file
    
    f = open(min_path, "r")
    lines = f.readlines()
    inf_patch_corresp = []
    inf_patch_value = []
    for i in lines:
        res = i.split()
        inf_patch_corresp.append(res[0])
        inf_patch_value.append(res[1])
    optimized_patch_value = np.array(inf_patch_value, dtype = "float64")
    optimized_patch_corresp =  np.array(inf_patch_corresp, dtype = "int")
    return(optimized_patch_value, optimized_patch_corresp)
    
    

def build_p_matrix(tmp, previous_patch, new_patch):
    # tmp = all_patches as ndarray
    # previous_patch, id of column of numpy ndarray corresponding to previous patch
    # new_patch, id of column of numpy ndarray corresponding to previous patch
    X_id = np.unique(tmp[:, previous_patch])
    Y_id = np.unique(tmp[:, new_patch])    
    MATRIX_P = np.zeros(shape = (len(Y_id), len(X_id)), dtype = 'int' )
    
    for i in X_id:
        tmp_group = np.where(tmp[:, previous_patch] == i)  
        new_group = tmp[tmp_group, new_patch]
        id_to_fill_new = np.unique(new_group)        
        
        MATRIX_P[id_to_fill_new-1,i-1] = 1
    return(MATRIX_P)







