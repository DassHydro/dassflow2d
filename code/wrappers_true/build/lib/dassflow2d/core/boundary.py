import dassflow2d
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt 
import h5py

from matplotlib.colors import ListedColormap

class Boundary(object):
    """
    Boundary class, contains tables of values, mesh correspondance and inhering meshings:
        - self.bc, np.ndarray: table of values of shape (number of values, 2)
        - self.corresp, dict:  ids of correspondance between table and self.mesh_fortran object
        - self.mesh_fortran, dtype: inherit **mesh_fortran**
        - self.mesh_pyvista, pv.unstructured grid: inherit **mesh_pysta** at the same time
        
    """
    def __init__(self, mesh_pyvista, mesh_fortran): # mesh_pyvista, mesh_fortran
        
        """
        load  mesh_pyvista, mesh_fortran from input, source mesh correspondance  and table of values from fortran kernel
            - self.mesh_fortran, dtype:  inherit **mesh_fortran**
            - self.mesh_pyvista, pyvista.unstructuredgrid: inherit corresponding **mesh_pysta**
        """
        # NECESSARY INITIALISATION
        self.corresp = dict()
        # store inputs
        self.mesh_pyvista = mesh_pyvista
        self.mesh_fortran = mesh_fortran
        
        # load boundary data 
        self.get_table()                # from bc 
        self.get_mesh_corresp()         # from bc + mesh
        self.metadata = self.get_metadata()
        
    def get_table(self):
        """
        Get from fortran kernel tables of value of boundaries
        
        **return**
        
        self.bc, a derived type (bcs), which is a copy of model's bc values
        """
        my_bc = dassflow2d.wrapping.m_model.get_bc()
        self.bc = dassflow2d.wrapping.call_model.boundaries_copy(my_bc)
        
        
        
        
        # calculate correspondance between table of value and edges on mesh
    def get_mesh_corresp(self):
        """
        calculate correspondance between table of value and edges on mesh
        
        **return**
        
        self.corresp
        """
        def find_indices(lst, condition):
               return np.array([i for i, elem in enumerate(lst) if elem == condition], dtype = "int")


        def find_indices_double(lst, condition1, condition2):
               return np.array([i for i, elem in enumerate(lst) if elem == condition1 and elem == condition2], dtype = "int")

        typlim  = [] # typlim of the edge
        group = [] # gorup of the edge
        index = [] # global edge index
        for i in list( self.mesh_fortran.edgeb.indices):
            typlim.append( self.mesh_fortran.edgeb[i].typlim.decode("utf8"))
            group.append(  self.mesh_fortran.edgeb[i].group)
            index.append(  self.mesh_fortran.edgeb[i].ind)
            
        conditions = [ "wall" , "discharg1", "discharg2", "neumann", "transm", "ratcurve" ,"hpresc" ,"zpresc", "internal_1D", "internal_2D"]
        
        res = dict()

        for my_condition in conditions:
            
            nb = typlim.count(my_condition)
            if nb>0:
                    res[my_condition] = dict()
                    my_ind = find_indices(typlim, my_condition )
                    all_group =[]
                    for tmp_ind in my_ind:
                        all_group.append(group[tmp_ind])
                    if my_condition == "wall": #group set to 0
                        tofind = find_indices(typlim, my_condition )
                        ind_final = []
                        for i in tofind :
                            ind_final.append(index[i])
                        res[my_condition][0]= dict()
                        res[my_condition][0]["id_edge"]=np.array(ind_final)
                        res[my_condition][0]["id_edgeb"]=np.array(tofind)
                        res[my_condition][0]["id_grp"]=0
                        res[my_condition][0]["id_cell"]= np.array([self.mesh_fortran.edge[id_edge-1].cell[0] 
                            for k, id_edge in enumerate( res[my_condition][0]["id_edge"] ) ])
                           
                    else:    
                        unique_groups = np.unique(all_group)
                        nb_grp = len(unique_groups)   
                        
                        for i in range(nb_grp) :
                            pointeur = np.where(all_group ==unique_groups[i]) # size of all_group ) size of my_ind
                            ind_grp = my_ind[pointeur]
                            res[my_condition][unique_groups[i]] = dict()
                            res[my_condition][unique_groups[i]]["id_edgeb"]=np.array(ind_grp)
                            id_edge=np.take(index, ind_grp)
                            res[my_condition][unique_groups[i]]["id_edge"]= np.array(id_edge)
                            res[my_condition][unique_groups[i]]["id_grp"]=np.array(unique_groups[i])
                            res[my_condition][unique_groups[i]]["id_cell"]= np.array([self.mesh_fortran.edge[id_edge-1].cell[0] 
                            for k, id_edge in enumerate( res[my_condition][unique_groups[i]]["id_edge"] ) ])
    
        #print(res["discharg1"])
        #print(res["ratcurve"])
        self.corresp = res
        return(self)


        
        
        # get metadata for calculations,,loops etc..
    def get_metadata(self):
        """get metadata for calculations,,loops etc..
        
        **return**
        a dictionary containing metadata
        
        """
        res = dict()
        for i in self.corresp.keys():
            res[i] = len(self.corresp.get(i, 0))
        #print("number of each bc type, (wall always = 1) : ")
        return(res)
        


    def plot(self, what = "values",title_plot = "Mesh and boundaries", save_plot= False, filename="boundaries_location", notebook=False):
        """
        Plot boundary localisation on mesh

        Parameters
        ----------
        self: dassflowmodel  object
        
        what: str,  either
            - "values": plot hydrograms, rating curve, or any table of value imposed
            - "meshing": plot the localisation of boundary condition imposed

        title_plot: str, default = "title",
                    title of the plot

        save_plot: logical, default = False
                Do you wish to save the plot or not
                
        filename: str, default = "tmp"
                Path to fsave the graphic file to. (see Plotter.save_graphic)

        Examples
        --------
         no example

        See Also
        --------
        plot_grid: similar method but where you can define your own variable to plot

        """        
        if what == 'meshing' :
            bc_corresp = self.corresp
            
            def find_indices(lst, condition):
                       return np.array([i for i, elem in enumerate(lst) if elem == condition], dtype = "int")


            def find_indices_double(lst, condition1, condition2):
                       return np.array([i for i, elem in enumerate(lst) if elem == condition1 and elem == condition2], dtype = "int")

            conditions = [        "wall" , "discharg1", "discharg2", "neumann","transm",  "ratcurve"    ,"hpresc" ,"zpresc", "internal_1D", "internal_2D"]
            colors =     ["black","white", "red"      , "orange"   , "yellow" , "yellow", "yellow"      ,"#779231", "#329C46", "blue","green" ]

            newcol=["black"]
            new_cond=["none"]
            for k in bc_corresp.keys():
                id = find_indices(conditions, k)[0]+1
                newcol.append(colors[id])
                new_cond.append(k)


            legend_entries = []
            # not a bc condition given first
            legend_entries.append(['none', 'black'])

            edges = self.mesh_pyvista.extract_all_edges()
            edge_color= np.zeros(edges.n_faces)


            for i, k in enumerate(bc_corresp.keys()):
                color_index= i+1# color index correspond also to the color calue given to edge_color
                for g in bc_corresp[k].keys():
#                    if k == 'wall':                        
#                        edge_color[[x for x in bc_corresp[k][g]]  ] = color_index
#                    else:
                        edge_color[[x-1 for x in bc_corresp[k][g]["id_edge"]]] = color_index
                # add legend entry
                col = newcol[color_index]
                #print(i,k,  col)
                legend_entries.append([f'{k}', f'{col}'])

            # create scalar typlim (within pyvista polydata object )
            edges["typlim"] = edge_color

            if notebook==True:
                plotter = pv.Plotter(notebook=True)
            else:
                plotter = pv.Plotter(notebook=False)

            plotter.add_mesh(edges,scalars = "typlim", cmap=ListedColormap(newcol), show_scalar_bar=False)
            plotter.add_legend(legend_entries, face = "line")
    #        plotter.show_bounds(
    #            self.mesh_pyvista='back',
    #            location='outer',
    #            all_edges=False)
    #        plotter.show_bounds(padding=0.2, location = 'outer')
            plotter.add_title(title_plot, font_size=18, color=None, font=None, shadow=False)

            if(save_plot == True):
                plotter.save_graphic(filename=f"{filename}.svg", title = title_plot)
            
            if notebook==False:
                plotter.show(cpos = "xy")
            
            return(plotter)
            

        if what == 'values' :
            # only arguments are sefl and save_plot (title_plot and filename treated automaticaly)

            def plot_hyd(hydrograph, title_plot= "Inflow discharge", subtitle = "group 0", save_plot = False):
                plt.plot(hydrograph.t, hydrograph.q)
                plt.xlabel("time [s]")
                plt.ylabel("Discharge [m3/s]")
                plt.suptitle(title_plot)
                plt.title(f"group {hydrograph.group}")
                if(save_plot == True):
                     plt.savefig(f'hydrograph_{hydrograph.group}.png')
                plt.show()
                plt.close()

            def plot_ratcurve(rat, title_plot= "Rating curve", save_plot = False):
                plt.plot(rat.q, rat.h)
                plt.xlabel("Discharge [m3/s]")
                plt.ylabel("Water heigth [m]")
                plt.suptitle(title_plot)
                plt.title(f"group {rat.group}")
                if(save_plot == True):
                     plt.savefig(f'ratcurve_{rat.group}.png')
                plt.show()
                plt.close()

            def plot_zpresc(zpresc, title_plot= "prescribed water surface elevation", save_plot = False):
                plt.plot(zpresc.t, zpresc.z)
                plt.xlabel("time [s]")
                plt.ylabel("Imposed surface elevation [m]")
                plt.suptitle(title_plot)
                plt.title(f"group {zpresc.group}")
                if(save_plot == True):
                     plt.savefig(f'zspresc_{zpresc.group}.png')
                plt.show()
                plt.close()

            def plot_hpresc(hpresc, title_plot= "prescribed water heigth", save_plot = False):
                plt.plot(hpresc.t, hpresc.h)
                plt.xlabel("time [s]")
                plt.ylabel("Water heigth [m]")
                plt.suptitle(title_plot)
                plt.title(f"group {hpresc.group}")
                if(save_plot == True):
                     plt.savefig(f'hpresc_{hpresc.group}.png')
                plt.show()
                plt.close()

            metadata_boundary=self.get_metadata()
            
            for key in metadata_boundary.keys():
                nb = metadata_boundary[key]
                if key == 'ratcurve':
                    for i in range(nb):
                        plot_ratcurve(self.bc.rat[i], save_plot = save_plot)
                if key == 'hpresc':
                    for i in range(nb):
                        plot_hpresc(self.bc.hpresc[i], save_plot = save_plot)
                if key == 'zpresc':
                    for i in range(nb):
                        plot_zpresc(self.bc.zpresc[i], save_plot = save_plot)
                if key == 'discharg1' or  key == 'discharge2':
                    for i in range(nb):
                        plot_hyd(self.bc.hyd[i], save_plot = save_plot)



    def save(self, hdf5_path):
            """
            save boundary condition in hdf5 file specified
            
            The following groups are created:
            
                - f["boundary"]:
                    - f["values"] : contains table of values
                    - f["corresp"] : contains geometry correspondance, saved as an attributes
                    
            the following attributes are made available:
            
                -  f["boundary"]: print number of each existing bc
                    - "wall" , 
                    - "discharg1", 
                    - "discharg2",
                    -"neumann", 
                    - "transm",
                    - "ratcurve"    ,
                    - "hpresc" ,
                    - "zpresc", 
                    - "internal_1D", 
                    - "internal_2D"

            
            **parameters**
            
            
            hdf5_path: str, path to hdf5 file to saveS
            
            """
            def save_hyd(hydrograph,hdf5_path, group):
                table = np.ndarray(shape=(len(hydrograph.t),2))
                table[:,0]= hydrograph.t 
                table[:,1]= hydrograph.q
                
                # create dataset containing table in hdf5 file
                f = h5py.File(hdf5_path, "a")
                
                discharge=f["input"]["boundary"]["values"]["discharge"] 
                discharge.create_dataset(f"{group}", table.shape, 'f')
                
                discharge[f"{group}"][:,:]= table
                discharge[f"{group}"].attrs['groups'] = hydrograph.group                
                f.close()
                
            def save_ratcurve( rat, hdf5_path):
                table = np.ndarray(shape=(len( rat.h),2))
                table[:,0] =  rat.h 
                table[:,1] =  rat.q
                # create dataset containing table in hdf5 file
                f = h5py.File(hdf5_path, "a")   
                
                ratcurve=f["input"]["boundary"]["values"]["ratcurve"] 
                ratcurve.create_dataset(f"{rat.group}", table.shape, 'f')
                
                ratcurve[f"{rat.group}"][:,:]= table
                ratcurve[f"{rat.group}"].attrs['groups'] = rat.group                
                f.close()
                
                
            def save_zpresc(zpresc,  hdf5_path):
                table = np.ndarray(shape=(len(zpresc.t),2))
                table[:,0] = zpresc.t
                table[:,1] = zpresc.z
                # create dataset containing table in hdf5 file
                f = h5py.File(hdf5_path, "a")   
                
                zpresc=f["input"]["boundary"]["values"]["zpresc"] 
                zpresc.create_dataset(f"{zpresc.group}", table.shape, 'f')
                
                zpresc[f"{rat.group}"][:,:]= table
                zpresc[f"{rat.group}"].attrs['groups'] = zpresc.group    

            def save_hpresc(hpresc,  hdf5_path):
                table = np.ndarray(shape=(len(hpresc.t),2))
                table[:,0] = hpresc.t
                table[:,1] = hpresc.h
                # create dataset containing table in hdf5 file
                f = h5py.File(hdf5_path, "a")   
                
                zpresc=f["input"]["boundary"]["values"]["hpresc"] 
                zpresc.create_dataset(f"{hpresc.group}", table.shape, 'f')
                
                zpresc[f"{hpresc.group}"][:,:]= table
                zpresc[f"{hpresc.group}"]            
                f.close()

            def generate_values_groups(metadata, hdf5_path):
                f = h5py.File(hdf5_path, "a")   
                try:    # create input group where will be the meshing data
                        input_group = f.create_group("input")
                except: # if exist, cant be created, then it is loaded
                        input_group = f["input"]
                            
                boundary = input_group.create_group("boundary")
                values = boundary.create_group("values")
                # create empty groups
                for key in metadata.keys():                    
                    if key == 'ratcurve':
                        values.create_group("ratcurve")
                    if key == 'hpresc':
                        values.create_group("hpresc")
                    if key == 'zpresc':
                        values.create_group("zpresc")
                    if key == 'discharg1' or key == 'discharge2':
                        values.create_group("discharge")
                f.close()
                
             # Load data to save
            self.get_mesh_corresp()         # from bc + mesh
            self.bc = dassflow2d.wrapping.m_model.get_bc() # from bc 
            self.metadata = self.get_metadata()
             # group
            # 1) Generate empty groups for values
            generate_values_groups(metadata = self.metadata, hdf5_path = hdf5_path)
            
            # 2) fill them by creating datasets            
            for key in self.metadata.keys():
                nb = self.metadata[key]
                
                if key == 'ratcurve':
                    for i in range(nb):
                        save_ratcurve(self.bc.rat[i], hdf5_path = hdf5_path)
                if key == 'hpresc':
                    for i in range(nb):
                        save_hpresc(self.bc.hpresc[i], hdf5_path = hdf5_path)
                if key == 'zpresc':
                    for i in range(nb):
                        save_hpresc(self.bc.zpresc[i], hdf5_path = hdf5_path)
                if key == 'discharg1' or  key == 'discharge2':
                    for i in range(nb):
                        save_hyd(self.bc.hyd[i], hdf5_path = hdf5_path, group = i)
            
            # save geometry correspondance            
            f = h5py.File(hdf5_path, "a")   
            input_group= f["input"]
            boundary_group = input_group["boundary"]
            corresp_group = boundary_group.create_group("corresp")
            
            for key in self.corresp.keys():
                    key_group = corresp_group.create_group(key)
                    
                    for key_2 in self.corresp[key].keys():                    
                        key2_str = str(key_2)
                        data_group = key_group.create_group(key2_str)
                        
                        
                        for key_3 in self.corresp[key][key_2].keys():
                            key3_str = str(key_3)
                            #if str(self.corresp["internal_1D"][5]["id_cell"].dtype) == "int32" or str(self.corresp["internal_1D"][5]["id_cell"].dtype) == "int64":
                            try:   
                                dataset_length = len(self.corresp[key][key_2][key_3])
                            except:
                                 dataset_length=1
                            
                            my_dataset = data_group.create_dataset(key3_str, 
                                                        (dataset_length), 'f')
 
                               
                               
                    
                        my_dataset[:] = np.array(self.corresp[key][key_2][key_3])
            # save metatdata     
            for key in self.metadata.keys():   
                    boundary_group.attrs[str(key)] = self.metadata[key]
            f.close()
            
            
            
