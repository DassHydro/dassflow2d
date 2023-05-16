import dassflow2d
import numpy as np
import pyvista as pv
import h5py


class Meshing(object):
    """
    Meshing is composed of:
               - self.mesh_fortran, from msh fortran derived type, fortran like manipulations.
               - self.mesh_pyvista, pyvista.unstructuredgrid object
    
    various plot method are available
    """
    def __init__(self, mesh_fortran) :
		
               self.mesh_fortran = mesh_fortran
               self.mesh_pyvista= self.build_grid()

	#=================================================================#
	# BUILD PYVISTA GRID                                              #
	#=================================================================#
	# my_mesh = DassFlowModel.kernel.mesh
	# self : dassflowmodel
	#    what = "grid" : build python_meshing["grid"] object
	#         = "python" : build python_meshing["fortran_like"] object
    def build_grid(self):
        """
        Build pyvista.unstructuredgrid object


        @param: kernel = model.kernel object
        @param what = what parameter to return
                - "pyvista" : pyvista grid
                - "arrays" : arrays needed by pyvista.grid()
                 - "None"    : just store in self.grid
        generate pyvista.Unstructuredgrid from self.model.mesh

        return:
        """

        my_mesh = self.mesh_fortran
        # The number of cells
        nc = my_mesh.nc
        # --------------------------------------------- #
        # define cell_type of each cell(triangle or quadrangle) --> Only quandangle considered
        # --------------------------------------------- #
        cell_type = np.empty(shape= (my_mesh.nc) )
        cell_type[:] = 9
        # --------------------------------------------- #
        # define cells and points
        # --------------------------------------------- #
        cells =  np.empty(shape= (0) )
        points = np.empty(shape= (0, 3) )

        existing_id_node = []
        for id_cell in range(nc ):
            id_nodes = my_mesh.cell[id_cell].node -1

            if id_nodes[3] == id_nodes[0]: # if this is a triangular mesh, 4th node = 1st one
                    id_nodes = id_nodes[0:3]
                    cell_type[id_cell] = 5

            all_nodes =np.empty(shape=(len(id_nodes),3) )
            i=-1
            # update necessary cell infomation
            cells = np.append(cells, len(id_nodes))
            cells = np.append(cells, id_nodes)
            cells = cells.astype(int)

            # update necessary node information
        for my_node in range(my_mesh.nn) :
                xyz= [0,0,0]
                xyz[0] = my_mesh.node[my_node].coord.x
                xyz[1] = my_mesh.node[my_node].coord.y
                xyz[2] = 0 # usefull if 3D grid
                points = np.append(points, [xyz], axis = 0)
        # --------------------------------------------- #
        # build mesh
        # --------------------------------------------- #
        grid = pv.UnstructuredGrid(cells, cell_type, points)
        return(grid)




#=================================================================#
# PLOT PYVISTA GRID                                               #
#=================================================================#

    def plot(self, my_scalar = None, title_scale_bar =" ", title_plot = " ", xlabel = " ", ylabel = " ", zlabel =" ", notebook=False):
        """
        Plot the mesh and define values on all cells for plot
        
        Parameters
        ----------
        
        my_scalar: numpy.array, default None
                  Numpy array corresponding to the values on each cell.
                  (the array must be correctly ordoned and of correct size : equal to the number of cells)
                  
        title_scale_bar: str, default = " "
                        Title of scale bar
                        
        title_plot: default = " "
                         title of the plot
                         
        Examples
        --------
         no example
         
        See Also
        --------
        plot_var: similar method but taking advantage of hdf5 file to get wished scalar to plot
        """
        my_mesh = self.mesh_pyvista
        
        # --------------------------------------------- #
        # plot only mesh
        # --------------------------------------------- #
        if my_scalar is None :
            my_scalar = np.empty(shape= (self.mesh_pyvista.number_of_cells) )
            my_scalar[:] = 0
            sargs = dict(vertical = True, title = title_scale_bar, position_y=0.3, position_x = 0.8  )
            # Plotter instance
            if notebook==True:
                plotter = pv.Plotter(notebook=True)
            else:
                plotter = pv.Plotter(notebook=False)

            # mesh plot
            plotter.add_mesh(my_mesh, scalars = my_scalar, scalar_bar_args=sargs,
                     show_edges = True)
            plotter.view_xy()
        #   plotter.show_axes()
            plotter.show_bounds()
            plotter.add_title(title_plot, font_size=18, color=None, font=None, shadow=False)
            if notebook==False:
                plotter.show()
            
        # --------------------------------------------- #
        # plot mesh with variable to show
        # --------------------------------------------- #
        else:
            sargs = dict(vertical = True, title = title_scale_bar, position_y=0.3, position_x = 0.8  )
            # Plotter instance
            if notebook==True:
                plotter = pv.Plotter(notebook=True)
            else:
                plotter = pv.Plotter(notebook=False)
            # mesh plot
            plotter.add_mesh(my_mesh, scalars = my_scalar, scalar_bar_args=sargs,
                     show_edges = True)
            plotter.add_axes(xlabel = xlabel, ylabel = ylabel)
            plotter.view_xy()
        #    plotter.show_axes()
            plotter.show_bounds()
            plotter.add_title(title_plot, font_size=18, color=None, font=None, shadow=False)
            if notebook==False:
                plotter.show()

        return(plotter)


#================================================================================================================================#
#  PLOT MESH
#================================================================================================================================#

    def plot_dev(self, what = "node"):
        """
        Plot the mesh and define values on all cells for plot

        Parameters
        ----------

        what: str, defknes fortran kernel values you want to access, either:           
            - "node"
            - "point"
            - "edge"
        """


        if what == "cell" :

            plotter = pv.Plotter(off_screen=False,  notebook=False)
            points = self.mesh_pyvista.cell_centers().points
            lab = points[:,(0,1)].tolist()

            all_cell_data = []
            scalar = np.zeros(shape=len(lab))
            for i in range(len(lab)) :
                j=i+1 # j is fortran index
                id_node = self.mesh_fortran.cell[i].node

                id_edge = self.mesh_fortran.cell[i].edge
                id_ngb_cell = self.mesh_fortran.cell[i].cell
                nb_ngb = self.mesh_fortran.cell[i].nbed
                is_boundary = self.mesh_fortran.cell[i].boundary ==1
                surf = self.mesh_fortran.cell[i].surf
                invsurf = self.mesh_fortran.cell[i].invsurf
                peri = self.mesh_fortran.cell[i].peri
                grav_x=self.mesh_fortran.cell[i].grav.x
                grav_y=self.mesh_fortran.cell[i].grav.y
                rain=self.mesh_fortran.cell[i].rain
                # add manualy (this is the rule), ghost cell's associated number is mesh.nc +

                lab[i] = f"""ID={j} \n COORD={lab[i]} \n NODES={id_node} \n EDGE={id_edge} \n NGB_CELL={id_ngb_cell}
                \n ID={j} neigbours={nb_ngb}
                \n boundary = {is_boundary}
                \n SURF={surf} invsurf={invsurf} peri={peri}
                \n grav_x={grav_x} grav_y={grav_y}     \n rain={rain}"""

                cell_data = {'id': j,
                         'id_node':id_node,
                         'id_edge':id_edge,
                         'id_ngb_cell':id_ngb_cell,
                         'nb_ngb':nb_ngb,
                         'is_boundary':is_boundary,
                         'surf':surf,
                         'invsurf':invsurf,
                         'peri':peri,
                         'grav_x':grav_x,
                         'grav_y':grav_y,
                         'rain':rain
                        }

                all_cell_data.append(cell_data)

                if is_boundary :
                    scalar[i] = 1
                else:
                    scalar[i] = 0

            plotter.add_mesh(self.mesh_pyvista, scalars = scalar , show_edges=True, color="tan")
            plotter.add_point_labels(points, lab, point_size=10, point_color ="red", font_size=15)
            plotter.show(cpos="xy")

            return(all_cell_data)



        if what == "edge":


            plotter = pv.Plotter(off_screen=False,  notebook=False)
            plotter.add_mesh(self.mesh_pyvista, show_edges=True, color="tan")
            # Add labels to points on the yz plane (where x == 0)
            #points = self.mesh_pyvista.cell_centers().points
            points = self.mesh_pyvista.extract_all_edges()
            #lab = points.tolist()
            edges = self.mesh_pyvista.extract_all_edges()
            edges.lines  # line connectivity stored here
            edges.plot(scalars=np.random.random(edges.n_faces),
                   line_width=10, cmap='jet', cpos="xy",
                   off_screen=False, notebook=False)


            all_edge_data = []
            lab = []
            points = np.ndarray(shape=(self.mesh_fortran.ne,3))

            for i in range(self.mesh_fortran.ne) :
                j=i+1 # j is fortran index
                id_node = self.mesh_fortran.edge[i].node

                id_cell = self.mesh_fortran.edge[i].cell
                #is_1D2D = self.mesh_fortran.edge[i].cell1D2D
                is_boundary = self.mesh_fortran.edge[i].boundary ==1
                subdomain   = self.mesh_fortran.edge[i].subdomain
                id_lim      = self.mesh_fortran.edge[i].lim
                length      = self.mesh_fortran.edge[i].length
                center_x    = self.mesh_fortran.edge[i].center.x
                center_y    = self.mesh_fortran.edge[i].center.y
                normal_x    = self.mesh_fortran.edge[i].normal.x
                normal_y    = self.mesh_fortran.edge[i].normal.y
                tangent_x     = self.mesh_fortran.edge[i].tangent.x
                tangent_y     = self.mesh_fortran.edge[i].tangent.y
                vcell_x      = self.mesh_fortran.edge[i].vcell.x
                vcell_y       = self.mesh_fortran.edge[i].vcell.y
                v_edge_cell_1x = self.mesh_fortran.edge[i].v_edge_cell[1].x
                v_edge_cell_1y = self.mesh_fortran.edge[i].v_edge_cell[1].y
                v_edge_cell_2x = self.mesh_fortran.edge[i].v_edge_cell[0].x
                v_edge_cell_2y = self.mesh_fortran.edge[i].v_edge_cell[0].y
                #1D2D={is_1D2D}
                lab.append( f"""ID={j}
                is_boundary={is_boundary}
                subdomain={subdomain}
                id_lim={id_lim}
                length={length}
                center_x={center_x}
                center_y={center_y}
                normal_x={normal_x}
                normal_y={normal_y}
                tangent_x={tangent_x}
                tangent_y={tangent_y}
                vcell_x={vcell_x}
                vcell_y={vcell_y}
                v_edge_cell_1_x={v_edge_cell_1x}
                v_edge_cell1__y={v_edge_cell_1y}
                v_edge_cell_2_x={v_edge_cell_2x}
                v_edge_cell_1_x={v_edge_cell_2y}   """ )

                points[i,0] = center_x
                points[i,1] = center_y
                points[i,2] = 0

                edge_data = {'id': j,
                         'id_node':id_node,
                         'id_cell':id_cell,
                         'is_boundary':is_boundary,
                         'subdomain':subdomain,
                         'id_lim':length,
                         'center_x':center_x,
                         'center_y':center_y,
                         'normal_x':normal_x,
                         'normal_y':normal_y,
                         'tangent_x':tangent_x,
                         'tangent_y':tangent_y,
                         'vcell_x':vcell_x,
                         'vcell_y':vcell_y,
                         'v_edge_cell_1x':v_edge_cell_1x,
                         'v_edge_cell_1y':v_edge_cell_1y,
                         'v_edge_cell_1x':v_edge_cell_1x,
                         'v_edge_cell_1y':v_edge_cell_1y
                         #'points':points
                        }
                all_edge_data.append(edge_data)

            plotter.add_point_labels(points, lab, point_size=10, point_color ="red", font_size=15)
            plotter.show(cpos="xy")

            return(all_edge_data)





        if what == "node" :


            plotter = pv.Plotter(off_screen=False,  notebook=False)
            plotter.add_mesh(self.mesh_pyvista, show_edges=True, color="tan")
            # Add labels to points on the yz plane (where x == 0)
            points = self.mesh_pyvista.points
            lab = points.tolist()

            all_node_data = [] # to be done


            for i in range(len(lab)) :
                j=i+1 # j is fortran index
                id_node = j
                id_edge = self.mesh_fortran.node[i].edge
                grav_x=self.mesh_fortran.node[i].coord.x
                grav_y=self.mesh_fortran.node[i].coord.y
                is_boundary=self.mesh_fortran.node[i].boundary
                id_lim=self.mesh_fortran.node[i].lim

                lab[i] = f"""ID_NODE={id_node} \n COORD={lab[i]} \n NODES={id_node} \n EDGE={id_edge}
                \n grav_x={grav_x} grav_y={grav_y}
                \n boundary = {is_boundary}     \n id_lim={id_lim}"""

                node_data = {'id': j,
                         'id_node':id_node,
                         'id_edge':id_edge,
                         'grav_x':grav_x,
                         'grav_y':grav_y,
                         'id_lim':id_lim,
                         'is_boundary':is_boundary
                        }
                all_node_data.append(node_data)


            plotter.add_point_labels(points, lab, point_size=10, point_color ="red", font_size=15)
            plotter.show(cpos="xy")

            return(all_node_data)
            
            # fmt : format hdf5/vtk

    def save(self, hdf5_path, fmt="hdf5" ):
        """
        save meshing
        
        **parameters**
        
            hdf5_path, str: path to the hdf5 file to save
            
            fmt, str: the format to save mesh
            - vtk, from msh fortran derived type, fortran like manipulations.
            - hdf5, from msh fortran dtype, save necessary information to write mesh.geo file
            (to write again mesh.geo file, boundary conditions are also necessary)
            
        """
        if fmt == "vtk":
            self.mesh_pyvista.save(f"{hdf5_path}_mesh.vtk", binary = False)
        if fmt == "hdf5":
#-------
# Save cells data
#-------
            
    # get some stats
            nb_cell = self.mesh_fortran.nc #number of cells  
                      
    # >>> create group and dataset in hdf5 file            
            f = h5py.File(hdf5_path, "a") 
            
            try:    # create input group where will be the meshing data
                        input = f.create_group("input")
            except: # if exist, cant be created, then it is loaded
                        input = f["input"]
                            
            meshing = input.create_group("meshing")
            # ---> create dataset cell
            # number of columns =  5:
            # ---->  id_fortran_cell, id_fortran_node1, id_fortran_node2,id_fortran_node3,id_fortran_node4, 
            # number of row:
            # variable nc
            cell = meshing.create_dataset("cell", (nb_cell,5), 'f')
            
    # >>> calculate + save group for cell data
    
             # these lines should be the same as in initialisation
            all_cell_data = []
            
            
            for i in range(nb_cell) :
				
                j=i+1 # j is fortran index
                id_node = self.mesh_fortran.cell[i].node    
                #id_edge = self.mesh_fortran.cell[i].edge
                #id_ngb_cell = self.mesh_fortran.cell[i].cell
                #nb_ngb = self.mesh_fortran.cell[i].nbed
                #is_boundary = self.mesh_fortran.cell[i].boundary ==1
                #surf = self.mesh_fortran.cell[i].surf
                #invsurf = self.mesh_fortran.cell[i].invsurf
                #peri = self.mesh_fortran.cell[i].peri
                #grav_x=self.mesh_fortran.cell[i].grav.x
                #grav_y=self.mesh_fortran.cell[i].grav.y
                #rain=self.mesh_fortran.cell[i].rain
                

                cell[i,:] = np.array([j, id_node[0],id_node[1],id_node[2],id_node[3]])
                
    
#-------
# Save node data
#-------
    
    # get some stats
            nb_node = self.mesh_fortran.nn #number of cells  
            
            
            # ---> create dataset cell
            # number of columns =  4:
            # ---->  id_fortran_node, x_coord, y_coord, bathymetry =0 
            # number of row:
            # variable nb_node
            node = meshing.create_dataset("node", (nb_node,4), 'f')
           
    
    
            for i in range(nb_node) :
                j=i+1 # j is fortran index
                
                grav_x=self.mesh_fortran.node[i].coord.x
                grav_y=self.mesh_fortran.node[i].coord.y
                #id_node = j
                #id_edge = self.mesh_fortran.node[i].edge
                #is_boundary=self.mesh_fortran.node[i].boundary
                #id_lim=self.mesh_fortran.node[i].lim

                node[i,:] = np.array([j, grav_x, grav_y, 0])
            
            f.close()






#######################
# usefull tools
# get_id_boundary_edge_from_cell :
#  return the id (from 1 to 4) of the edge on which the boundary is positionned  ( for mesh.geo file in inlet and outlet lines )

# STEP 1) extract edge data
# - on all edges, identify wich are boundary (and not wall),  extract all the id
            
# STEP 2) treatment for all cell data  + return

# loop on cells: for each cells ,
#   - get the id of cells which are boundary
#   - if the cell is boundary:
#   	-  extract the ordered list of of edges
#   	- find eventual match wich ids obteined in step 1 (if not this is a boundary wall)
#       	- if match -->  return id of match (from 1 to 4) of connection
#       	- if not match -->  print no inflow or outflow bondary connected

 
    def get_id_edgeb(self, id_cell, id_edge):  
        """
        # get id of boundary edge from onecell :
        #        
        #  return the id (from 1 to 4) of the edge on which the boundary is positionned  ( for mesh.geo file in inlet and outlet lines )
        #
        # STEP 1) extract edge data
        # - on all edges, identify wich are boundary (and not wall),  extract all the id
        #
        # STEP 2) treatment for all cell data  + return
        #
        # loop on cells: for each cells ,
        #   - get the id of cells which are boundary
        #   - if the cell is boundary:
        #       -  extract the ordered list of of edges
        #       - find eventual match wich ids obteined in step 1 (if not this is a boundary wall)
        #           - if match -->  return id of match (from 1 to 4) of connection
        #           - if not match -->  print no inflow or outflow bondary connected

        """

        cell = self.mesh_fortran.cell[id_cell].edge
        
        res = np.where(cell == id_edge)
        #print(res)
        if len(res) >0:
            return(res)
        else:
            print(f"cell {id_cell} and edge {id_edge} are not connected",)
        
        
