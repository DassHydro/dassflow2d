import numpy as np
import matplotlib.pyplot as plt
import re
import pyvista as pv
import os
import keyboard

def plot_dat(graphe,display_ref,mesh_name,ts,h,u,v,h0,u0,v0):

    # INPUTS :

    # graphe       : display decision vector  
    # => meaning   : 0 = to not display / 1 = to display
    # => order     : [ 'h' , 'u' , 'v' , 'qx' , 'qy' ]
    # display_ref  : display or not reference
    # mesh_name    : name of the mesh used
    # ts           : time simulation
    # h , u , v    : calculated solutions
    # h0 , u0 , v0 : initial data


    # IMPORTANT :
    
    # To ensure that the graphs are generated, the mesh name must be in the following format :
    #   'mesh_dx_L.geo' with dx = step space and L = length of the domain

    # Mesh information extraction

    extract = re.findall(r'\d+\.?\d*', mesh_name)
    extract = [float(num) if '.' in num else int(num) for num in extract]
    dx , L = extract
    L = int(L)
    x = np.arange(0,L,dx)

    save_step = round(1/dx)

    # Resizing to improve visibility
    
    if dx < 1. :
        x = x[::save_step]
        h0 = h0[::save_step]
        u0 = u0[::save_step]
        v0 = v0[::save_step]
        h = h[::save_step]
        u = u[::save_step]
        v = v[::save_step]

    q0x = h0[:] * u0[:]
    q0y = h0[:] * v0[:]
    qx = h[:] * u[:]
    qy = h[:] * v[:]

    len_x = len(x)

    # Reference data extraction

    if (display_ref) :

        file_name = 'reference.dat'
        path_bin = os.getcwd()
        path_file_ref = os.path.join(path_bin,file_name)

        data = np.loadtxt(path_file_ref, comments='#')

        h_ref = data[:,5]
        u_ref = data[:,7]
        v_ref = data[:,8]

        ratio = round(len(h_ref)/len_x)

        h_ref = h_ref[::ratio]
        u_ref = u_ref[::ratio]
        v_ref = v_ref[::ratio]
        qx_ref = h_ref[:] * u_ref[:]
        qy_ref = h_ref[:] * v_ref[:]


    # Water depth

    comp = 0

    if (graphe[comp]) :

        calc, = plt.plot(x,h,"b+",label='Solveur HLLC')
        init, = plt.plot(x,h0,"k.",label='Profil initial',markersize=2.5,alpha = 0.3)

        if (display_ref) :
            refe,  = plt.plot(x,h_ref,"k-",label='Profil de référence',alpha = 0.65)
            handles = [init,refe,calc]

        else : 
            handles = [init,calc]
        
        plt.title("Hauteur d'eau à ts = {} s".format(ts))
        plt.xlabel("x (m)")
        plt.ylabel("h (m)")
        plt.xlim(0,L)
        plt.legend(handles=handles,fontsize=9)
    
        plt.show()

    comp = comp + 1

    # x-velocity

    if (graphe[comp]) :

        calc, = plt.plot(x,u,"b+",label='Solveur HLLC')
        init, = plt.plot(x,u0,"k.",label='Profil initial',markersize=2.5,alpha = 0.3)

        if (display_ref) :
            refe,  = plt.plot(x,u_ref,"k-",label='Profil de référence',alpha = 0.65)
            handles = [init,refe,calc]

        else : 
            handles = [init,calc]
        
        plt.title("Vitesse selon x à ts = {} s".format(ts))
        plt.xlabel("x (m)")
        plt.ylabel("u (m/s)")
        plt.xlim(0,L)
        plt.legend(handles=handles,fontsize=9)
    
        plt.show()

    comp = comp + 1

    # y-velocity

    if (graphe[comp]) :

        calc, = plt.plot(x,v,"b+",label='Solveur HLLC')
        init, = plt.plot(x,v0,"k.",label='Profil initial',markersize=2.5,alpha = 0.3)

        if (display_ref) :
            refe,  = plt.plot(x,v_ref,"k-",label='Profil de référence',alpha = 0.65)
            handles = [init,refe,calc]

        else : 
            handles = [init,calc]
        
        plt.title("Vitesse selon y à ts = {} s".format(ts))
        plt.xlabel("x (m)")
        plt.ylabel("v (m/s)")
        plt.xlim(0,L)
        plt.legend(handles=handles,fontsize=9)
    
        plt.show()

    comp = comp + 1

    # Unit discharge along x

    if (graphe[comp]) :

        calc, = plt.plot(x,qx,"b+",label='Solveur HLLC')
        init, = plt.plot(x,q0x,"k.",label='Profil initial',markersize=2.5,alpha = 0.3)

        if (display_ref) :
            refe,  = plt.plot(x,qx_ref,"k-",label='Profil de référence',alpha = 0.65)
            handles = [init,refe,calc]

        else : 
            handles = [init,calc]
        
        plt.title("Débit unitaire selon x à ts = {} s".format(ts))
        plt.xlabel("x (m)")
        plt.ylabel("q (m²/s)")
        plt.xlim(0,L)
        plt.legend(handles=handles,fontsize=9)
    
        plt.show()

    comp = comp + 1

    # Unit discharge along y

    if (graphe[comp]) :

        calc, = plt.plot(x,qy,"b+",label='Solveur HLLC')
        init, = plt.plot(x,q0y,"k.",label='Profil initial',markersize=2.5,alpha = 0.3)

        if (display_ref) :
            refe,  = plt.plot(x,qy_ref,"k-",label='Profil de référence',alpha = 0.65)
            handles = [init,refe,calc]

        else : 
            handles = [init,calc]
        
        plt.title("Débit unitaire selon y à ts = {} s".format(ts))
        plt.xlabel("x (m)")
        plt.ylabel("q (m²/s)")
        plt.xlim(0,L)
        plt.legend(handles=handles,fontsize=9)
    
        plt.show()

    comp = comp + 1

    return

def plot_dat_visu(graphe,display_ref,mesh_name,ts,h,u,v,h0,u0,v0):

    # INPUTS :

    # graphe       : display decision vector  
    # => meaning   : 0 = to not display / 1 = to display
    # => order     : [ 'h' , 'u' , 'v' , 'qx' , 'qy' ]
    # display_ref  : display or not reference
    # mesh_name    : name of the mesh used
    # ts           : time simulation
    # h , u , v    : calculated solutions
    # h0 , u0 , v0 : initial data


    # IMPORTANT :
    
    # To ensure that the graphs are generated, the mesh name must be in the following format :
    #   'mesh_dx_L.geo' with dx = step space and L = length of the domain

    # Mesh information extraction

    extract = re.findall(r'\d+\.?\d*', mesh_name)
    extract = [float(num) if '.' in num else int(num) for num in extract]
    dx , L = extract
    L = int(L)
    x = np.arange(0,L,dx)

    save_step = round(1/dx)

    # Resizing to improve visibility
    
    if dx < 1. :
        x = x[::save_step]
        h0 = h0[::save_step]
        u0 = u0[::save_step]
        v0 = v0[::save_step]
        h = h[::save_step]
        u = u[::save_step]
        v = v[::save_step]

    q0x = h0[:] * u0[:]
    q0y = h0[:] * v0[:]
    qx = h[:] * u[:]
    qy = h[:] * v[:]

    len_x = len(x)

    # Reference data extraction

    if (display_ref) :

        file_name = 'reference.dat'
        path_bin = os.getcwd()
        path_file_ref = os.path.join(path_bin,file_name)

        data = np.loadtxt(path_file_ref, comments='#')

        h_ref = data[:,5]
        u_ref = data[:,7]
        v_ref = data[:,8]

        ratio = round(len(h_ref)/len_x)

        h_ref = h_ref[::ratio]
        u_ref = u_ref[::ratio]
        v_ref = v_ref[::ratio]
        qx_ref = h_ref[:] * u_ref[:]
        qy_ref = h_ref[:] * v_ref[:]


    # Water depth

    comp = 0

    if (graphe[comp]) :

        calc, = plt.plot(x,h,"b+",label='Profil calculé sans porosité')

        if (display_ref) :
            refe,  = plt.plot(x,h_ref,"r+",label='Profil calculé avec porosité',alpha = 0.65)
            handles = [calc,refe]

        else : 
            handles = [calc]
        
        plt.title("Hauteur d'eau à ts = {} s".format(ts))
        plt.xlabel("x (m)")
        plt.ylabel("h (m)")
        plt.xlim(0,L)
        plt.legend(handles=handles,fontsize=9)
    
        plt.show()

    comp = comp + 1

    # x-velocity

    if (graphe[comp]) :

        calc, = plt.plot(x,u,"b+",label='Profil calculé sans porosité')

        if (display_ref) :
            refe,  = plt.plot(x,u_ref,"r+",label='Profil calculé avec porosité',alpha = 0.65)
            handles = [calc,refe]

        else : 
            handles = [calc]
        
        plt.title("Vitesse selon x à ts = {} s".format(ts))
        plt.xlabel("x (m)")
        plt.ylabel("u (m/s)")
        plt.xlim(0,L)
        plt.legend(handles=handles,fontsize=9)
    
        plt.show()

    comp = comp + 1

    # y-velocity

    if (graphe[comp]) :

        calc, = plt.plot(x,v,"b+",label='Profil calculé sans porosité')

        if (display_ref) :
            refe,  = plt.plot(x,v_ref,"r+",label='Profil calculé avec porosité',alpha = 0.65)
            handles = [calc,refe]

        else : 
            handles = [calc]
        
        plt.title("Vitesse selon y à ts = {} s".format(ts))
        plt.xlabel("x (m)")
        plt.ylabel("v (m/s)")
        plt.xlim(0,L)
        plt.legend(handles=handles,fontsize=9)
    
        plt.show()

    comp = comp + 1

    # Unit discharge along x

    if (graphe[comp]) :

        calc, = plt.plot(x,qx,"b+",label='Profil calculé sans porosité')

        if (display_ref) :
            refe,  = plt.plot(x,qx_ref,"r+",label='Profil calculé avec porosité',alpha = 0.65)
            handles = [calc,refe]

        else : 
            handles = [calc]
        
        plt.title("Débit unitaire selon x à ts = {} s".format(ts))
        plt.xlabel("x (m)")
        plt.ylabel("q (m²/s)")
        plt.xlim(0,L)
        plt.legend(handles=handles,fontsize=9)
    
        plt.show()

    comp = comp + 1

    # Unit discharge along y

    if (graphe[comp]) :

        calc, = plt.plot(x,qy,"b+",label='Profil calculé sans porosité')

        if (display_ref) :
            refe,  = plt.plot(x,qy_ref,"r+",label='Profil calculé avec porosité',alpha = 0.65)
            handles = [calc,refe]

        else : 
            handles = [calc]
        
        plt.title("Débit unitaire selon y à ts = {} s".format(ts))
        plt.xlabel("x (m)")
        plt.ylabel("q (m²/s)")
        plt.xlim(0,L)
        plt.legend(handles=handles,fontsize=9)
    
        plt.show()

    comp = comp + 1

    return


def plot_vtk(test_case_path,scalar,time,ts,save):

    # INPUTS :

    # test_case_path : path to the test case
    # scalar         : charactere specifying the scalar to display 
    # time           : string specifying the date to display
    # ts             : time simulation
    # save           : boolean to save or not the graph


    # Path to 'result_*.vtk'

    file_vtk_name = 'result_{}.vtk'.format(time)

    file_vtk_path = os.path.join(test_case_path,'bin_A','res',file_vtk_name)

    # Title definition

    ts_sauv = ts

    if time == 'initial' : ts = 0

    if scalar == 'h' :
        title = "Water depth at ts = {} s".format(ts)
    elif scalar == 'u' :
        title = "X-velocity at ts = {} s".format(ts)
    elif scalar == 'v' :
        title = "Y-velocity at ts = {} s".format(ts)
    elif scalar == 'zs' :
        title = "Surface elevation at ts = {} s".format(ts)
    elif scalar == 'porosity' :
        title = "Porosity at ts = {} s".format(ts)
    else :
        return print("Error in the scalar to display : 'h','u','v','zs' or 'porosity' available. Scalar '{}' is not defined".format(scalar))
    
    # Data retrieval 

    data = pv.read(file_vtk_path)
    plotter = pv.Plotter()
    bounds = data.bounds


    plotter.add_mesh(data,scalars = scalar, cmap = 'viridis',
        scalar_bar_args={
        'title':scalar,
        'height' : 0.1,
        'position_x' : 0.5 ,
        'position_y' : 0.25 ,
        'width' : 0.4 },
        show_edges = False )
    plotter.add_text(title, position = ( 175 , 550 , 0 ), font_size=20)
    plotter.add_axes(interactive = True, xlabel = 'x (m)', ylabel='y (m)', zlabel='z')

    plotter.view_xy()
    plotter.show_bounds(font_size=13,
                        xtitle='X Axis ',
                        ytitle='Y Axis ',
                        ztitle='Z Axis ',
                        fmt = '%.0f')


    # Creation of the new image folder

    new_folder = 'Saved_pictures'

    new_folder_path = os.path.join(test_case_path,'bin_A',new_folder)

    if not os.path.exists(new_folder_path) :
        os.makedirs(new_folder_path)
    
    else :
        end_of_name = "{}.pdf".format(ts)
        files = os.listdir(new_folder_path)

        for file in files :
            file_path = os.path.join(new_folder_path,file)

            if (os.path.isfile(file_path) and (not file.endswith(end_of_name))) :
                os.remove(file_path)

    # Save the created file with a .pdf extension

    picture_name = "{}_ts={}.pdf".format(scalar,ts)

    picture_path = os.path.join(new_folder_path , picture_name)

    if save : 
        os.chdir(new_folder_path)
        plotter.save_graphic(picture_name)
        plotter.show_axes()
        plotter.show()

    else : 
        plotter.show_axes()
        plotter.show()

    ts = ts_sauv

    return 

def plot_vtk_box(test_case_path,scalar,time,ts,save):

    # INPUTS :

    # test_case_path : path to the test case
    # scalar         : charactere specifying the scalar to display 
    # time           : string specifying the date to display
    # ts             : time simulation
    # save           : boolean to save or not the graph


    # Path to 'result_*.vtk'

    file_vtk_name = 'result_{}.vtk'.format(time)

    file_vtk_path = os.path.join(test_case_path,'bin_A','res',file_vtk_name)

    # Title definition

    ts_sauv = ts

    if time == 'initial' : ts = 0

    if scalar == 'h' :
        title = "Water depth at ts = {} s".format(ts)
    elif scalar == 'u' :
        title = "X-velocity at ts = {} s".format(ts)
    elif scalar == 'v' :
        title = "Y-velocity at ts = {} s".format(ts)
    elif scalar == 'zs' :
        title = "Surface elevation at ts = {} s".format(ts)
    elif scalar == 'porosity' :
        title = "Porosity at ts = {} s".format(ts)
    else :
        return print("Error in the scalar to display : 'h','u','v','zs' or 'porosity' available. Scalar '{}' is not defined".format(scalar))
    
    # Data retrieval 

    data = pv.read(file_vtk_path)
    plotter = pv.Plotter()
    bounds = data.bounds


    plotter.add_mesh(data,scalars = scalar, cmap = 'viridis',
        scalar_bar_args={
        'title':scalar,
        'height' : 0.5,
        'position_x' : 0.85 ,
        'position_y' : 0.25 ,
        'width' : 0.1,
        'vertical':True },
        show_edges = False )
    plotter.add_text(title, position = ( 175 , 675 , 0 ), font_size=20)
    plotter.add_axes(interactive = True, xlabel = 'x (m)', ylabel='y (m)', zlabel='z')

    plotter.view_xy()
    plotter.show_bounds(font_size=13,
                        xtitle='X Axis ',
                        ytitle='Y Axis ',
                        ztitle='Z Axis ',
                        fmt = '%.0f')


    # Creation of the new image folder

    new_folder = 'Saved_pictures'

    new_folder_path = os.path.join(test_case_path,'bin_A',new_folder)

    if not os.path.exists(new_folder_path) :
        os.makedirs(new_folder_path)
    
    else :
        end_of_name = "{}.pdf".format(ts)
        files = os.listdir(new_folder_path)

        for file in files :
            file_path = os.path.join(new_folder_path,file)

            if (os.path.isfile(file_path) and (not file.endswith(end_of_name))) :
                os.remove(file_path)

    # Save the created file with a .pdf extension

    picture_name = "{}_ts={}.pdf".format(scalar,ts)

    picture_path = os.path.join(new_folder_path , picture_name)

    if save : 
        os.chdir(new_folder_path)
        plotter.save_graphic(picture_name)
        plotter.show_axes()
        plotter.show()

    else : 
        plotter.show_axes()
        plotter.show()

    ts = ts_sauv

    return 

