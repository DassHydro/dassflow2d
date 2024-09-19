import numpy as np
import random
import os


###########################################################################

#      Génération du maillage

###########################################################################


def gen_mesh(type: str,L: float,dx: float) :

    # INPUTS :

    # type // type de maillage
    # L    // longueur du domaine
    # dx   // pas d'espace

    # OUTPUT :

    # Création du maillage souhaité sous la forme "type_dx=*_L=*.geo"


    #      ---------------------------------------
    #       Définition des noeuds et des éléments
    #      ---------------------------------------

    h = 10

    types = ['channel','slope','bump','flat&slope','flat_square','box']

    if type not in types :
        print('************************************\n')
        print(f"'{type}' is not a valid mesh type\n")
        print('************************************\n')
        print("Type mesh available : 'channel','slope','bump','flat&slope','flat_square' and 'box'\n")
        print('************************************\n')
        print("Stop of mesh generation \n")
        print('************************************')
        return
    
    # Maillage carré

    if ( type == 'box' or type == 'flat_square') :

        n =  int(np.floor(L/dx)) + 1    # Noeud par ligne
        nnode = n * n
        ncell = ( n - 1 ) * ( n - 1 )


        nodes = np.zeros((nnode,2))
        element =np.zeros((ncell,4),dtype=int)

        x_loc = np.zeros(n)

        for i in range(0,n) :

            x_loc[i] = i * L / ( n - 1 )

        for j in range(0,n) :
            for i in range(0,n) :

                # Définition des coordonnées des noeuds

                nodes[j * n + i, 0] = x_loc[i]
                nodes[j * n + i, 1] = x_loc[j]

        
        for j in range(1,n) :
            for i in range(1,n) :

                # Définition des éléments

                node1 = (j-1) * n + i
                node2 = node1 + 1
                node3 = node1 + n
                node4 = node3 + 1
                element[ (j-1) * (n-1) + i-1 , : ] = [node1,node2,node4,node3]


    # Maillage "rectangulaire"

    else : 
        
        ncell = int(np.floor(L/dx))
        nnode = 2 * ( ncell + 1 )

        # Définition des coordonnées des noeuds

        nodes = np.zeros((nnode,2))

        for i in range(0,nnode) :

            if ( i%2 == 0 ) :
                nodes[i,1] = 0.
                j = i

            else :
                nodes[i,1] = h
                j = i - 1

            nodes[i,0] = j * 0.5 * dx


        # Définition des éléments

        element = np.zeros((ncell,4),dtype=int)

        for i in range(0,ncell) :

            node1 = 1 + 2 * i
            node2 = node1 + 1
            node3 = node1 + 3
            node4 = node1 + 2
            
            element[i,:] = [node1,node2,node3,node4]



    #      ---------------------
    #       Conditions de bords
    #      ---------------------

    bc_file = 'bc.txt'

    bin_path = os.getcwd()

    bc_path = os.path.join(bin_path,bc_file)

    if not os.path.exists(bc_path) :
        print('************************************\n')
        print(" 'bc.txt' file is missing\n")
        print('************************************\n')
        print(" BC equals 0\n")
        print('************************************\n')

    else :

        n_bc_data = []

        with open(bc_path,'r') as bc_file:
            lines = bc_file.readlines()

            for line in lines :
                if line.startswith('!'):
                    continue

                parts = line.strip().split('\t')

                if (len(parts)==1) :
                    n_bc_data.append(parts)

        n_bc_data = np.array(n_bc_data)
        chaine = n_bc_data[0][0]
        n_bc_data = int(chaine)


    #      ---------------------------------
    #       Création et écriture du fichier
    #      ---------------------------------


    mesh_name = '{}_dx={}_L={}.geo'.format(type,dx,L)

    landtype = 1

    file = open(mesh_name,"w")

    file.write('# Generated mesh with my program gen_channel() ||| number of nodes | number of cells | mesh scale == 0 always\n')
    file.write(f'{nnode:d} ')
    file.write(f'{ncell:d} ')
    file.write(f'{0.:4e} \n')

    # Ecriture des noeuds

    file.write('# Nodes||| id node, x coord, y coord, bathy (x,y)\n')

    for i in range(0,nnode) :

        file.write(f'{i+1:6d} ')
        file.write(f'{nodes[i,0]:4e} ')
        file.write(f'{nodes[i,1]:4e} ')
        file.write(f'{0.:4e} \n')

    # Ecriture des cellules

    file.write('# Cells|||  id of cell, id node 1, id node 2, id  node 3,  id node 4 , land_type, bathymetry\n')

    for i in range(0,ncell) :

        file.write(f'{i+1:6d} ')
        file.write(f'{element[i,0]:6d} ')
        file.write(f'{element[i,1]:6d} ')
        file.write(f'{element[i,2]:6d} ')
        file.write(f'{element[i,3]:6d} ')
        file.write(f'{landtype:6d} ')
        file.write(f'{bathy(type,nodes[element[i,0],0],nodes[element[i,1],1],L):4e} \n')

    # Ecriture des conditions de bords

    file.write('# Boundaries\n')

    if ( n_bc_data == 0 ) :
        file.write('INLET ')
        file.write(f'{0:6d}')
        file.write(f'{0:6d}\n')

        file.write('OUTLET ')
        file.write(f'{0:6d}')
        file.write(f'{0:6d}\n')

    else :
        file.write('INLET ')
        file.write(f'{1:6d}') 
        file.write(f'{1:6d}\n')

        file.write(f'{1:6d} ')
        file.write(f'{1:6d} ')
        file.write(f'{1:6d} ')
        file.write(f'{0.:4e}')
        file.write(f'{1:6d}\n')

        file.write('OUTLET ')
        file.write(f'{1:6d}') 
        file.write(f'{1:6d}\n')

        file.write(f'{ncell:6d} ')
        file.write(f'{3:6d} ')
        file.write(f'{1:6d} ')
        file.write(f'{0.:4e} ')
        file.write(f'{1:6d}\n')

    file.close()

    return mesh_name



###########################################################################

#       Définition de la bathymétrie

###########################################################################


def bathy(type,x,y,L) :

    mil = L*0.5

    if ( type == 'channel' or type == 'flat_square') :

        res = 0.

    elif ( type == 'slope' ) :

        res = 0.01*x

    elif ( type == 'bump' ) :

        res = max(0.,-0.08*(x-(mil-10))*(x-(mil-10)) + 1.)

    elif ( type == 'flat&slope' ) :

        z_av = 0.5
        z_am = 1.
        ecart = 10

        if (x <= mil-ecart) :

            res = z_av

        elif ((x > mil-ecart) and (x < mil+ecart)) :

            a = (z_am - z_av) / ( 2 * ecart )

            b = z_av - a * ( mil - ecart )

            res = a * x + b

        else :

            res = z_am
    
    
    elif ( type == 'box' ) :

        res = random.random() 

    return res



###########################################################################

#       Destruction du maillage

###########################################################################


def delete_mesh(mesh_file: str) :

    bin_path = os.getcwd()

    file_path = os.path.join(bin_path,mesh_file)

    os.remove(file_path)

    return
