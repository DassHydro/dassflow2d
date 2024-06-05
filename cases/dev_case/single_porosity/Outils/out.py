import numpy as numpy
import matplotlib.pyplot as plt
import re

def affichage(graphe,mesh,nc,ts,h,u,v):

    # Extraction des informations du maillage

    extract = re.findall(r'\d+', mesh)
    extract = [int(num) for num in extract]
    dx , L = extract

    qx = h[:] * u[:]
    qy = h[:] * v[:]

    # Trace de la hauteur d'eau

    comp = 0

    if (graphe[comp]) :

        plt.plot(range(0,L,dx),h,"b+",label='HLLC')
        plt.title("Hauteur d'eau à ts = {} s".format(ts))
        plt.xlabel("x (m)")
        plt.ylabel("h (m)")
        plt.legend()
    
        plt.show()
    
    comp = comp + 1

    # Trace de la vitesse u

    if (graphe[comp]) :

        plt.plot(range(0,L,dx),u,"r+",label='HLLC')
        plt.title("Vitesse selon x à ts = {} s".format(ts))
        plt.xlabel("x (m)")
        plt.ylabel("u (m/s)")
        plt.legend()
    
        plt.show()

    comp = comp + 1

    # Trace de la vitesse v

    if (graphe[comp]) :

        plt.plot(range(0,L,dx),v,"g+",label='HLLC')
        plt.title("Vitesse selon y à ts = {} s".format(ts))
        plt.xlabel("x (m)")
        plt.ylabel("v (m/s)")
        plt.legend()
    
        plt.show()

    comp = comp + 1

    # Trace du débit unitaire selon x

    if (graphe[comp]) :

        plt.plot(range(0,L,dx),qx,"k+",label='HLLC')
        plt.title("Débit unitaire selon x à ts = {} s".format(ts))
        plt.xlabel("x (m)")
        plt.ylabel("q (m²/s)")
        plt.legend()
    
        plt.show()

    comp = comp + 1

    # Trace du débit unitaire selon y

    if (graphe[comp]) :

        plt.plot(range(0,L,dx),qy,"c+",label='HLLC')
        plt.title("Débit unitaire selon x à ts = {} s".format(ts))
        plt.xlabel("x (m)")
        plt.ylabel("q (m²/s)")
        plt.legend()
    
        plt.show()

    comp = comp + 1
