#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 18:59:18 2022

@author: livillenave
"""
import numpy as np
import dassflow2d
import os 
import h5py
import matplotlib.pyplot as plt

class Min(object):
        """
        Min Class
        
        contains:
            - self.nb_ite, integer: the number of iteration performed
            - self.all_ite, np.array, size(nb_ite):  
                self.all_ite[i]=i+1, the id of each iteration
            - self.j, np.array, size(nb_ite):  
                contains for each iteration the cost function estimation
            - self.gradj, np.ndarray of size(nb_ite,nb_param_infered):  
                contains for each iteration the gradient of cost function  for each parameter infered
            - self.gradj_name, list of strings size (nb_param_infered): 
                gives the name coresponding to the column of min.gradj 
            - self.param,  dictionary infered parameters at each iteration. 
                the dictionary is accesible as: min.param[name_param][id_param] 
                id param is used in the case of multiple hydrograph for exemple, else by default id_parm=1 and is unique
        """
        
        def __init__(self, bin_dir):
            self.bin_dir = bin_dir
            self.min_dir = f"{bin_dir}/min/"
            self.gradj_name=[]
            self.param=dict()
            print("")
            
        def source_all(self):        
            ### load min cost.txt
            min_cost = np.loadtxt(f"{self.min_dir}/min_cost.txt")    
            # list all files in /min_dir/ repertory
            files = os.listdir(self.min_dir)
            ### set all values that can be set with the file
            self.nb_ite =  int(min_cost[:,0][-1] )
            self.all_ite =  min_cost[:,0] 
            self.nb_param = min_cost.shape[1]-2
            self.j = min_cost[:,1] 
            self.gradj = min_cost[:, -self.nb_param:]        
            
            ### define to inferable  params
            # warning (TODO): the ordering must be the same as in fortran var_2_control (and control_2_var because its the same order)
            all_param_inferable = ['manning', 'manning_beta', 'bathy', 'hydrograph', 'ratcurve', 'rain', 'ic']        
                    
            
            self.param=dict()
            
            # loop on all inferable param
            for param_name in all_param_inferable: 
                ite_file = [x for x in files if param_name in x]
                
                if param_name == "hydrograph":
                    id_start_nb_file = -7
                    id_end_nb_file = -4
                elif param_name == "manning":                    
                    id_start_nb_file = -4
                    id_end_nb_file = -1
                else:                    
                    id_start_nb_file = -4
                    id_end_nb_file = -1
                    
                if len(ite_file)>0: 
                    if param_name == 'hydrograph':
                    # get nb unique files
                        nb_file=len(np.unique(np.array([x[id_start_nb_file:id_end_nb_file] for x in ite_file])))
                        fast_param_name='Q' 
                        for i in range(nb_file):
                            self.gradj_name.append(fast_param_name+f"_{i+1}"      )                  
                    else:
                            nb_file =1
                            fast_param_name = param_name                            
                            self.gradj_name.append(fast_param_name)
                    
                    #store for each iteration
                    all_ite_files = dict()
                    for id_ite in range(1, self.nb_ite+1):
                        tmp_files=dict()
                        # store if parameter is divised (as hydrograph)
                        for id_file in range(1,nb_file+1):
                            print(id_file, id_ite) 
                            if param_name== 'hydrograph':
                                tmp_files[id_file] =  np.loadtxt(f"{self.min_dir}/"+str(''.join([x for x in files if  f"hydrograph_{id_file:03d}.{id_ite:03d}" in x])))
                            elif param_name == "manning":
                                   tmp_files[id_file] =  np.loadtxt(f"{self.min_dir}/"+str(''.join([x for x in files if  f"manning.{id_ite:03d}" in x])))
                            else:
                                print("load min for", param_name, "to implement")
                        all_ite_files[id_ite] =  tmp_files
            
                    self.param[param_name]=all_ite_files
            
        
        def save(self, hdf5_path):
            # open hdf5
            f = h5py.File(hdf5_path, "a")
            
            min_group = f.create_group("min")
            min_group.create_dataset(f"gradj", self.gradj.shape, 'f')
            min_group.create_dataset(f"j", self.j.shape, 'f')
            min_group.create_dataset(f"all_ite", np.asanyarray(self.all_ite).shape, 'f')
            
            
            min_group["gradj"][...] = self.gradj
            min_group["j"][...] = self.j
            min_group["all_ite"][...] =  np.asanyarray(self.all_ite)
            
            for i in range(len(self.gradj_name)):
                min_group.attrs[f'name_{i}'] =  self.gradj_name[i]
                
            f.close()
            
            
            
        def plot(self):
            
            plt.plot(self.all_ite, self.j     ,'b',linewidth=1,label=r'$J_{hy}$')
            
            my_colors = ["r-", "y-", "o-", "g-"]
            for i in range(self.gradj.shape[1]) :
                print(i)
                tex = r"$|| \nabla_{" +fr"{self.gradj_name[i]}" + "} J_{hy}||$"
                tex = fr"$|| \nabla {{{self.gradj_name[i]}}} ||$"
                # name update necesary soon
                plt.plot(self.all_ite,self.gradj[:,i],my_colors[i],linewidth=0.5, label= rf"{tex}")
            #plt.plot(ite,grad_cost2,'r',linewidth=1,label=r'$||\nabla_{Q_2} J_{hy}||$')
            
            plt.xlabel('Iterations')
            plt.legend()
            #plt.xscale('log')
            plt.yscale('log')
#           
            plt.show()
            
#            #titlefig = r'Evolution of $J$ and $||\nabla J||$ in function of iterations.'
#            titlefig = r'$J_{hy}$ and $||\nabla_Q J_{hy}||$ vs iterations.' 
#            plt.title(titlefig)
#            
#            pathfile=pathOutFile+'min_cost.png'
#            plt.savefig(pathfile)
#            os.system('eog '+str(pathfile) + '&' ) 
#            
