#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 14:59:59 2023

@author: livillenave
"""


tmp = np.ndarray(shape = (len(all_patches[0]), len(all_patches)+1), dtype = 'int')
tmp[:,0] = np.arange(len(all_patches[0]))+1

for i in range(len(all_patches)):
        tmp[:,i+1] = all_patches[i][:]



reorder=  np.argsort(tmp[:,4])
tmp = tmp[reorder,:]


# mapping_patch1_to_patch2 = 

values_new_patch = np.unique(tmp[:,2])



X = np.zeros(shape = (1))
Y = np.zeros(shape = (10))
P = np.ones(shape = (len(Y),len(X)) )

X[:] = 0.05
Y = P @ X

tmp[:, (2,3)]
b = tmp[:, (2,3)]

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

previous_patch = 2
new_patch= 3

X = np.zeros(shape = (10))
X[:] = 0.05 + np.random.random(10)

MATRIX_P = build_p_matrix(tmp, 2, 3)
Y = MATRIX_P @ X




