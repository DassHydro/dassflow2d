#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 10:36:00 2022

@author: livillenave
"""

import dassflow2d as df2d
#import numpy as np
#import sys
#import os

#=======================================================#
#  run fortran model
#=======================================================#
bin_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A"
# initialise fortran instance, and python corrponding data
my_model = df2d.DassFlowModel(bin_dir = bin_dir, run_type = "direct") # run_type can be min or direct (grad ?)