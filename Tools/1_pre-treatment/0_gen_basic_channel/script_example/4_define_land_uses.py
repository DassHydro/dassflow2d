#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 11:20:57 2022

@author: livillenave
"""


os.chdir("/home/livillenave/Documents/distant/dassflow2d-wrap/Tools/2_gen_basic_channel")


import gen_channel_case
gen_channel_case.gen_land_uses(		manning_alpha  = 0.003, 
                                   manning_beta = 0)