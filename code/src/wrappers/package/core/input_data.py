import dassflow2d
import numpy as np
import pandas as pd
import datetime
from osgeo import gdal
import matplotlib.pyplot as plt

from datetime import datetime #, timedelta
import glob
from tqdm import tqdm

#from typing import List

import warnings

#class Spatial_data (object):
    
     ##* numpy array (nb cellule) qui correspond id patch
     ##* numpy ndarray(nb_patch, 2)
    

#class Spatio_temporal_data(object):

     ## add temporal dim
     ## add machin stocker des raster et les tracer... 
     
#class Spatial_param (object)
# 


class Rain(object): 
    """
    Rain is composed of:
            self.mesh_fortran = mesh_fortran                  : msh derived type fortran
            self.path_prcp_directory = path_prcp_directory    : str
            self.dx_rain = dx_rain                            : integer
            self.dy_rain = dy_rain                            : integer
            self.dt_rain = dt_rain                            : integer
            self.start_date = start_date                      : str : date format : "%Y-%m-%d %H:%M:%S"
            self.end_date = end_date                          : str : date format : "%Y-%m-%d %H:%M:%S"
            self.rainfall_data = self.source_rainfall()       : (output) Panda data frame of rainfall data ...
    
    various plot method are available
    """
    
    def __index_containing_substring(self,the_list: list, substring: str):
        for i, s in enumerate(the_list):
            if substring in s:
                return i
        return -1


    def __read_windowed_raster(self,path: str, spatial_window) -> np.ndarray:
        
        # From User defined spatial window for reading external data map  
        xmin_window = spatial_window[0]
        ymin_window = spatial_window[1]
        xmax_window = spatial_window[2]
        ymax_window = spatial_window[3]
        dx_window   = spatial_window[4]
        dy_window   = spatial_window[5]
        
        nrow_window = int((ymax_window - ymin_window)/dy_window)+1
        ncol_window = int((xmax_window - xmin_window)/dx_window)+1

        # Raster data
        ds = gdal.Open(path)
        
        transform = ds.GetGeoTransform()

        xmin_raster = transform[0]
        ymax_raster = transform[3]
        dx_raster = transform[1]
        dy_raster = -transform[5]
        
        # Offset computation and data reading
        col_off = (xmin_window - xmin_raster) / dx_raster
        row_off = (ymax_raster - ymax_window) / dy_raster

        band = ds.GetRasterBand(1)
        
        nodata = band.GetNoDataValue()
        # print("col_off, row_off, mesh.ncol, mesh.nrow",col_off, row_off, mesh.ncol, mesh.nrow)
        arr = band.ReadAsArray(col_off, row_off, ncol_window, nrow_window)

        arr = np.where(arr == nodata, -99, arr)
        #print("arr",arr)
        return arr#, xmin, ymax, xres, yres

    def __read_rain_data(self,path_prcp_directory, mesh_data, spatial_window, time_window, prcp_format = "tiff"):

        # TODO in fortran : implement id_cell (to ensure data reading and mapping for parallel runs)
        #id_dass,xcoord,xcoord,bathy,h,zs,Manning,u,v=np.genfromtxt(path_mesh,unpack=True)

        # From User defined spatial window for reading external data map    
        xmin_window = spatial_window[0]
        ymin_window = spatial_window[1]
        xmax_window = spatial_window[2]
        ymax_window = spatial_window[3]
        dx_window   = spatial_window[4]
        dy_window   = spatial_window[5]
        
        nrow_window = int((ymax_window - ymin_window)/dy_window)+1
        ncol_window = int((xmax_window - xmin_window)/dx_window)+1

        x_cell = mesh_data[0]
        y_cell = mesh_data[1]
        id_cell = mesh_data[2]
        nxd = len(x_cell)

        # From User defined temporal window for reading external data map
        # start_date  = time_window[0]
        # end_date    = time_window[1]
        # dt_data     = time_window[2] # In seconds
        

        date_range = pd.date_range(
                                start = time_window[0],
                                end   = time_window[1],
                                freq  = f"{int(time_window[2])}s",
                                )[1:]

        ntimes = len(date_range)
        
        prcp_matrix_all_times = []

        #Prepare objects for output dataframe
        prcp_cell_dass = np.zeros((nxd,ntimes))
        # id_dasst=np.zeros((nxd))
        id_col_ras=np.zeros((nxd))
        id_row_ras=np.zeros((nxd))
        id_dass_tile=np.zeros((nxd))

        #if prcp_format == "tif":
        files = sorted(glob.glob(f"{path_prcp_directory}/**/*tif*", recursive=True))
            #files = _adjust_left_files(files, date_range)
        # other formats
        # python array from model    
        # ind = -1

        for i, date in enumerate(tqdm(date_range, desc="</> Reading precipitation")):
            date_strf = date.strftime("%Y%m%d%H%M")

            ind = self.__index_containing_substring(files, date_strf)
            prcp_matrix_current_time = []

            if ind == -1:
                prcp_matrix_current_time = [[-99] * ncol_window] * nrow_window #CHECK THIS wrt row col transpose
                prcp_matrix_all_times.append([[-99] * ncol_window] * nrow_window) #input_data.prcp[..., i] = -99.0
                warnings.warn(f"Missing precipitation file for date {date}")

            else:
                prcp_matrix_current_time = (
                    self.__read_windowed_raster(files[ind], spatial_window) #* setup.prcp_conversion_factor
                )
                
                prcp_matrix_all_times.append(prcp_matrix_current_time) #input_data.prcp[..., i] = matrix

                files.pop(ind)
                
            for j in range(0,nxd): #len(mesh1.x_cell)
                # "projection" with integer portion, i.e. for mesh1 cell center within each raster pixel
                # count += 1
                
                # id_dasst[j] = count
                
                id_col = int((x_cell[j] - xmin_window) / dx_window)
                id_row = int((ymax_window - y_cell[j]) / dy_window)

                #print("id_col",id_col,"id_row", id_row)
                
                prcp_cell_dass[j,i] = prcp_matrix_current_time[id_col][id_row]# prcp_raster[-1][id_col][id_row]
                #print("matrix(t) : ", prcp_matrix_current_time[id_col][id_row])
                
                id_col_ras[j] = id_col
                id_row_ras[j] = id_row
                
                id_dass_tile[j] = int((id_row * ncol_window) + id_col)
                
        colnames=["x_cell_dass","y_cell_dass", "id_dass", "id_row_raster","id_col_raster","id_dass_tile"]
        
        for i in range(1,ntimes+1):
            dtname = f"P_t_{i}"
            colnames.append(dtname)
            
        alldata = np.column_stack((x_cell,y_cell,id_cell,id_row_ras,id_col_ras,id_dass_tile,prcp_cell_dass))
        
        df_rain = pd.DataFrame(alldata, index = range(0,nxd), columns=colnames)
        
        return prcp_matrix_all_times, df_rain, ncol_window, nrow_window, ntimes
    
    
    def __init__(self, mesh_fortran, path_prcp_directory, dx_rain, dy_rain, dt_rain, start_date, end_date) :

            self.mesh_fortran = mesh_fortran
            self.path_prcp_directory = path_prcp_directory
            self.dx_rain = dx_rain
            self.dy_rain = dy_rain
            self.dt_rain = dt_rain
            self.start_date = start_date
            self.end_date = end_date
            self.rainfall_data = self.source_rainfall()
    
    def __get_mesh_centers_from_mesh_fortran(self, mesh_fortran):
        """
        Build an array of cell centers coordinates and cell_id from mesh
        
        Input : fortran mesh
        
        Output : array(x_center[1:nc], y_center[1:nc], id_cell[1:nc], x_min, y_min, x_max, y_max) with nc the number of mesh cells 
        """
        my_mesh = mesh_fortran
        
        x_center =[]
        y_center = []
        id_cell  = []
        
        for i in range(0,my_mesh.nc):
            x_center.append(my_mesh.cell[i].grav.x)
            y_center.append(my_mesh.cell[i].grav.y) 
            id_cell.append(i) # TOBE ADDED into the fortran (for parallel safety)
        
        # Mesh "Extent" : outter rectangle
        x_min = min(x_center)
        x_max = max(x_center)
        y_min = min(y_center)
        y_max = max(y_center)
        
        return(x_center, y_center, id_cell, x_min, y_min, x_max, y_max)

    
    def source_rainfall(self):
        
        my_mesh = self.mesh_fortran
        path_prcp_directory = self.path_prcp_directory
        
        date_format = "%Y-%m-%d %H:%M:%S"
        start_date = datetime.strptime(self.start_date,date_format) # start_date = datetime.strptime(raw_segm_events.iloc[run_id,2], date_format)
        end_date = datetime.strptime(self.end_date,date_format)   # end_date = datetime.strptime(raw_segm_events.iloc[run_id,3], date_format)
        
        (x_cell, y_cell, id_cell, xmin, ymin, xmax, ymax) = self.__get_mesh_centers_from_mesh_fortran(my_mesh)
        
        # spatial data
        mesh_data = [x_cell, y_cell, id_cell]
        time_window = [start_date,end_date,self.dt_rain]
        spatial_window = [xmin-self.dx_rain, ymin-self.dy_rain, xmax+self.dx_rain, ymax+self.dy_rain, self.dx_rain, self.dy_rain] #enlarge extent of one rain pixel to ensure covering the mesh

        ########################################################################################################################################
        #read data maps
        ########################################################################################################################################
        # read rainfall maps (spatio-temporal, evenly spaced)
        prcp = []
        prcp_raster, df_rain, ncol_rain, nrow_rain, ntimes_rain = self.__read_rain_data(self.path_prcp_directory, mesh_data, spatial_window, time_window, prcp_format = "tiff")
        
        return df_rain
    
    
class Input_data(object):
    """
    Input_data is composed of:
        
        -  df_rain, a panda data frame containing : "x_cell_dass","y_cell_dass", "id_dass", "id_row_raster","id_col_raster","id_dass_tile, f"P_t_{i}" "
        -
    
    various plot method are available
    """
    
    def __init__(self, path_prcp_directory, mesh_fortran, dassflow_dir, ) :
		
            self.mesh_fortran = mesh_fortran
            self.dassflow_dir = dassflow_dir
            self.rain = Rain(mesh_fortran        = mesh_fortran, 
                             path_prcp_directory = path_prcp_directory, 
                             dx_rain             = dx_rain, 
                             dy_rain             = dy_rain, 
                             dt_rain             = dt_rain, 
                             start_date          = start_date,  
                             end_date            = end_date)    
    
    

