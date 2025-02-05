#=======================================================#
# Source librairies
#=======================================================#

import dassflow2d as df2d
import numpy as np
import os
import re
import sys
import matplotlib.pyplot as plt

#=======================================================#
# copy existing case files
#=======================================================#

# working directory (USER-DEFINED)
dassflow_dir = "/home/yts/dassflow2d-devs_yuri"   

os.chdir(dassflow_dir)

print("DassFlow directory is: ", dassflow_dir)

# Define directory where case is run 
# (its name 'bin_A' is imposed in  {dassflow_dir}/code/makefile.inc : CASEDIR='bin_A')
run_dir = f"{dassflow_dir}/code/bin_A/" 

# Define directory containing case data
case_data_dir = f"{dassflow_dir}/cases/dev_case/steady_uniform/bin_A/"

# Clean run directory
os.system(f"rm -r {run_dir}*") 

# Copy case data to runing directory
os.system(f"cp -r {case_data_dir}* {run_dir}") # Copy of case files from existing case to bin_A

# Move to code directory and clean bin directory
os.chdir( f"{dassflow_dir}/code/")
os.system("make cleanres cleanmin") # Clean forward run and minimization results

#=======================================================#
# initialise + run + save results
#=======================================================#

# input file reading (simulation settings)
df2d.wrapping.read_input(f"{run_dir}/input.txt")

# Creation of dassflowmodel object using case data: 
df2d.wrapping.m_mpi.init_mpi()
my_model = df2d.dassflowmodel(bin_dir =  run_dir, hdf5_path = f"{dassflow_dir}/code/bin_A/res/simu.hdf5" , 
                              run_type = "direct", clean = True, custom_config = None)
my_model.config.get()
my_model.init_all()
my_model.run() # run model

#=======================================================#
# calculate analytical solution
#=======================================================#

from scipy.optimize import fsolve

g = my_model.config['g']
rho = my_model.config['rho']       			  # density from the input.txt
tau_c = my_model.config['tau_c']              # yield stress from input.txt
Kn = my_model.config['K_index']               # consistency index from input.txt
n = 1/my_model.config['m_powerlaw_index']     # power law index from input.txt (inverse)

##### USER DEFINED PARAMETERS TO CALCULATE ANALYTICAL SOLUTION #####

Q = 0.01  ## flow rate from the constant hydrograph   
theta = 4       ## slope (in degrees) from the mesh
w = 0.3         ## width of the channel, from the mesh

def implicit_h_function(h, tau_c=tau_c, Kn=Kn, n=n, rho=rho, g = g, w=w, Q=Q, theta=theta):

    part1 = (rho*g*np.sin(np.deg2rad(theta))/Kn)**(1/n)
    part2 = ((rho*g*h*np.sin(np.deg2rad(theta)) - tau_c)/(rho*g*np.sin(np.deg2rad(theta))))**(1+1/n)
    part3 = 1 - (n/(2*n+1))*((rho*g*h*np.sin(np.deg2rad(theta)) - tau_c)/(rho*g*h*np.sin(np.deg2rad(theta))))
    return (n/(n+1))*part1*part2*part3*h*w - Q
    

h_solved = fsolve(implicit_h_function, 1)
print(h_solved)

#=======================================================#
# get the last .dat file from the simulation
#=======================================================#

# Directory containing the .dat files
data_dir = f"{run_dir}/res/"

# Regular expression to extract time step from filenames
file_pattern = re.compile(r"result_(\d+\.\d+E[+-]\d+)\.dat")

# Get all .dat files and extract their time steps
dat_files = [f for f in os.listdir(data_dir) if file_pattern.match(f)]

if not dat_files:
    raise FileNotFoundError("No .dat files found in the directory.")

# Sort files based on the extracted time step (convert to float for sorting)
dat_files.sort(key=lambda f: float(file_pattern.match(f).group(1)))

# Select the latest file
latest_file = dat_files[-1]
latest_filepath = os.path.join(data_dir, latest_file)

print(f"Processing latest file: {latest_file}")

# Initialize an empty list to store the h values
h_values, x_values, u_values, bathy_values = [], [], [], []

# closest coordinate from the center of the channel (USER-DEFINED)
y_monitor = 1.65000000E-01  

# Open the .dat file and process it line by line
with open(latest_filepath, 'r') as file:
    for line in file:
        # Skip lines that start with '#'
        if line.startswith('#'):
            continue
        
        # Strip any leading/trailing whitespace and split the line into columns
        columns = line.strip().split()
        
        # Ensure there are enough columns to process (i.e., it's not a blank line)
        if len(columns) >= 8:
            try:
                # Convert the relevant columns to floats
                y_value = float(columns[2])  # y is the 3rd column (index 2)
                h_value = float(columns[4])  # h is the 5th column (index 4)
                u_value = float(columns[7])  # v is the 8th column (index 7)
                bathy_value = float(columns[3])
                x_value = float(columns[1])
                # Check if y equals 0.15
                if y_value == y_monitor:
                    h_values.append(h_value)
                    x_values.append(x_value)
                    u_values.append(u_value)
                    bathy_values.append(bathy_value)
            except ValueError:
                continue

plt.figure()
plt.xlabel('x (m)')
plt.ylabel('h (m)')
plt.plot(x_values, np.asarray(h_values) ,color='blue', label='Numerical solution')
plt.hlines(h_solved, 0, 10, color='black', label='Analytical solution')
plt.show()
