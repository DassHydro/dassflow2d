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

dassflow_dir = "/home/yts/dassflow2d-devs_yuri"
os.chdir(dassflow_dir)
print("DassFlow directory is: ", dassflow_dir)

# Define directory where case is run 
# (its name 'bin_A' is imposed in  {dassflow_dir}/code/makefile.inc : CASEDIR='bin_A')
run_dir = f"{dassflow_dir}/code/bin_A/" 

# Define directory containing case data
case_data_dir = f"{dassflow_dir}/cases/dev_case/roll_waves/bin_A/"

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
# start data filtering for roll waves in spatial domain
#=======================================================#

# Directory containing the .dat files
data_dir = f"{run_dir}/res/"

# Regular expression to extract time step from filenames
file_pattern = re.compile(r"result_(\d+\.\d+E[+-]\d+)\.dat")

# Get all .dat files and extract their time steps
dat_files = [
    f for f in os.listdir(data_dir) if file_pattern.match(f)
]

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

# the values below vary with the mesh

x_monitor = 2.00192300E+000   ## fixed point to evaluate RW in time domain
y_monitor = 1.55000000E-001   ## the closest to the center of the channel

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
                # Handle any conversion errors gracefully
                continue

## plot data in spatial domain ##

plt.figure()
plt.title('Roll Waves in space domain')
plt.xlabel('x (m)')
plt.ylabel('h (m)')
plt.plot(x_values, np.asarray(h_values) ,color='blue')
plt.show()

#=======================================================#
# start data filtering for roll waves in time domain
#=======================================================#

folder_path = data_dir

# Target coordinates
target_x = x_monitor
target_y = y_monitor

# Arrays to store timesteps and corresponding h values
timesteps = []
h_values = []
u_values = []

# Iterate through files in the folder
for filename in os.listdir(folder_path):
    if filename.endswith(".dat"):
        print(filename)
        # Extract timestep from filename
        match = re.search(r"result_(\d+\.\d+E[+-]\d+)\.dat", filename)
        if match:
            timestep = float(match.group(1))
            timesteps.append(timestep)
            
            # Open and read the file
            with open(os.path.join(folder_path, filename), "r") as file:
                # Skip the first two rows
                next(file)
                next(file)
                
                for line in file:
                    # Split line into columns
                    cols = line.split()
                    x, y, h, u = float(cols[1]), float(cols[2]), float(cols[4]), float(cols[7])
                    
                    # Check if x and y match the target coordinates
                    if abs(x - target_x) < 1e-8 and abs(y - target_y) < 1e-8:
                        h_values.append(h)
                        u_values.append(u)
                        break

# Ensure the arrays are sorted by timestep
sorted_indices = np.argsort(timesteps)
timesteps = np.array(timesteps)[sorted_indices]
h_values = np.array(h_values)[sorted_indices]
u_values = np.array(u_values)[sorted_indices]

## export reference results ##

TVD_MC = np.loadtxt(run_dir + 'cunha2024b_TVDMC_time.txt')
experimental = np.loadtxt(run_dir + 'cunha2024b_experimental.txt')

## plot data in time domain, for a fixed point of x = 2m ##

## obs.: Roll Waves from DassFlow and TVD-MacCormack are moved in time to match the peaks.

plt.figure()
plt.title('Roll Waves in time domain, for $x = 2$ m')
plt.xlim(0,3)
plt.plot(timesteps - 4.053, np.asarray(h_values), color='blue', label='DassFlow')
plt.plot(TVD_MC[:,0] - 5.760153e+00, TVD_MC[:,1], color='red', label='Numerical (TVD-MacCormack)')
plt.plot(experimental[:,0], experimental[:,1], color='black', label='Experimental')
plt.xlabel('t (s)')
plt.ylabel('h (m)')
plt.legend(loc='upper right')
plt.show()
