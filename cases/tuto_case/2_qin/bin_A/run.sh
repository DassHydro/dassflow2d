#Directory script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Clean all ( bin, cpp, tap)
make cleanall


# Copy source on the bin directory
cp $DIR/* ./bin


# Copy the direct Makefile
cp ./souk/MakeFile/Makefile_dir_seq ./Makefile

# Compile, link and launch direct case
make
make runexe

# Clean source
make clean

# Copy the adjoint Makefile
cp ./souk/MakeFile/Makefile_inv_seq ./Makefile

# Create the obs directory 
mkdir ./bin/obs

# Copy results on the obs directory
cp ./bin/res/* ./bin/obs/

# Copy first guess hydrograph
cp ./bin/hydrograph_first_guess.txt ./bin/hydrograph.txt

# Copy inverse input.txt file
cp ./bin/input_hydro.txt ./bin/input.txt


# Remove restart.bin, clean source, generate differentiate and adjoint code, compile, link and compute cost function
rm -rf bin/restart.bin
make clean
make tap_files
make clean
rm -rf bin/restart.bin
make runexe

# Compute gradient
make clean
rm -rf bin/restart.bin
make rungrad

# Identification
make clean
rm -rf bin/restart.bin
make runmin

