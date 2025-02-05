Open the file "steady_uniform_test_case.py";

Update the directory "dassflow_dir";

In the terminal, execute "steady_uniform_test_case.py" with python3;

Note that the non-Newtonian variables are all included within the file "input.txt" inside bin_A. 

Analytical solution will be calculated according to the physical characteristics of the fluid (rho, tau_c, K_index and m_powerlaw_index),
flow rate Q and channel characteristics (slope and width).

Please note that there are 3 user-defined parameters (Q, theta and w), used to calculate the analytical solution.
They must be coherent with the hydrograph.txt file (for Q) and the mesh file (theta and w).
