# -*- coding: utf-8 -*-
"""
main program of coupling 1D models
"""
from telapy.coupling.class_mod_1d import ClassMod1D


def mascaret_cpl():
    """
    Parallel 1D instances of Mascaret in a coupling

    """
    # Initialize an instance of a 1D coupled model
    mod = ClassMod1D()
    # Loop on the time integration
    for mod.cplstp in range(mod.cplsteps):

        # Initialise a coupling step
        mod.cpl.init_cplstp()
        # Loop on the coupling updates
        for mod.iter in range(mod.maxiter):

            # Update the boundary conditions at the coupling sections
            mod.cpl.update_bc()

            # Integrate the model for the current iterate
            mod()

            # Trasmit the boundary conditions to the 2D model
            mod.cpl.transmit_bc()

            # Check convergence
            mod.cpl.check_convergence()
            if mod.converged:
                break

        # Output at the given frequency
        mod.selective_output()

    # Finalization
    mod.finalize()
    del mod


if __name__ == "__main__":
    mascaret_cpl()
