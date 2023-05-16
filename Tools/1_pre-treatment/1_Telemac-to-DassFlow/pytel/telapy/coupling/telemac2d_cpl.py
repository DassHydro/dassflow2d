# -*- coding: utf-8 -*-
"""
main program of coupling 2D model
"""
import sys
import time
from os import path, environ

from telapy.coupling.class_mod_2d import ClassMod2D

if (not path.exists(path.join(environ.get('HOMETEL', ''),
                              'builds',
                              environ.get('USETELCFG', ''),
                              'wrap_api'))):
    print("  -> TelApy not available doing nothing")
    sys.exit(0)


def telemac2d_cpl():
    """
    Main function of script
    @return nothing
    """
    # Initialize an instance of a 2D coupled model
    start_time = time.time()
    cpl = ClassMod2D()
    if cpl.coupler.rank == 0:
        print("\n--- Initialisation {} seconds ---\n".format((time.time()
              - start_time)))

    # Loop on the time integration
    for cpl.stp in range(cpl.cplsteps):
        start_uloop_time = time.time()

        # Loop on the coupling updates
        for cpl.iter in range(cpl.maxiter):

            # Integrate the model for the current iterate
            cpl.run()

            # Exhchange the boundary conditions with the 1D models
            cpl.transmit_step()

            # Check convergence
            if cpl.converged:
                break

        # Output at the given frequency
        cpl.write_output()
        if cpl.coupler.rank == 0:
            print("\n--- {} s, Time loop {} s ---\n"
                  .format(cpl.time, time.time() - start_uloop_time))

    # Finalization
    cpl.finalize()
    if cpl.coupler.rank == 0:
        print("\n--- {} seconds ---\n".format((time.time() - start_time)))


if __name__ == "__main__":
    telemac2d_cpl()
