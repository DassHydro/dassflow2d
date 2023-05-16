# coding: utf-8
"""
Telemac2DStudy and Telemac2DTestCase classes
============================================
"""
import io
import os
import stat
import shutil
import time

import numpy as np

from telapy.api.t2d import Telemac2d
from data_manip.extraction.telemac_file import TelemacFile
from postel.plot2d import *
from telapy.tools.driven_utils import shell_cmd, mpirun_cmd
import matplotlib
import matplotlib.animation as manimation
import matplotlib.pyplot as plt
matplotlib.use('Agg')


class Telemac2DTestCase(object):
    """ Class defining a test case for a study with Telemac2D. """

    def __init__(self, steering_file, user_fortran, ks_area="ks_area.txt",
                 path=".", name=None):
        """ Constructor

        @param steering_file (str) steering file
        @param user_fortran (str) user Fortran file
        @param name (str) name of the test case
        @param ks_area (str) Name of file containing strickler area
        @param path (str) path of the test case
        """
        # Configuration files
        self.steering_file = os.path.realpath(steering_file)
        self.user_fortran = os.path.realpath(user_fortran)
        self.ks_area = os.path.realpath(ks_area)
        if os.path.isfile(path):
            self.path = os.path.dirname(os.path.realpath(path))
        if os.path.isdir(path):
            self.path = os.path.realpath(path)
        self.name = name

    def __call__(self, t2d, x_val):
        """ Change configuration of a Telemac2d instance

        @param t2d (Telemac2d) instance of Telemac2d
        @param x_val (list) inputs
        """
        pass


class Telemac2DStudy(object):
    """Telemac2DStudy API."""

    def __init__(self, points, test_case, results_file, work_dir=".", stdout=6,
                 nproc=1, mesh_state=False):
        """ Constructor
        @param points (list) points of interest [[x1, y1], ..., [xN, yN]]
        @param test_case (Telemac2DTestCase) test case for the study
        @param results_file (str) result file (SELAFIN)
        @param work_dir (str) working directory with files specified
         by test_case
        @param stdout (int) standard output (defaut = 6 [console];
                           if 666 => file 'fort.666')
        @param nproc (int) number of processors
        @param mesh_state (bool) store the final hydraulic state at each node
        """
        # Working directory
        self.work_dir = os.path.realpath(work_dir)
        if not os.path.exists(self.work_dir):
            os.makedirs(self.work_dir)

        # Enter working directory
        caller_path = os.getcwd()
        os.chdir(self.work_dir)
        # Number of processors
        self.nproc = nproc
        # Test case
        self.test_case = test_case
        # Result_file
        self.results_file = results_file
        # Coordinates of interest
        self.points = points
        # Standard output
        self.stdout = stdout
        # Initialize 'output' at None
        self.output = None
        # Store all mesh state
        self.mesh_state = mesh_state
        # Back to caller directory
        os.chdir(caller_path)

    def __call__(self, x_val=None, finalize=False, out_dir="."):
        """ Run Telemac2D

        @param x_val (list) inputs [Ks1, Ks2, Ks3, Ks4, CDZ, Q2, Q3]
        @param finalize (bool) delete the Telemac 2D instance after execution
        @param out_dir (str) relative path of the directory of results
        """

        # Enter working directory
        caller_path = os.getcwd()
        os.chdir(self.work_dir)

        out_dir = os.path.realpath(out_dir)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        if x_val is not None:
            with open(os.path.join(out_dir, 'inputs.txt'), 'w') as thefile:
                thefile.write("%s" % x_val)

        self.run_simulation(x_val)

        with open(os.path.join(out_dir, 'outputs.txt'), 'w') as thefile:
            thefile.write("%s" % self.z_node)

        shutil.move(self.results_file,
                    os.path.join(out_dir, self.results_file))
        if self.stdout == 666:
            shutil.move("fort.666", os.path.join(out_dir, "log"))

        # Back to caller directory
        os.chdir(caller_path)

    def run_simulation(self, x_val, filename="run_launcher.py"):
        """ Run Telemac2D Simulation

        @param x_val (list) inputs [Ks1, Ks2, Ks3, Ks4, CDZ, Q2, Q3]
        @param filename (str) Python file for the simulation
        """
        t_0 = time.time()
        self.create_launcher_file(x_val, filename)
        filename = os.path.join(self.work_dir, filename)

        if self.nproc == 1:
            print("Sequential Run")
            shell_cmd(filename)
        else:
            print("Parallel Run on {} procs".format(self.nproc))
            self.mpirun(filename)

        self.get_water_depth()
        if self.mesh_state:
            self.get_mesh_water_depth()

        t_1 = time.time()
        self.time = t_1 - t_0

    def create_launcher_file(self, x_val, filename):
        """ Create the Python file for the execution

        @param x_val (list) inputs [Ks1, Ks2, Ks3, Ks4, CDZ, Q2, Q3]
        @param filename (str) Python file for the MPI execution
        """
        with io.FileIO(filename, "w") as file:
            file.write(self.cmd2str("header"))
            file.write(self.cmd2str("newline"))
            file.write(self.cmd2str("commworld"))
            file.write(self.cmd2str("newline"))
            file.write(self.cmd2str("create"))
            if self.nproc > 1:
                file.write(self.cmd2str("barrier"))
            file.write(self.cmd2str("newline"))
            file.write(self.cmd2str("setcase"))
            file.write(self.cmd2str("newline"))
            file.write(self.cmd2str("resultsfile"))
            file.write(self.cmd2str("newline"))
            file.write(self.cmd2str("init"))
            file.write(self.cmd2str("newline"))
            file.write(self.cmd2str("testcase"))
            file.write(self.cmd2str("newline"))
            file.write(self.cmd2str("setX", x_val))
            if self.nproc > 1:
                file.write(self.cmd2str("barrier"))
            file.write(self.cmd2str("newline"))
            file.write(self.cmd2str("run"))
            file.write(self.cmd2str("newline"))
            if self.nproc > 1:
                file.write(self.cmd2str("barrier"))
            file.write(self.cmd2str("newline"))
            file.write(self.cmd2str("finalize"))
            file.write(self.cmd2str("newline"))
            file.write(self.cmd2str("del"))
        file.close()
        os.chmod(filename, os.stat(filename).st_mode | stat.S_IEXEC)

    def reset(self, l_t=0., a_t=0.):
        """ Reset the provider

        @param l_t (float) current time step
        @param a_t (float) current time
        """
        self.provider.set_state(self.initial_state[0][:],
                                self.initial_state[1][:],
                                self.initial_state[2][:])
        self.provider.set("MODEL.LT", l_t)
        self.provider.set("MODEL.AT", a_t)
        self.provider.set("MODEL.COMPLEO", l_t)
        self.provider.set("MODEL.NTIMESTEPS", self.ntimesteps)

    def get_water_depth(self):
        """ Extract the water depth at the points of interest."""
        self.output_node = []
        self.output_points = []
        results = TelemacFile(self.results_file)
        if "HAUTEUR D'EAU" in results.varnames:
            varname = "HAUTEUR D'EAU"
        if "WATER DEPTH" in results.varnames:
            varname = "WATER DEPTH"
        t_val = results.ntimestep - 1
        self.output_node = [results.get_data_value(varname, t_val)
                            [results.get_closest_node(point)]
                            for point in self.points]
        self.output_points = results.get_data_on_points(varname, t_val,
                                                        self.points)
        del results

    def get_mesh_water_depth(self):
        """ Extract the water depth at each node."""

        results = TelemacFile(self.results_file)
        if "HAUTEUR D'EAU" in results.varnames:
            varname = "HAUTEUR D'EAU"
        if "WATER DEPTH" in results.varnames:
            varname = "WATER DEPTH"
        t_val = results.ntimestep - 1

        self.mesh_water_depth = {'x-coord': results.meshx,
                                 'y-coord': results.meshy,
                                 'waterDepth':
                                 results.get_data_value(varname, t_val)}
        del results

    def mpirun(self, filename):
        """ Launch the Python script 'filename' in parallel

        @param filename (str) Python file for the MPI execution
        """
        cmd = mpirun_cmd()

        cmd = cmd.replace('<ncsize>', str(self.nproc))
        cmd = cmd.replace('<exename>', filename)
        cmd = cmd + ' 1> mpi.out 2> mpi.err'

        _, return_code = shell_cmd(cmd)
        if return_code != 0:
            raise Exception("\nERROR IN PARALLEL RUN.\n"
                            "THE COMMAND IS : {} \n"
                            " PROGRAM STOP.\n".format(cmd))

    def cmd2str(self, keyword, x_val=None):
        """ Convert a keyword into Python lines for writing the Python script
        used by the function 'self.mpirun(filename)'

        @param keyword (str) keyword to convert into Python lines
        @param x_val (list) inputs [Ks1, Ks2, Ks3, Ks4, CDZ, Q2, Q3]
        """
        if keyword == "header":
            string = ("#!/usr/bin/env python3\n"
                      "# Class Telemac2d import\n"
                      "import sys\n"
                      "sys.path.append('"+self.test_case.path+"')\n"
                      "from telapy.api.t2d import Telemac2d\n"
                      "from testcases_t2d_driven import " +
                      self.test_case.name + "\n")
        elif keyword == "commworld":
            string = ("try:\n" +
                      "    from mpi4py import MPI\n" +
                      "    comm = MPI.COMM_WORLD\n" +
                      "except:\n" +
                      "    comm = None\n")
        elif keyword == "testcase":
            string = ("test_case = " + self.test_case.name + "('" +
                      self.test_case.steering_file + "', " +
                      "user_fortran='" + self.test_case.user_fortran + "', " +
                      "ks_area='" + self.test_case.ks_area + "', " +
                      "name='" + self.test_case.name + "', " +
                      "config=" + str(self.test_case.config) + ")\n")
        elif keyword == "create":
            string = ("t2d = Telemac2d('" + self.test_case.steering_file +
                      "', " +
                      "user_fortran='" + self.test_case.user_fortran + "', " +
                      "comm=comm, " +
                      "stdout=" + str(self.stdout) + ")\n")
        elif keyword == "barrier":
            string = "comm.Barrier()\n"
        elif keyword == "setcase":
            string = "t2d.set_case()\n"
        elif keyword == "init":
            string = "t2d.init_state_default()\n"
        elif keyword == "setX":
            string = "test_case(t2d, " + str(x_val) + ")\n"
        elif keyword == "run":
            string = "t2d.run_all_time_steps()\n"
        elif keyword == "finalize":
            string = "t2d.finalize()\n"
        elif keyword == "del":
            string = "del(t2d)\n"
        elif keyword == "resultsfile":
            string = "t2d.set('MODEL.RESULTFILE', '" + \
                self.results_file + "')\n"
        elif keyword == "newline":
            string = "\n"
        return string.encode()

    def finalize(self):
        """ Delete the Telemac 2D instance"""
        if self.nproc >= 1:
            self.provider.finalize()
            del self.provider

    @property
    def z_node(self):
        """ Get the output at final time and at points of interest

        Use the method "nearest node"
        """
        return self.output_node

    @property
    def z_points(self):
        """ Get the output at final time and at points of interest

        Use the method "Triangulation from element"
        """
        return self.output_points

    @property
    def npoints(self):
        """ Get the number of points of interest """
        return len(self.points)

    def plot_final_state(self, out_dir=".", filename="water_level.pdf",
                         plot_mesh=False):
        """
        Plot the hydraulic state at final time from the .slf output file into
        an output PDF file.

        @param out_dir (str) Output directory
        @param filename (str) output PDF file
        @param plot_mesh (boolean) Display the mesh
        """
        # Enter working directory
        caller_path = os.getcwd()
        os.chdir(self.work_dir)

        out_dir = os.path.realpath(out_dir)
        results = TelemacFile(os.path.join(out_dir, self.results_file))
        final_time = results.ntimestep - 1
        if "HAUTEUR D'EAU" in results.varnames:
            varname = "HAUTEUR D'EAU"
        if "WATER DEPTH" in results.varnames:
            varname = "WATER DEPTH"
        values = results.get_data_value(varname, final_time)

        fig, ax = plt.subplots(1, 1, figsize=(8, 5.8), dpi=300)

        if plot_mesh:
            # Define wet and dry zones mask
            mask_wet = mask_triangles(results.tri, values, relation='geq',
                                      threshold=0.01)
            mask_dry = mask_triangles(results.tri, values, relation='leq',
                                      threshold=0.01)

            # Plotting mesh only on dry zones
            results.tri.set_mask(mask_wet)
            plot2d_triangle_mesh(ax, results.tri, color='k', linewidth=0.01)
            # Plotting scalar map only on wet zones
            results.tri.set_mask(mask_dry)

        plot2d_scalar_map(fig, ax, results.tri, values,
                          shading='gouraud', cmap_name='winter',
                          data_name='water depth (m)')

        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')

        fig.savefig(filename, bbox_inches='tight')
        plt.close('all')
        del results

        # Back to caller directory
        os.chdir(caller_path)

    def film_state(self, fps=15, out_dir=".", filename="water_level.mp4",
                   plot_mesh=False):
        """
        film_state function
        @param fps (int) frame per second
        @param out_dir (str) Output directory
        @param filename (str) output file name
        @param plot_mesh (boolean) Display the mesh
        """
        matplotlib.use("Agg")
        ffmpeg_writer = manimation.writers['ffmpeg']
        metadata = dict(title='Movie Test', artist='Matplotlib',
                        comment='Movie support!')
        writer = ffmpeg_writer(fps=fps, metadata=metadata)

        results = TelemacFile(os.path.join(out_dir, self.results_file))
        ntimes = results.ntimestep
        if "HAUTEUR D'EAU" in results.varnames:
            varname = "HAUTEUR D'EAU"
        if "WATER DEPTH" in results.varnames:
            varname = "WATER DEPTH"
        fig, ax = plt.subplots(1, 1, figsize=(8, 5.8), dpi=300)

        with writer.saving(fig, filename, ntimes):
            for idt in range(ntimes):
                values = results.get_data_value(varname, idt)
                if plot_mesh:
                    # Define wet and dry zones mask
                    mask_wet = mask_triangles(results.tri, values,
                                              relation='geq', threshold=0.01)
                    mask_dry = mask_triangles(results.tri, values,
                                              relation='leq', threshold=0.01)

                    # Plotting mesh only on dry zones
                    results.tri.set_mask(mask_wet)
                    plot2d_triangle_mesh(ax, results.tri, color='k',
                                         linewidth=0.01)
                    # Plotting scalar map only on wet zones
                    results.tri.set_mask(mask_dry)

                plot2d_scalar_map(fig, ax, results.tri, values,
                                  shading='gouraud', cmap_name='winter',
                                  data_name='water depth (m)')

                ax.set_xlabel('X (m)')
                ax.set_ylabel('Y (m)')
                plt.show()
                writer.grab_frame()

        plt.close(fig)
        del results
