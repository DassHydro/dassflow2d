# coding: utf-8
"""
Tools for driven studies parallel run and cmd executions
========================================================
"""
import os
import subprocess
from configuration.cfg import Config
from execution.get import get_mpi_cmd
import re


def shell_cmd(cmd):
    """
    execute the command cmd, trapping the possible errors
    @param cmd (str) the command
    """
    # abort = False
    process = subprocess.Popen(cmd, shell=True, stdin=None,
                               stdout=subprocess.PIPE)

    (stdout, stderr) = process.communicate()
    del stderr
    return stdout, process.returncode


def mpirun_cmd():
    """
    returns the mpirun command to be launched in a shell
    """
    CFGS = Config()
    python_dir = os.path.dirname(os.path.realpath(__file__))
    root_dir = os.environ.get('HOMETEL')
    cfg_name = os.environ.get('USETELCFG')
    cfg_file = os.environ.get('SYSTELCFG')
    try:
        CFGS.parse_cfg_file(cfg_file, cfg_name, root_dir, python_dir)
        CFGS.compute_execution_info()
        cfgmpi = CFGS.configs[CFGS.cfgname]['MPI']
    except RuntimeError:
        cfgmpi = {}

    if cfgmpi != {}:
        cmd = get_mpi_cmd(cfgmpi)
    else:
        try:
            shcmd = "ompi_info -V -parsable 2>/dev/null"
            shout, _ = shell_cmd(shcmd)
            openmpi_series = int(re.sub(r'^.*Open MPI v([0-9]+).*$', r'\1',
                                 str(shout)))

            if openmpi_series == 1:
                cmd = "mpirun -np <ncsize> <exename>"
            elif openmpi_series == 4:
                cmd = "mpirun --mca btl ^openib " \
                      "--mca mca_component_show_load_errors 0 " \
                      "--use-hwthread-cpus --oversubscribe " \
                      "-np <ncsize> <exename>"
            else:
                cmd = "mpirun --oversubscribe " \
                      "-np <ncsize> <exename>"
        except ValueError:
            shcmd = "mpichversion -v"
            shout, return_code = shell_cmd(shcmd)
            if return_code != 0:
                cmd = "mpirun -n <ncsize> <exename>"
            else:
                # MPICH not detected, trying intelMPI
                cmd = "mpirun -n <ncsize> <exename>"

    return cmd
