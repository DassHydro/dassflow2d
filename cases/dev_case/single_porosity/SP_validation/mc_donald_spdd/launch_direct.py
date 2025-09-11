import subprocess
from multiprocessing import Pool, Manager
import os
import json
import time
from mpi4py import MPI
import pandas as pd
import numpy as np 

def run_script(script_name, result_dict, num_procs, user_params):
    try:
        user_params_str = json.dumps(user_params)
        str_user_params = {key: str(value) for key, value in user_params.items()}
        env = {**os.environ, 'USER_PARAMS': user_params_str, **str_user_params}
        env.pop('OMPI_MCA_plm_rsh_agent', None)
        command = f"xterm -hold -e python3 {script_name} {num_procs} '{user_params_str}'"
        subprocess.run(command, shell=True, check=True, env=env)
        result_dict[script_name] = "Success"
    except subprocess.CalledProcessError as e:
        result_dict[script_name] = f"Error: {e}"

if __name__ == "__main__":
    

    code_dir = os.getcwd() 

    script_names = ["direct_sim.py"]
    result_dict = Manager().dict()
    num_procs = 1

    common_params = {}

    user_params_list = [     

    {"bin_folder": f"{code_dir}"
     , "mesh_name":'channel.geo'
     , "porosity_land_file_name":f"{code_dir}/land_uses.txt"
     , "porosity_values_file_name":f"{code_dir}/porosity_values.csv"
     , "porosity_gamma_file_name":f"{code_dir}/porosity_beta.csv"
     , "porosity_hbanks_file_name":f"{code_dir}/porosity_hbanks.csv"
     , "friction_prior":0.025
     #, "ts":1296000
     ,"ts": 86400
     , "dtw":3600},
    ]

    for params in user_params_list:
        params.update(common_params)

    with Pool() as pool:
        pool.starmap(run_script, [(script, result_dict, num_procs, params) for params in user_params_list for script in script_names])

    time.sleep(0.5)

    for script, result in result_dict.items():
        print(f"Script '{script}' execution result: {result}")


