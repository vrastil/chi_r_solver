"""
'sim_chameleon.py' modules serves for running and sorting simulations
"""

import subprocess
import os
import sys
import time
import copy

import numpy as np

from .plot_chameleon import plot_generic

NOT_KWARGS = [
    'label'
]

class Simulation(object):
    def __init__(self, run_cmd, **kwargs):
        self.run_cmd = run_cmd
        self.sim_kwargs = kwargs
        self.default_dir = None

    def set_kwargs(self, **kwargs):
        for arg, val in kwargs.items():
            if arg not in NOT_KWARGS:
                self.sim_kwargs[arg] = val

    def set_run_cmd(self, run_cmd):
        self.run_cmd = run_cmd

    def get_parallel_subdirectory(self, i):
        return f"run_{i}"

    def set_parallel_dir(self, i, key="out_dir"):
        if self.default_dir is None:
            self.default_dir = self.sim_kwargs[key]
        self.sim_kwargs[key] = os.path.join(self.default_dir, self.get_parallel_subdirectory(i), "")

    def set_parallel(self, i,  key="out_dir"):
        self.set_parallel_dir(i, key=key)
    
    def run(self, stdout=subprocess.PIPE, parallel=False):
        cmd = self.run_cmd
        for arg, val in self.sim_kwargs.items():
            cmd += f" --{arg}={val}"
        if stdout is not None:
            print(f"Running command '{cmd}' from {os.getcwd()}")
        process  = subprocess.Popen(cmd, shell=True, stdout=stdout)
        if not parallel:
            process.wait()
        if stdout is not None:
            print(process.stdout.read().decode('utf-8'))
        return process

    def get_out_dir(self, i, key="out_dir"):
        return os.path.join(self.default_dir, self.get_parallel_subdirectory(i))

    def get_data(self, a_file, i):
        a_file = os.path.join(self.get_out_dir(i), a_file)
        return np.loadtxt(a_file).transpose()

def get_all_kwargs(kwargs_sims, all_kwargs=None):
    if all_kwargs is None:
        all_kwargs = []

    __kwargs_sim = kwargs_sims.copy()
    for key, val_list in kwargs_sims.items():
        if isinstance(val_list, list):
            for val in val_list:
                __kwargs_sim[key] = val
                get_all_kwargs(__kwargs_sim, all_kwargs)
            break
    else:
        all_kwargs.append(__kwargs_sim)
    return all_kwargs


def run_many_sims(kwargs_dflt, kwargs_sims, stdout=subprocess.PIPE, parallel=False):
    # set default arguments
    sim = Simulation("../build/main.a", **kwargs_dflt)
    
    # get all combinations of parameters
    if isinstance(kwargs_sims, dict):
        all_kwargs = get_all_kwargs(kwargs_sims)
    elif isinstance(kwargs_sims, list):
        all_kwargs = copy.deepcopy(kwargs_sims)
    else:
        print(f"ERROR! Type of kwargs_sims cannot be {type(kwargs_sims)}")
    
    # for every combination run simulation and save results
    results_all = []
    processes_all = []
    for i, kwargs_sim in enumerate(all_kwargs):
        # print some info
        print(f"Running simulation {i+1}/{len(all_kwargs)}")

        # run simulation
        sim.set_kwargs(**kwargs_sim)
        sim.set_parallel(i)
        # run will wait if not parallel
        processes_all.append(sim.run(stdout=stdout, parallel=parallel))

    # wait, ignore exist codes
    if parallel:
        msg = "Waiting for subprocesses to finish"
        print(msg, end="\r")
        for i, p in enumerate(processes_all):
            p.wait()
            print(msg + f" ({i+1}/{len(processes_all)})", end="\r")

    print("\nSaving results") 
    # save results
    for i, kwargs_sim in enumerate(all_kwargs):
        results = {
            'params' : kwargs_sim,
            'potential' : sim.get_data('potential.dat', i),
            'forces' : sim.get_data('forces.dat', i),
        }
        
        results_all.append(results)
    return results_all

def find_simulations(results_all, **kwargs):
    results = []
    for sim in results_all:
        for key, val in kwargs.items():
            if key in sim['params'] and sim['params'][key] != val:
                break
        else:
            results.append(sim)
    return results

def plot_simulations(results_all, p_type='potential', label_by='Ys', data_col=1, **kwargs):
    # filter
    sims_to_plot = find_simulations(results_all, **kwargs)

    # extract data, set labels
    data_all = {}
    for sim in sims_to_plot:
        data = sim[p_type][0], sim[p_type][data_col]
        if isinstance(label_by, str):
            label = f"{label_by} = {sim['params'][label_by]}"
        elif isinstance(label_by, list):
            label = ", ".join([f"{label_by__} = {sim['params'][label_by__]}" for label_by__ in label_by])
        data_all[label] = data

    # plot
    plot_generic(data_all, **kwargs)