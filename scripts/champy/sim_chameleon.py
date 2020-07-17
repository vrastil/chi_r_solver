"""
'sim_chameleon.py' modules serves for running and sorting simulations
"""

import subprocess
import os
import sys
import time

import numpy as np

from .plot_chameleon import plot_generic


class Simulation(object):
    def __init__(self, run_cmd, **kwargs):
        self.run_cmd = run_cmd
        self.sim_kwargs = kwargs

    def set_kwargs(self, **kwargs):
        for arg, val in kwargs.items():
            self.sim_kwargs[arg] = val

    def set_run_cmd(self, run_cmd):
        self.run_cmd = run_cmd
    
    def run(self, stdout=subprocess.PIPE):
        cmd = self.run_cmd
        for arg, val in self.sim_kwargs.items():
            cmd += f" --{arg}={val}"
        if stdout is not None:
            print(f"Running command '{cmd}' from {os.getcwd()}")
        process  = subprocess.Popen(cmd, shell=True, stdout=stdout)
        process.wait()
        if stdout is not None:
            print(process.stdout.read().decode('utf-8'))

    def get_out_dir(self, key="out_dir"):
        return self.sim_kwargs[key]

    def get_data(self, a_file):
        a_file = os.path.join(self.get_out_dir(), a_file)
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


def run_many_sims(kwargs_dflt, kwargs_sims, stdout=subprocess.PIPE):
    # set default arguments
    sim = Simulation("../build/main.a", **kwargs_dflt)
    
    # get all combinations of parameters
    all_kwargs = get_all_kwargs(kwargs_sims)
    
    # for every combination run simulation and save results
    results_all = []
    for i, kwargs_sim in enumerate(all_kwargs):
        # print some info
        print(f"Running simulation {i+1}/{len(all_kwargs)}")

        # run simulation
        sim.set_kwargs(**kwargs_sim)
        sim.run(stdout=stdout)
        
        # save results
        results = {
            'params' : kwargs_sim,
            'potential' : sim.get_data('potential.dat'),
            'forces' : sim.get_data('forces.dat'),
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

        data_all[label] = data

    # plot
    plot_generic(data_all, **kwargs)