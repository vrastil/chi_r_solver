import subprocess
import os
import sys


class Simulation(object):
    def __init__(self, run_cmd, **kwargs):
        self.run_cmd = run_cmd
        self.sim_kwargs = kwargs

    def set_kwargs(self, **kwargs):
        for arg, val in kwargs.items():
            self.sim_kwargs[arg] = val

    def set_run_cmd(self, run_cmd):
        self.run_cmd = run_cmd
    
    def run(self):
        cmd = self.run_cmd
        for arg, val in self.sim_kwargs.items():
            cmd += f" --{arg}={val}"
        print(f"Running command '{cmd}'")
        process  = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        print(process.stdout.read().decode('utf-8'))

