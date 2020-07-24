import numpy as np

from .sim_chameleon import run_many_sims, find_simulations
from .plot_chameleon import plot_generic

kwargs_dflt = {
    "err": 1E-12,
    "n": 0.5,
    "print_par" : 0,
    "out_dir" : "../output/",
    "step": 0.01,
}


def run_plot_star_pot(out_file='starlike.png', parallel=True):
    print(f"Running STAR-like simulation (potential)")

    kwargs_sims = {
        "mod" : 0,
        "Omega_m" : 1e6,
        "M200_sun" : 3E-16,
        "Ys" : [1E-13, 1E-14, 1E-15],
        "R": 1E-1,
    }

    results_all = run_many_sims(kwargs_dflt, kwargs_sims, stdout=None, parallel=parallel)
    print(f"Plot and save.")

    # extract data, set labels
    data_all = {}
    for sim in results_all:
        data = sim['potential'][0], sim['potential'][2]
        Ys = int(np.log10(sim['params']['Ys']))
        label = f"$\Phi_{{scr}} = 10^{{{Ys}}}$"
        data_all[label] = data

    plot_generic(data_all, ymin=1E-16, ymax=1E-13, xlabel=u'$r/R_s$', ylabel=u'$\tilde\chi$',
                 out_file=out_file)

def run_plot_star_for(out_file='starlike_forces.png', parallel=True):
    print(f"Running STAR-like simulation (forces)")

    kwargs_sims = {
        "mod" : 0,
        "Omega_m" : 1e6,
        "M200_sun" : 3E-16,
        "Ys" : 0,
        "R": 1,
        "R_eq" : [-1, 0.2, 0.5, 0.8, 1, 2]
    }

    results_all = run_many_sims(kwargs_dflt, kwargs_sims, stdout=None, parallel=parallel)
    print(f"Plot and save.")

    # extract data, set labels
    data_all = {}
    for sim in results_all:
        data = sim['forces'][0], sim['forces'][1]
        R_eq = sim['params']['R_eq']
        label = f"$R_{{eq}} = {R_eq} R_s$"
        data_all[label] = data


    plot_generic(data_all, yscale='linear', ymin=0, ymax=1.05, xlabel=u'$r/R_s$',
                 ylabel=u'$F_\chi/(2\\beta^2F_N)$', out_file=out_file)


def run_plot_nfw_for(out_file='nfwlike_forces.png', parallel=True):
    print(f"Running NFW-like simulation (forces)")

    kwargs_sims = {
        "mod" : 1,
        "Omega_m" : 1,
        "c": 16,
        "M200_sun": 1E-5,
        "R_eq" : [-1, 0.1, 0.5, 1, 10],
        "Ys" : 0,
    }

    results_all = run_many_sims(kwargs_dflt, kwargs_sims, stdout=None, parallel=parallel)
    print(f"Plot and save.")

    # extract data, set labels
    data_all = {}
    for sim in results_all:
        data = sim['forces'][0], sim['forces'][1]
        R_eq = sim['params']['R_eq']
        label = f"$R_{{eq}} = {R_eq} R_s$"
        data_all[label] = data


    plot_generic(data_all, yscale='linear', ymin=0, ymax=1.05, xmax=1E2, xlabel=u'$r/R_s$',
                 ylabel=u'$F_\chi/(2\\beta^2F_N)$', out_file=out_file)

def run_plot_nfw_pot_eff( out_file='nfwlike_pot_eff.png', parallel=True):
    print(f"Running NFW-like simulation (forces)")

    kwargs_sims = {
    "mod" : 1,
    "n": [0.1, 0.5, 0.7],
    "Omega_m" : 1,
    "c": 4,
    "M200_sun": 1E2,
    "R_eq" : [-100, -1, 1],
    "Ys" : 0,
    }

    results_all = run_many_sims(kwargs_dflt, kwargs_sims, stdout=None, parallel=parallel)
    print(f"Plot and save.")

    # extract data, set labels
    sims_to_plot = find_simulations(results_all, n=0.5)
    data_all = {}
    for sim in sims_to_plot:
        data = sim['potential'][0], sim['potential'][3]
        R_eq = sim['params']['R_eq']
        label = f"$R_{{eq}} = {R_eq} R_s$"
        data_all[label] = data


    plot_generic(data_all, yscale='linear', ymin=0, xmin=1E-1, xmax=1E2, xlabel=u'$r/R_s$',
                 ylabel=u'$\Phi_{scr,eff}/\Phi_{scr}$', out_file=out_file)

    # extract data, set labels
    out_file = out_file.replace('.png', '_n.png')
    sims_to_plot = find_simulations(results_all, R_eq=1)
    data_all = {}
    for sim in sims_to_plot:
        data = sim['potential'][0], sim['potential'][3]
        n = sim['params']['n']
        label = f"$n = {n}$"
        data_all[label] = data


    plot_generic(data_all, yscale='log', xmin=1E-1, xmax=1E2, xlabel=u'$r/R_s$',
                 ylabel=u'$\Phi_{scr,eff}/\Phi_{scr}$', out_file=out_file)

def main():
    run_plot_star_pot(out_file='starlike.png')
    run_plot_star_for(out_file='starlike_forces.png')
    run_plot_nfw_for(out_file='nfwlike_forces.png')
    run_plot_nfw_pot_eff(out_file='nfwlike_pot_eff.png')


if __name__ == "__main__":
    main()