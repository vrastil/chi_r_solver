import numpy as np

from .sim_chameleon import run_many_sims
from .plot_chameleon import plot_generic

kwargs_dflt = {
    "err": 1E-12,
    "n": 0.5,
    "print_par" : 0,
    "out_dir" : "../output/",
    "step": 0.01,
}

def run_plot_star_pot():
    print(f"Running STAR-like simulation.")

    kwargs_sims = {
        "mod" : 0,
        "Omega_m" : 1e6,
        "M200_sun" : 3E-16,
        "Ys" : [1E-13, 1E-14, 1E-15],
        "R": 1E-1,
    }

    results_all = run_many_sims(kwargs_dflt, kwargs_sims, stdout=None)
    print(f"Plot and save.")

    # extract data, set labels
    data_all = {}
    for sim in results_all:
        data = sim['potential'][0], sim['potential'][2]
        Ys = int(np.log10(sim['params']['Ys']))
        label = f"$\Phi_s = 10^{{{Ys}}}$"
        data_all[label] = data


    plot_generic(data_all, ymin=1E-16, ymax=1E-13, xlabel=u'$r/R_s$', ylabel=u'$\chi$',
                 out_file='starlike.png')

def main():
    run_plot_star_pot()

if __name__ == "__main__":
    main()