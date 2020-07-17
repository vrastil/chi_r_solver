"""
'plot_chameleon.py' modules serves for plotting results of simulations
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import animation
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib.colors import SymLogNorm
from matplotlib.colors import ListedColormap
from matplotlib.ticker import FormatStrFormatter

# color-blind
CB_colors = ['#377eb8', '#ff7f00', '#4daf4a',
            '#f781bf', '#a65628', '#984ea3',
            '#999999', '#e41a1c', '#dede00']
# newcmp = ListedColormap(CB_colors, name='ColorBlind')


# default values
matplotlib.rcParams['legend.numpoints'] = 1
matplotlib.rcParams['lines.linewidth'] = 4.0
matplotlib.rcParams['lines.markersize'] = 6.0
# matplotlib.rcParams['axes.prop_cycle'] = cycler(color=CB_colors)
matplotlib.rcParams['axes.labelsize'] = 30
matplotlib.rcParams['xtick.labelsize'] = 20
matplotlib.rcParams['ytick.labelsize'] = 20
matplotlib.rcParams['legend.fontsize'] = 25
matplotlib.rcParams['font.size'] = 25

dflt_suptitle_size = 25
dflt_fig_size = (14, 9)
dflt_fig_size_map = (14, 14)

def plot_generic(data_all, **kwargs):
    # extract arguments, set defaults
    fig_size = kwargs.get('figsize', dflt_fig_size)
    xscale = kwargs.get('xscale', 'log')
    yscale = kwargs.get('yscale', 'log')
    show = kwargs.get('show', True)
    xmin = kwargs.get('xmin', 0.1)
    xmax = kwargs.get('xmax', 10)
    ymin = kwargs.get('ymin', None)
    ymax = kwargs.get('ymax', None)
    out_dir = kwargs.get('out_dir', '../../dizertace/img/spherical_cham/')
    out_file = kwargs.get('out_file', None)
    xlabel = kwargs.get('xlabel', None)
    ylabel = kwargs.get('ylabel', None)
    
    # create figure
    fig = plt.figure(figsize=fig_size)
    ax = plt.gca()

    # plot all data
    for label, data in data_all.items():
        x, y = data
        ax.plot(x, y, label=label)

    # set scales
    ax.set_yscale(yscale)
    ax.set_xscale(xscale)

    # set limits
    ax.set_xlim(xmin=xmin, xmax=xmax)
    ax.set_ylim(ymin=ymin, ymax=ymax)

    # labels
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

    # legend
    plt.legend()

    # save, show, close
    if show:
        plt.show()

    if out_file:
        plt.savefig(out_dir + out_file)

    plt.close('all')