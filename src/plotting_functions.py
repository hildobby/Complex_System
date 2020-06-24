#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Tuesday June 16 2020
This code was implemented by
Louis Weyland, Hildebert Mouilé, Philippe Nicolau & Binjie Zhou.
"""
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap


def plot_setting():
    """
    Change the default matplotlib parameters to make sure the graphs are unifrom
    """
    params = {'legend.fontsize': 12,
              'figure.figsize': (10, 6),
              'axes.labelsize': 14,
              'axes.titlesize': 20,
              'xtick.labelsize': 12,
              'ytick.labelsize': 12}
    plt.rcParams.update(params)

def plt_color(i):
    """
    Return MatPlotLib default color from cycle at the provided index
    """

    colors = np.array(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    return colors[np.array(i) % len(colors)]