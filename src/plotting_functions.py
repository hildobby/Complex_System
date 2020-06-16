#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Tuesday June 16 2020
This code was implemented by
Louis Weyland, Hildebert Mouilé, Philippe Nicolau & Binjie Zhou.
"""

import matplotlib.pyplot as plt
import numpy as np

def plot_setting():
    """
    Change the default matplotlib parameters to make sure the graphs are unifrom
    """
    params = {'legend.fontsize': 'x-large',
              'figure.figsize': (8, 6),
              'axes.labelsize': 'x-large',
              'axes.titlesize': 'x-large',
              'xtick.labelsize': 'x-large',
              'ytick.labelsize': 'x-large'}
    plt.rcParams.update(params)

def plt_color(i):
    """
    Return MatPlotLib default color from cycle at the provided index
    """
    
    colors = np.array(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    return colors[np.array(i) % len(colors)]