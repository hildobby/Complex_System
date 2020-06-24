#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Tuesday June 16 2020
This code was implemented by
Louis Weyland, Hildebert Mouilé, Philippe Nicolau & Binjie Zhou.
"""

#from os import path
#import matplotlib.pyplot as plt
#import matplotlib.animation
#from plotting_functions import plot_setting
from lattice import Lattice
from plotting_functions import plot_setting
import matplotlib.pyplot as plt
import time
import numpy as np
#import seaborn as sb
import math

# Get a comparison between the different random distribution

iterations = 2000
t0 = time.time()
# if rand_dist take 1 arg, rand_dist=('uniform',) !! Comma needed here



# =============================================================================
# repetition = 10
# 
# uniform_list = []
# gaussian_list = []
# exponential_list = []
# 
# for i in range(repetition):
#     print("Repetiton {}".format(i))
# 
#     uniform = Lattice(size=(20, 20), torus_mode=True, rand_dist=('uniform',), free_percent=0, iterations=iterations,
#                       age_fraction=1 / 10)
#     gaussian = Lattice(size=(20, 20), torus_mode=True, rand_dist=('gaussian', 0.5, 0.2), free_percent=0,
#                       iterations=iterations, age_fraction=1 / 10)
#     exponential = Lattice(size=(20, 20), torus_mode=True, rand_dist=('exponential', 1), free_percent=0,
#                           iterations=iterations, age_fraction=1 / 10)
# 
# 
#     uniform.run(["mutation"])
#     gaussian.run("mutation")
#     exponential.run("mutation")
# 
#     uniform_list.append(uniform.average_fit_list)
#     gaussian_list.append(gaussian.average_fit_list)
#     exponential_list.append(exponential.average_fit_list)
# 
# 
# 
# t1 = time.time()
# print("TOTAL TIME NEEDED {}".format(t1 - t0))
# 
# 
# =============================================================================



def is_free_variation(i_min,i_max,i_iter):
    '''
    runs several instances of the lattice 
    with different percentages of empty nodes in the lattice.
    avalanche time and thresholds are then compared between runs.
    ____
    
    i_min: lower percentage of the range (within [0,1])
    i_max: upper percentage (within [0,1])
    i_iter: number of steps between i_min and i_max to be taken (Integer)
    '''
# =============================================================================
#     #list of thresholds
#     free_thresh = []
#     #list of avalanche times
#     free_avalanche = []
# =============================================================================
    
    
    #figure & settings for plots
    plot_setting()
    plt.figure(1) 
    plt.title('Avalanche times for different percentages of empty space in the lattice')
    plt.xlabel('Avalanche Time')
    plt.ylabel('Instances')
    plt.figure(2)
    plt.title('Threshold values for different percentages of empty space in the lattice')
    

    #looping over the different percentages
    for i in np.linspace(i_min,i_max,i_iter):
        
        i = round(i,1)
        free_iter = Lattice(size=(20, 20), torus_mode=True, rand_dist=('uniform',), free_percent=i, iterations=iterations,
                      age_fraction=1 / 10)
        free_iter.run("all")
        
        av_times = free_iter.avalanche_time_list['avalanche_time']
        thresholds = free_iter.threshold_list
        
        avalanche_bins = range(min(av_times),max(av_times)+1)
        threshold_bins = np.linspace(min(thresholds),max(thresholds),len(thresholds))
        
        #plt.xscale('log')
        #plt.yscale('log')
        #sb.distplot(av_times, label= str(i), hist=True)
        plt.figure(1)
        plt.hist(av_times, avalanche_bins, label= str(i))
        plt.legend(loc='upper right')
        
        plt.figure(2)
        plt.plot(range(len(thresholds)), thresholds, label= str(i))
        plt.legend(loc='upper right')
        
    plt.show()
        
        
    
        
        
# =============================================================================
#         free_thresh.append([i,free_iter.threshold_list])
#         free_avalanche.append([i,free_iter.avalanche_time_list])
# =============================================================================
        
    

    
    
    
    