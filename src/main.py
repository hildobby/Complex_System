#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Tuesday June 16 2020
This code was implemented by
Louis Weyland, Hildebert Mouilé, Philippe Nicolau & Binjie Zhou.
"""

from os import path
import matplotlib.pyplot as plt
import matplotlib.animation
if path.isdir("src"):
    from src.plotting_functions import plot_setting
    from src.lattice import Lattice
else:
    from plotting_functions import plot_setting
    from lattice import Lattice

import time
import numpy as np


def comp_average_fitness(size=(20, 20),iteration = 2000,repetition = 10 , std = 0.3):
    """
    Plots the average fintess for different distribution and the threshold
    :param : number of iterations, number of repetition and standard deviation for gaussian distribution

    """
    # Get a comparison between the different random distribution
    iterations = iteration
    repetition = repetition

    uniform_list= []
    gaussian_list = []
    threshold_uniform = []
    threshold_gaussian = []

    for i in range(repetition):

        uniform = Lattice(size=size, torus_mode=True, rand_dist=('uniform',), free_percent=0, iterations=iterations,
                          age_fraction=1 / 10)
        gaussian = Lattice(size=size, torus_mode=True, rand_dist=('gauss', 0.5, std), free_percent=0,
                          iterations=iterations, age_fraction=1/10)

        uniform.run(["mutation"])
        gaussian.run(["mutation"])

        uniform_list.append(uniform.average_fit_list)
        gaussian_list.append(gaussian.average_fit_list)
        threshold_uniform = threshold_uniform + uniform.threshold_list
        threshold_gaussian = threshold_gaussian + gaussian.threshold_list

    # get the average
    average_uniform = np.average(uniform_list,axis=0)
    average_gaussian = np.average(gaussian_list,axis=0)

    threshold_bar_uniform = np.ones(len(average_uniform))*max(threshold_uniform)
    threshold_bar_gaussian = np.ones(len(average_gaussian))*max(threshold_gaussian)

    # plot for comparision
    plot_setting()

    plt.plot(average_uniform,label='Uniform Distribution')
    plt.plot(average_gaussian,label='Gaussian Distribution')
    plt.plot(np.linspace(0,len(average_uniform),len(average_uniform)),threshold_bar_uniform,label='Threshold for Uniform Distribution')
    plt.plot(np.linspace(0,len(average_gaussian),len(average_gaussian)),threshold_bar_gaussian,
             label='Threshold for Gaussian Distribution with std ={}'.format(std))

    plt.legend()
    plt.title("Average Fitness over the different time step")
    plt.xlabel("Time steps (a.u.)")
    plt.ylabel("Fitness (a.u.)")

    plt.show()

def comp_avalanche_time(size=(20, 20),iteration = 2000,repetition = 10 , std = 0.2):
    """
    Plots the avalanche distribution in a log-log plot
    :param : number of iterations, number of repetition and standard deviation for gaussian distribution

    """
    # Get a comparison between the different random distribution
    iterations = iteration
    repetition = repetition

    avalanche_uniform_list = []
    avalanche_gaussian_list = []

    for i in range(repetition):
        uniform = Lattice(size=size, torus_mode=True, rand_dist=('uniform',), free_percent=0, iterations=iterations,
                          age_fraction=1 / 10)
        gaussian = Lattice(size=size, torus_mode=True, rand_dist=('gauss', 0.5, std), free_percent=0,
                          iterations=iterations, age_fraction=1/10)

        uniform.run(["mutation","avalanche_time"])
        gaussian.run(["mutation","avalanche_time"])

        avalanche_uniform_list = avalanche_uniform_list+ uniform.avalanche_time_list['avalanche_time']
        avalanche_gaussian_list = avalanche_gaussian_list+ gaussian.avalanche_time_list['avalanche_time']

    n_uniform, bins_uniform = np.histogram(avalanche_uniform_list, density=True)
    n_gaussian, bins_gaussian = np.histogram(avalanche_gaussian_list, density=True)

    # plot for comparision
    plot_setting()
    plt.plot(bins_uniform[:-1], n_uniform,label='Uniform Distribution')
    plt.plot(bins_gaussian[:-1], n_gaussian,label='Gaussian Distribution')

    #plt.plot(average_gaussian,label='Gaussian Distribution')
    plt.legend()
    plt.title("Avalanche sizes")
    plt.xlabel("Probability (a.u.)")
    plt.ylabel("Avalanche sizes (a.u.)")
    plt.yscale('log')
    plt.xscale('log')
    plt.show()


def comp_mutation_dist(size=(20, 20),iteration = 2000,repetition = 10 , std = 0.2):
    """
    Plots the distribution between distances between mutations
    :param : size of the grid,number of iterations, number of repetition and standard deviation for gaussian distribution
    """


    iterations = iteration
    repetition = repetition

    mutation_dist_uniform_list = []
    mutation_dist_gaussian_list = []

    for i in range(repetition):
        uniform = Lattice(size=size, torus_mode=True, rand_dist=('uniform',), free_percent=0, iterations=iterations,
                          age_fraction=1 / 10)
        gaussian = Lattice(size=size, torus_mode=True, rand_dist=('gauss', 0.5, std), free_percent=0,
                           iterations=iterations, age_fraction=1 / 10)
    
        uniform.run(["mutation", "get_dist_btw_mutation"])
        gaussian.run(["mutation", "get_dist_btw_mutation"])
    
        mutation_dist_uniform_list =  mutation_dist_uniform_list + uniform.distance_btw_mutation_list
        mutation_dist_gaussian_list = mutation_dist_gaussian_list + gaussian.distance_btw_mutation_list

    n_uniform, bins_uniform = np.histogram(mutation_dist_uniform_list, density=True)
    n_gaussian, bins_gaussian = np.histogram(mutation_dist_gaussian_list, density=True)

    # plot for comparision
    plot_setting()
    plt.plot(bins_uniform[:-1], n_uniform, label='Uniform Distribution')
    plt.plot(bins_gaussian[:-1], n_gaussian, label='Gaussian Distribution')

    # plt.plot(average_gaussian,label='Gaussian Distribution')
    plt.legend()
    plt.title("Distribution of the distances between consecutive mutations")
    plt.xlabel("Probability (a.u.)")
    plt.ylabel("Distances between consecutive mutations (a.u.)")
    plt.yscale('log')
    plt.xscale('log')
    plt.show()


