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
import powerlaw
import warnings
warnings.filterwarnings("ignore")


def comp_average_fitness(size=(20, 20), iteration=2000, repetition=10, std=0.3):
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
        threshold_uniform = threshold_uniform + uniform.threshold_list['threshold']
        threshold_gaussian = threshold_gaussian + gaussian.threshold_list['threshold']

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
    plt.grid()
    plt.tight_layout()

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

    result_uniform = powerlaw.Fit(avalanche_uniform_list, discrete = True)
    R_unifrom, p_uniform = result_uniform.distribution_compare('power_law', 'lognormal', normalized_ratio = True)
    result_gaussian = powerlaw.Fit(avalanche_gaussian_list, discrete = True)
    R_gaussian, p_gaussian = result_gaussian.distribution_compare('power_law', 'lognormal', normalized_ratio = True)
    print("The slope with a uniform distribution is {}".format(result_uniform.power_law.alpha))
    print("If {} > 0, the distribution of the data with uniform distribution resembles more a powerlaw than exponential distribution \n"
          "with a p value of {}".format(R_unifrom,p_uniform))
    print("The slope with a Gaussian distribtion is {}".format(result_gaussian.power_law.alpha))
    print("If {} > 0, the distribution of the data with a gaussian distribution resembles more a powerlaw than exponential distribution \n"
          "with a p value of {}".format(R_gaussian,p_gaussian))


    # plot for comparision
    plot_setting()
    powerlaw.plot_pdf(avalanche_gaussian_list, color='b',label='Gaussian Random Distribution')
    powerlaw.plot_pdf(avalanche_uniform_list,color='r',label='Uniform Random Distribution')


    #plt.plot(average_gaussian,label='Gaussian Distribution')
    plt.legend()
    plt.title("Avalanche sizes")
    plt.xlabel("Probability (a.u.)")
    plt.ylabel("Avalanche sizes (a.u.)")
    plt.grid()
    plt.tight_layout()

    plt.show()


def comp_mutation_dist(size=(20, 20),iteration = 2000,repetition = 10 ,std = 0.2):
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

    result_uniform = powerlaw.Fit(mutation_dist_uniform_list, discrete = True)
    R_unifrom, p_uniform = result_uniform.distribution_compare('power_law', 'lognormal', normalized_ratio = True)
    result_gaussian = powerlaw.Fit(mutation_dist_gaussian_list, discrete = True)
    R_gaussian, p_gaussian = result_gaussian.distribution_compare('power_law', 'lognormal', normalized_ratio = True)
    print("The slope with a uniform distribution is {}".format(result_uniform.power_law.alpha))
    print("If {} > 0, the distribution of the data with uniform distribution resembles more a powerlaw than exponential distribution \n"
          "with a p value of {}".format(R_unifrom,p_uniform))
    print("The slope with a Gaussian distribtion is {}".format(result_gaussian.power_law.alpha))
    print("If {} > 0, the distribution of the data with a gaussian distribution resembles more a powerlaw than exponential distribution \n"
          "with a p value of {}".format(R_gaussian,p_gaussian))

    n_uniform, bins_uniform = np.histogram( mutation_dist_uniform_list, density=True)
    n_gaussian, bins_gaussian = np.histogram( mutation_dist_gaussian_list, density=True)

    # plot for comparision
    plot_setting()
    plt.plot(bins_uniform[:-1], n_uniform,label='Uniform Distribution')
    plt.plot(bins_gaussian[:-1], n_gaussian,label='Gaussian Distribution')

    # plt.plot(average_gaussian,label='Gaussian Distribution')
    plt.legend()
    plt.title("Distribution of the distances between consecutive mutations")
    plt.xlabel("Probability (a.u.)")
    plt.ylabel("Distances between consecutive mutations (a.u.)")
    plt.yscale('log')
    plt.xscale('log')
    plt.grid()
    plt.tight_layout()

    plt.show()



def comp_diff_neighbours(size=(20, 20),iteration = 2000,repetition = 10):
    """
    Plots the avalanche distribution in a log-log plot
    :param : number of iterations, number of repetition and standard deviation for gaussian distribution

    """
    # Get a comparison between the different random distribution
    iterations = iteration
    repetition = repetition

    mutation_dist_vonNeumann_list = []
    mutation_dist_moore_list = []
    avalanche_moore_list = []
    avalanche_vonNeumann_list = []

    for i in range(repetition):
        moore = Lattice(size=size, torus_mode=True,neighbourhood='Moore', rand_dist=('uniform',), free_percent=0, iterations=iterations,
                        )
        vonNeumann = Lattice(size=size, torus_mode=True,neighbourhood='vonNeumann', rand_dist=('uniform',), free_percent=0,
                          iterations=iterations,)

        moore.run(["mutation","avalanche_time","get_dist_btw_mutation"])
        vonNeumann.run(["mutation","avalanche_time","get_dist_btw_mutation"])

        avalanche_moore_list = avalanche_moore_list + moore.avalanche_time_list['avalanche_time']
        avalanche_vonNeumann_list = avalanche_vonNeumann_list + vonNeumann.avalanche_time_list['avalanche_time']

        mutation_dist_moore_list =  mutation_dist_moore_list + moore.distance_btw_mutation_list
        mutation_dist_vonNeumann_list = mutation_dist_vonNeumann_list + vonNeumann.distance_btw_mutation_list

    result_moore = powerlaw.Fit(avalanche_moore_list,discrete = True)
    R_moore, p_moore = result_moore.distribution_compare('power_law', 'exponential', normalized_ratio = True)
    result_vonNeumann = powerlaw.Fit(avalanche_vonNeumann_list, discrete = True)
    R_vonNeumann, p_vonNeumann = result_vonNeumann.distribution_compare('power_law', 'exponential', normalized_ratio=True)
    print("The slope with a unifrom distribtion is {}".format(result_moore.power_law.alpha))
    print("If {} > 0, the distribution of the data with Moore neighbourhood resembles more a powerlaw than exponential distribution \n"
          "with a p value of {}".format(R_moore,p_moore))
    print("The slope with a Gaussian distribtion is {}".format(result_vonNeumann.power_law.alpha))
    print("If {} > 0, the distribution of the data with Moore neighbourhood resembles more a powerlaw than exponential distribution \n"
          "with a p value of {}".format(R_vonNeumann,p_vonNeumann))

    # plot for comparision
    plot_setting()
    powerlaw.plot_pdf(avalanche_moore_list, color='b',label='Moore')
    powerlaw.plot_pdf(avalanche_vonNeumann_list,color='r',label='vonNeuman')

    plt.legend()
    plt.title("Avalanche sizes")
    plt.xlabel("Probability (a.u.)")
    plt.ylabel("Avalanche sizes (a.u.)")
    plt.yscale('log')
    plt.xscale('log')
    plt.grid()
    plt.tight_layout()
    
    plt.show()

    # new figure
    result_moore = powerlaw.Fit(mutation_dist_moore_list,discrete = True)
    R_moore, p_moore = result_moore.distribution_compare('power_law', 'exponential', normalized_ratio = True)
    result_vonNeumann = powerlaw.Fit(mutation_dist_vonNeumann_list, discrete = True)
    R_vonNeumann, p_vonNeumann = result_vonNeumann.distribution_compare('power_law', 'exponential', normalized_ratio=True)
    print("The slope with a unifrom distribtion is {}".format(result_moore.power_law.alpha))
    print("If {} > 0, the distribution of the data with Moore neighbourhood resembles more a powerlaw than exponential distribution \n"
          "with a p value of {}".format(R_moore,p_moore))
    print("The slope with a Gaussian distribtion is {}".format(result_vonNeumann.power_law.alpha))
    print("If {} > 0, the distribution of the data with Moore neighbourhood resembles more a powerlaw than exponential distribution \n"
          "with a p value of {}".format(R_vonNeumann,p_vonNeumann))

    #n_moore, bins_moore = np.histogram(avalanche_moore_list, density=True)
    #n_gaussian, bins_gaussian = np.histogram(avalanche_gaussian_list, density=True)

    # plot for comparision
    plot_setting()
    #plt.plot(bins_uniform[:-1], n_uniform,label='Uniform Distribution')
    #plt.plot(bins_gaussian[:-1], n_gaussian,label='Gaussian Distribution')
    powerlaw.plot_pdf(mutation_dist_moore_list, color='b',label='Moore')
    powerlaw.plot_pdf(mutation_dist_vonNeumann_list,color='r',label='vonNeuman')

    #plt.plot(average_gaussian,label='Gaussian Distribution')
    plt.legend()
    plt.title("Avalanche sizes")
    plt.xlabel("Probability (a.u.)")
    plt.ylabel("Avalanche sizes (a.u.)")
    plt.yscale('log')
    plt.xscale('log')
    plt.show()

comp_average_fitness(iteration=2000)