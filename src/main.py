#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Tuesday June 16 2020
This code was implemented by
Louis Weyland, Hildebert Mouilé, Philippe Nicolau & Binjie Zhou.
"""
import sys
import os
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
import matplotlib.patches as mpatches

# Automatically setting the local path to this repo for easy file writing and saving
dir_path = path.dirname(path.realpath(__file__))
# Only keep the warnings printed in the output
warnings.filterwarnings("ignore")


def comp_average_fitness(size=(20, 20), iteration=2000, repetition=10, std=0.3):
    """
    Plots the average fitness for different distribution and the threshold
    :param : number of iterations, number of repetition and standard deviation for gaussian distribution
    """

    # Get a comparison between the different random distribution
    iterations = iteration
    repetition = repetition

    uniform_list = []
    gaussian_list = []
    threshold_uniform = []
    threshold_gaussian = []

    for i in range(repetition):

        uniform = Lattice(size=size, torus_mode=True, rand_dist=('uniform',), free_percent=0, iterations=iterations,
                          age_fraction=1 / 10)
        gaussian = Lattice(size=size, torus_mode=True, rand_dist=('gauss', 0.5, std), free_percent=0,
                           iterations=iterations, age_fraction=1 / 10)

        uniform.run(["mutation"])
        gaussian.run(["mutation"])

        uniform_list.append(uniform.average_fit_list)
        gaussian_list.append(gaussian.average_fit_list)
        threshold_uniform = threshold_uniform + uniform.threshold_list['threshold']
        threshold_gaussian = threshold_gaussian + gaussian.threshold_list['threshold']

    # get the average
    average_uniform = np.average(uniform_list, axis=0)
    average_gaussian = np.average(gaussian_list, axis=0)

    threshold_bar_uniform = np.ones(len(average_uniform)) * max(threshold_uniform)
    threshold_bar_gaussian = np.ones(len(average_gaussian)) * max(threshold_gaussian)

    # plot for comparision
    plot_setting()

    plt.plot(average_uniform, label='Uniform Distribution')
    plt.plot(average_gaussian, label='Gaussian Distribution')
    plt.plot(np.linspace(0, len(average_uniform), len(average_uniform)),
             threshold_bar_uniform, label='Threshold for Uniform Distribution')
    plt.plot(np.linspace(0, len(average_gaussian), len(average_gaussian)), threshold_bar_gaussian,
             label='Threshold for Gaussian Distribution with std ={}'.format(std))

    plt.legend()
    plt.title("Average Fitness over the different time step")
    plt.xlabel("Time steps (a.u.)")
    plt.ylabel("Fitness (a.u.)")
    plt.grid()
    plt.tight_layout()
    plt.savefig(path.join(dir_path,
                          'figures/average_fitness_s={}_itr={}_rep={}_std={}.png'.format(size,
                                                                                         iteration,
                                                                                         repetition,
                                                                                         std)),
                dpi=300)
    plt.show()


def comp_avalanche_time(size=(20, 20), iteration=2000, repetition=10, std=0.2):
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
                           iterations=iterations, age_fraction=1 / 10)

        uniform.run(["mutation", "avalanche_time"])
        gaussian.run(["mutation", "avalanche_time"])

        avalanche_uniform_list = avalanche_uniform_list + uniform.avalanche_time_list['avalanche_time']
        avalanche_gaussian_list = avalanche_gaussian_list + gaussian.avalanche_time_list['avalanche_time']

    result_uniform = powerlaw.Fit(avalanche_uniform_list, discrete=True, verbose=False)
    R_unifrom, p_uniform = result_uniform.distribution_compare('power_law', 'lognormal', normalized_ratio=True)
    result_gaussian = powerlaw.Fit(avalanche_gaussian_list, discrete=True, verbose=False)
    R_gaussian, p_gaussian = result_gaussian.distribution_compare('power_law', 'lognormal', normalized_ratio=True)

    # plot for comparision
    plot_setting()
    powerlaw.plot_pdf(avalanche_gaussian_list, color='b', label='Gaussian Random Distribution')
    powerlaw.plot_pdf(avalanche_uniform_list, color='r', label='Uniform Random Distribution')

    #plt.plot(average_gaussian,label='Gaussian Distribution')
    plt.legend()
    plt.title("Avalanche sizes")
    plt.xlabel("Probability (a.u.)")
    plt.ylabel("Avalanche sizes (a.u.)")
    plt.grid()
    plt.tight_layout()
    plt.savefig(path.join(dir_path,
                          'figures/avalanche_time_s={}_itr={}_rep={}_std={}.png'.format(size,
                                                                                        iteration,
                                                                                        repetition,
                                                                                        std)),
                dpi=300)
    plt.show()

    print("The slope with a uniform distribution is {}".format(result_uniform.power_law.alpha))
    print(
        "If {} > 0, the distribution of the data with uniform distribution resembles more a powerlaw than exponential distribution \n"
        "with a p value of {}".format(
            R_unifrom,
            p_uniform))
    print("The slope with a Gaussian distribtion is {}".format(result_gaussian.power_law.alpha))
    print(
        "If {} > 0, the distribution of the data with a gaussian distribution resembles more a powerlaw than exponential distribution \n"
        "with a p value of {}".format(
            R_gaussian,
            p_gaussian))


def comp_mutation_dist(size=(20, 20), iteration=2000, repetition=10, std=0.2):
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

        mutation_dist_uniform_list = mutation_dist_uniform_list + uniform.distance_btw_mutation_list
        mutation_dist_gaussian_list = mutation_dist_gaussian_list + gaussian.distance_btw_mutation_list

    result_uniform = powerlaw.Fit(mutation_dist_uniform_list, discrete=True, verbose=False)
    R_unifrom, p_uniform = result_uniform.distribution_compare('power_law', 'lognormal', normalized_ratio=True)
    result_gaussian = powerlaw.Fit(mutation_dist_gaussian_list, discrete=True, verbose=False)
    R_gaussian, p_gaussian = result_gaussian.distribution_compare('power_law', 'lognormal', normalized_ratio=True)

    n_uniform, bins_uniform = np.histogram(mutation_dist_uniform_list, density=True)
    n_gaussian, bins_gaussian = np.histogram(mutation_dist_gaussian_list, density=True)

    # plot for comparision
    plot_setting()
    plt.plot(bins_uniform[:-1], n_uniform, label='Uniform Distribution')
    plt.plot(bins_gaussian[:-1], n_gaussian, label='Gaussian Distribution')

    plt.legend()
    plt.title("Distribution of the distances between consecutive mutations")
    plt.xlabel("Probability (a.u.)")
    plt.ylabel("Distances between consecutive mutations (a.u.)")
    plt.yscale('log')
    plt.xscale('log')
    plt.grid()
    plt.tight_layout()
    plt.savefig(path.join(dir_path,
                          'figures/mutation_distance_s={}_itr={}_rep={}_std={}.png'.format(size,
                                                                                           iteration,
                                                                                           repetition,
                                                                                           std)),
                dpi=300)
    plt.show()

    print("The slope with a uniform distribution is {}".format(result_uniform.power_law.alpha))
    print(
        "If {} > 0, the distribution of the data with uniform distribution resembles more a powerlaw than exponential distribution \n"
        "with a p value of {}".format(
            R_unifrom,
            p_uniform))
    print("The slope with a Gaussian distribtion is {}".format(result_gaussian.power_law.alpha))
    print(
        "If {} > 0, the distribution of the data with a gaussian distribution resembles more a powerlaw than exponential distribution \n"
        "with a p value of {}".format(
            R_gaussian,
            p_gaussian))


def comp_diff_neighbours(size=(20, 20), iteration=2000, repetition=10):
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
        moore = Lattice(
            size=size,
            torus_mode=True,
            neighbourhood='Moore',
            rand_dist=(
                'uniform',
            ),
            free_percent=0,
            iterations=iterations,
        )
        vonNeumann = Lattice(
            size=size,
            torus_mode=True,
            neighbourhood='vonNeumann',
            rand_dist=(
                'uniform',
            ),
            free_percent=0,
            iterations=iterations,
        )

        moore.run(["mutation", "avalanche_time", "get_dist_btw_mutation"])
        vonNeumann.run(["mutation", "avalanche_time", "get_dist_btw_mutation"])

        avalanche_moore_list = avalanche_moore_list + moore.avalanche_time_list['avalanche_time']
        avalanche_vonNeumann_list = avalanche_vonNeumann_list + vonNeumann.avalanche_time_list['avalanche_time']

        mutation_dist_moore_list = mutation_dist_moore_list + moore.distance_btw_mutation_list
        mutation_dist_vonNeumann_list = mutation_dist_vonNeumann_list + vonNeumann.distance_btw_mutation_list

    result_moore = powerlaw.Fit(avalanche_moore_list, discrete=True, verbose=False)
    R_moore, p_moore = result_moore.distribution_compare('power_law', 'exponential', normalized_ratio=True)
    result_vonNeumann = powerlaw.Fit(avalanche_vonNeumann_list, discrete=True, verbose=False)
    R_vonNeumann, p_vonNeumann = result_vonNeumann.distribution_compare(
        'power_law', 'exponential', normalized_ratio=True)

    # plot for comparision
    plot_setting()
    powerlaw.plot_pdf(avalanche_moore_list, color='b', label='Moore')
    powerlaw.plot_pdf(avalanche_vonNeumann_list, color='r', label='vonNeuman')

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
    result_moore = powerlaw.Fit(mutation_dist_moore_list, discrete=True, verbose=False)
    R_moore, p_moore = result_moore.distribution_compare('power_law', 'exponential', normalized_ratio=True)
    result_vonNeumann = powerlaw.Fit(mutation_dist_vonNeumann_list, discrete=True, verbose=False)
    R_vonNeumann, p_vonNeumann = result_vonNeumann.distribution_compare(
        'power_law', 'exponential', normalized_ratio=True)
    print("The slope with a unifrom distribtion is {}".format(result_moore.power_law.alpha))
    print(
        "If {} > 0, the distribution of the data with Moore neighbourhood resembles more a powerlaw than exponential distribution \n"
        "with a p value of {}".format(
            R_moore,
            p_moore))
    print("The slope with a Gaussian distribtion is {}".format(result_vonNeumann.power_law.alpha))
    print(
        "If {} > 0, the distribution of the data with Moore neighbourhood resembles more a powerlaw than exponential distribution \n"
        "with a p value of {}".format(
            R_vonNeumann,
            p_vonNeumann))

    n_moore, bins_moore = np.histogram(mutation_dist_moore_list, density=True)
    n_vonNeumann, bins_vonNeumann = np.histogram(mutation_dist_vonNeumann_list, density=True)

    # plot for comparision
    plot_setting()
    plt.plot(bins_moore[:-1], n_moore, label='Moore Neighbourhood')
    plt.plot(bins_vonNeumann[:-1], n_vonNeumann, label='vonNeumann Neighbourhood')

    plt.legend()
    plt.title("Distribution of the distances between consecutive mutations")
    plt.xlabel("Probability (a.u.)")
    plt.ylabel("Distances between consecutive mutations (a.u.)")
    plt.yscale('log')
    plt.xscale('log')
    plt.grid()
    plt.tight_layout()
    plt.savefig(path.join(dir_path, 'figures/diff_neighbours_s={}_itr={}_rep={}.png'.format(size, iteration, repetition)), dpi=300)
    plt.show()

    print("The slope with a unifrom distribtion is {}".format(result_moore.power_law.alpha))
    print(
        "If {} > 0, the distribution of the data with Moore neighbourhood resembles more a powerlaw than exponential distribution \n"
        "with a p value of {}".format(
            R_moore,
            p_moore))
    print("The slope with a Gaussian distribtion is {}".format(result_vonNeumann.power_law.alpha))
    print(
        "If {} > 0, the distribution of the data with Moore neighbourhood resembles more a powerlaw than exponential distribution \n"
        "with a p value of {}".format(
            R_vonNeumann,
            p_vonNeumann))


def is_free_variation(i_min=0, i_max=1, i_iter=6, iterations = 2000):
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

    # figure & settings for plots
    plot_setting()
    plt.figure(1)
    plt.title('Avalanche times for different percentages of empty space in the lattice')
    plt.xlabel('Avalanche Time')
    plt.ylabel('Instances')
    plt.figure(2)
    plt.title('Threshold values for different percentages of empty space in the lattice')

    # looping over the different percentages
    for i in np.linspace(i_min, i_max, i_iter):

        i = round(i, 1)
        free_iter = Lattice(size=(20, 20), torus_mode=True, rand_dist=('uniform',),
                            free_percent=i, iterations=iterations, age_fraction=1 / 10)
        free_iter.run("all")

        av_times = free_iter.avalanche_time_list['avalanche_time']
        thresholds = free_iter.threshold_list['threshold']
        thresh_time = free_iter.threshold_list['time_step']

        avalanche_bins = range(min(av_times), max(av_times) + 1)
        threshold_bins = np.linspace(min(thresholds), max(thresholds), len(thresholds))

        # plt.xscale('log')
        # plt.yscale('log')
        #sb.distplot(av_times, label= str(i), hist=True)
        plt.figure(1)
        powerlaw.plot_pdf(av_times)
        #plt.hist(av_times, avalanche_bins, label= str(i))
        blue = mpatches.Patch(color='b', label='0.0')
        orange = mpatches.Patch(color='orange', label='0.2')
        green = mpatches.Patch(color='green', label='0.4')
        red = mpatches.Patch(color='red', label='0.6')
        purple = mpatches.Patch(color='purple', label='0.8')
        plt.legend(handles=[blue, orange, green, red, purple])

        plt.figure(2)
        plt.plot(thresh_time, thresholds, label=str(i))
        plt.xlabel('Iteration Number')
        plt.ylabel('Threshold Fitness Level')
        plt.legend(loc='upper right')

    plt.show()
    plt.grid()
    plt.tight_layout()
    plt.savefig(path.join(dir_path, 'figures/free_variation_imin={}_imax={}_iterations={}.png'.format(i_min, i_max, i_iter)), dpi=300)
    plt.show()
