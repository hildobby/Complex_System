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
import argparse
from termcolor import colored

# Automatically setting the local path to this repo for easy file writing and saving
dir_path = path.dirname(path.realpath(__file__))
# Only keep the warnings printed in the output
warnings.filterwarnings("ignore")




parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
description='Uses the lattice object from lattice.py to perform different experiments.\n\
Call functions by their name!\n\
func \'comp_average_fitness\' runs the model and plots the evolution of average fitness and the threshold at the end\n\
func \'comp_avalanche_time\'  compares the avalanche time between fintess generated unifromly or by the gaussin random generator\n\
func \'comp_mutation_dist\'  compares the distribution of the distance between mutations for models using gaussin/uniform random generator\n\
func \'is_free_variation\' compares the impact of the density on the avalanche time   \n\
func \'comp_cluster_sizes\' compares the cluster size distribution for different grid sizes \n\
func \'comp_moving_vs_stationary\' compares the cluster sizes  and avalanche time between a \n\
func \'comp_diff_dim\' the cluster size for 2D/3D model \n\
func \'get_fitness_dist\' the fitness distribution after n iterations ')



parser.add_argument("-func",type = str, help='Defines which function to execute')
parser.add_argument('-itr', type=int,default=200, help='Number of iteration for the model (default : 50)')
parser.add_argument('-rep', type=int,default=10, help='Number of repetition (default : 10)')
parser.add_argument('-std', type=float,default=0.3, help='Standard deviation for gaussian random generator (default : 0.3)')
parser=parser.parse_args()



def print_statement(alpha, r, p, name):

    print("The slope of {} disrtibution is {}".format(name, round(alpha, 4)))
    if r < 0:
        print("{} follows a lognormal distribution with a p = {} ".format(name, round(p, 4)))
    if r > 0:
        print("{} follows a powerlaw distribution with a p = {} ".format(name, round(p, 4)))


def comp_average_fitness(size=(20, 20), iteration=2000, repetition=10, std=0.3):
    """
    Plots the average fitness for different distribution and the threshold for a model
    using uniform/gaussian random generator
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
    plt.xlabel("Time steps ")
    plt.ylabel("Fitness ")
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
    Plots the avalanche distribution in a log-log plot for a model using uniform/gaussian random generator
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
    plt.title("Avalanche Time")
    plt.ylabel("Probability ")
    plt.xlabel("Avalanche Time")
    plt.grid()
    plt.tight_layout()
    plt.savefig(path.join(dir_path,
                          'figures/avalanche_time_s={}_itr={}_rep={}_std={}.png'.format(size,
                                                                                        iteration,
                                                                                        repetition,
                                                                                        std)),
                dpi=300)
    plt.show()

    print_statement(result_uniform.power_law.alpha, R_unifrom, p_uniform, "uniform")
    print_statement(result_gaussian.power_law.alpha, R_gaussian, p_gaussian, "gaussian")


def comp_mutation_dist(size=(20, 20), iteration=2000, repetition=10, std=0.2):
    """
    Plots the distribution between distances between mutations for a model using uniform/gaussian random generator
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
    plt.ylabel("Probability ")
    plt.xlabel("Distances between consecutive mutations")
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

    print_statement(result_uniform.power_law.alpha, R_unifrom, p_uniform, "uniform")
    print_statement(result_gaussian.power_law.alpha, R_gaussian, p_gaussian, "gaussian")


def comp_diff_neighbours(size=(20, 20), iteration=2000, repetition=10):
    """
    Plots the avalanche distribution in a log-log plot  for a model using Moore/van Neumann Neighbourhood
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
            neighbourhood='von Neumann',
            rand_dist=(
                'uniform',
            ),
            free_percent=0,
            iterations=iterations,
        )
        moore = Lattice(size=size, torus_mode=True, neighbourhood='Moore', rand_dist=('uniform',), free_percent=0, iterations=iterations,
                        )
        vonNeumann = Lattice(size=(50, 50), torus_mode=True, neighbourhood='von Neumann', rand_dist=('uniform',), free_percent=0,
                             iterations=iterations,)

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
    powerlaw.plot_pdf(avalanche_vonNeumann_list, color='r', label='von Neumann')

    plt.legend()
    plt.title("Avalanche Time")
    plt.ylabel("Probability ")
    plt.xlabel("Avalanche Time")
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

    print_statement(result_moore.power_law.alpha, R_moore, p_moore, "More Neighbour")
    print_statement(result_vonNeumann.power_law.alpha, R_vonNeumann, p_vonNeumann, "von Neumann Neighbourhood")

    n_moore, bins_moore = np.histogram(mutation_dist_moore_list, density=True)
    n_vonNeumann, bins_vonNeumann = np.histogram(mutation_dist_vonNeumann_list, density=True)

    # plot for comparision
    plot_setting()
    plt.plot(bins_moore[:-1], n_moore, label='Moore Neighbourhood')
    plt.plot(bins_vonNeumann[:-1], n_vonNeumann, label='von Neumann Neighbourhood')

    plt.legend()
    plt.title("Distribution of the distances between consecutive mutations")
    plt.ylabel("Probability")
    plt.xlabel("Distances between consecutive mutations")
    plt.yscale('log')
    plt.xscale('log')
    plt.grid()
    plt.tight_layout()
    plt.savefig(path.join(dir_path, 'figures/diff_neighbours_s={}_itr={}_rep={}.png'.format(size, iteration, repetition)), dpi=300)
    plt.show()

    print_statement(result_moore.power_law.alpha, R_moore, p_moore, "More Neighbour")
    print_statement(result_vonNeumann.power_law.alpha, R_vonNeumann, p_vonNeumann, 'von Neumann')


def is_free_variation(i_min=0, i_max=1, i_iter=6, iterations=2000):
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


def comp_cluster_sizes(iterations=2000):
    """
     Compares the cluster sizes of different sizes of grid  using uniform random generator for the fitness
    :param iterations: number of steps to run for the Bak-Sneppen model
    """


    print(colored("Warning this function might take long (2000 itr ~ 40 min)",'red'))

    small = Lattice(size=(20, 20), torus_mode=True, rand_dist=('uniform',), free_percent=0, iterations=iterations,
                    age_fraction=1 / 10)
    medium = Lattice(size=(50, 50), torus_mode=True, rand_dist=('uniform',), free_percent=0, iterations=iterations,
                     age_fraction=1 / 10)

    large = Lattice(size=(70, 70), torus_mode=True, rand_dist=('uniform',), free_percent=0, iterations=iterations,
                    age_fraction=1 / 10)

    small.run(["mutation", "update_age", "get_cluster"])
    medium.run(["mutation", "update_age", "get_cluster"])
    large.run(["mutation", "update_age", "get_cluster"])

    small_hist = np.concatenate([small.cluster_size[x] for x in small.cluster_size])
    medium_hist = np.concatenate([medium.cluster_size[x] for x in medium.cluster_size])
    large_hist = np.concatenate([large.cluster_size[x] for x in large.cluster_size])

    # get the power law
    small_results = powerlaw.Fit(small_hist, discrete=True, verbose=False)
    medium_results = powerlaw.Fit(medium_hist, dicsrete=True, verbose=False)
    large_resutls = powerlaw.Fit(large_hist, discrete=True, verbose=False)

    r_small, p_small = small_results.distribution_compare('power_law', 'exponential', normalized_ratio=True)
    r_medium, p_medium = medium_results.distribution_compare('power_law', 'exponential', normalized_ratio=True)
    r_large, p_large = large_resutls.distribution_compare('power_law', 'exponential', normalized_ratio=True)

    # plot the power law
    plt.figure()
    powerlaw.plot_pdf(small_hist, label='20X20 Grid')
    powerlaw.plot_pdf(medium_hist, label='50X50 Grid')
    powerlaw.plot_pdf(large_hist, label='70X70 Grid')

    plt.title("Compare cluster size for different grid sizes")
    plt.xlabel("Cluster size ")
    plt.ylabel("Probability")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(path.join(dir_path, 'figures/cluster-sizes_rep={}.png'.format(iterations)), dpi=300)
    plt.show()

    print_statement(small_results.power_law.alpha, r_small, p_small, "the 20X20 grid's")
    print_statement(medium_results.power_law.alpha, r_medium, p_medium, "the 50X50 grid's")
    print_statement(large_resutls.power_law.alpha, r_large, p_large, "the 70X70 grid's")


def comp_moving_vs_stationary(size=(20, 20), iteration=2000, repetition=10):
    """
    Compares the cluster sizes  and avalanche time between a stationary model and a model where the nodes
    can move to free space
    :param : number of iterations, number of repetition and standard deviation for gaussian distribution
    """
    # Get a comparison between the different random distribution
    iterations = iteration
    repetition = repetition

    avalanche_move_list = []
    avalanche_stationary_list = []

    for i in range(repetition):
        stationary = Lattice(
            size=size,
            torus_mode=True,
            neighbourhood='Moore',
            rand_dist=(
                'uniform',
            ),
            free_percent=0,
            iterations=iterations,
        )
        move = Lattice(size=(50, 50), torus_mode=True, neighbourhood='Moore', rand_dist=('uniform',), free_percent=0.3,
                       iterations=iterations,)

        stationary.run(["mutation", "avalanche_time"])
        move.run(["moving", "avalanche_time"])

        avalanche_move_list = avalanche_move_list + move.avalanche_time_list['avalanche_time']
        avalanche_stationary_list = avalanche_stationary_list + stationary.avalanche_time_list['avalanche_time']

    result_move = powerlaw.Fit(avalanche_move_list, discrete=True, verbose=False)
    R_move, p_move = result_move.distribution_compare('power_law', 'exponential', normalized_ratio=True)
    result_stationary = powerlaw.Fit(avalanche_stationary_list, discrete=True, verbose=False)
    R_stationary, p_stationary = result_stationary.distribution_compare(
        'power_law', 'exponential', normalized_ratio=True)

    # plot for comparision
    plot_setting()
    powerlaw.plot_pdf(avalanche_move_list, color='b', label='Migration')
    powerlaw.plot_pdf(avalanche_stationary_list, color='r', label='No Migration')

    plt.legend()
    plt.title("Avalanche sizes")
    plt.ylabel("Probability ")
    plt.xlabel("Avalanche sizes ")
    plt.grid()
    plt.tight_layout()
    plt.savefig(path.join(dir_path, 'figures/moving-vs-stationary_size={}_itr{}_rep={}.png'.format(size, iteration, repetition)), dpi=300)
    plt.show()

    print_statement(result_move.power_law.alpha, R_move, p_move, "migration")
    print_statement(result_stationary.power_law.alpha, R_stationary, p_stationary, "no migration")


def comp_diff_dim(iterations=2000):
    """
    Compares the cluster size distribution for 2 Dimensions and 3 Dimensions
    """
    # Compares the cluster sizes of different sizes of grid

    grid = Lattice(size=(20, 20), torus_mode=True, rand_dist=('uniform',), free_percent=0, iterations=iterations,
                   age_fraction=1 / 10)
    cube = Lattice(size=(20, 20, 3), torus_mode=True, rand_dist=('uniform',), free_percent=0, iterations=iterations,
                   age_fraction=1 / 10)

    grid.run(["mutation", "update_age", "get_cluster"])
    cube.run(["mutation", "update_age", "get_cluster"])

    grid_hist = np.concatenate([grid.cluster_size[x] for x in grid.cluster_size])
    cube_hist = np.concatenate([cube.cluster_size[x] for x in cube.cluster_size])

    # get the power law
    grid_results = powerlaw.Fit(grid_hist, discrete=True, verbose=False)
    cube_results = powerlaw.Fit(grid_hist, dicsrete=True, verbose=False)

    r_grid, p_grid = grid_results.distribution_compare('power_law', 'exponential', normalized_ratio=True)
    r_cube, p_cube = cube_results.distribution_compare('power_law', 'exponential', normalized_ratio=True)

    # plot the power law
    plot_setting()
    powerlaw.plot_pdf(grid_hist, label='2 Dimensions')
    powerlaw.plot_pdf(cube_hist, label='3 Dimensions')
    plt.title("Cluster Distribution for 2D and 3D")
    plt.xlabel("Cluster size ")
    plt.ylabel("Probability ")
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(path.join(dir_path, 'figures/different_dimentions_itr={}.png'.format(iterations)), dpi=300)
    plt.show()

    print_statement(grid_results.power_law.alpha, r_grid, p_grid, "2D")
    print_statement(cube_results.power_law.alpha, r_cube, p_cube, "3D")


def get_fitness_dist(iterations=20000):
    """
    Get the distribution of fitness and situate it with respect to the threshold
    """

    lattice = Lattice(size=(40, 40), torus_mode=True, rand_dist=('uniform',),
                      free_percent=0, iterations=iterations, age_fraction=1 / 10)
    lattice.run(["mutation"])

    plot_setting()
    plt.figure()
    plt.hist(list(lattice.fitness_dict.values()), label='Fitness Distribution')
    plt.axvline(x=max(lattice.threshold_list['threshold']), color='red', label='Threshold')
    plt.xlim((0, 1))
    plt.xlabel('Fitness')
    plt.ylabel('Probability')
    plt.tight_layout()
    plt.legend()
    plt.savefig(path.join(dir_path, 'figures/fitness_distance_itr={}.png'.format(iterations)), dpi=300)
    plt.show()




if len(sys.argv) >= 1:

    if parser.func == 'comp_average_fitness':
        comp_average_fitness(iteration = parser.itr,repetition =parser.rep,std= parser.std)

    elif parser.func ==  'comp_average_fitness':
        comp_average_fitness(iteration = parser.itr,repetition =parser.rep,std= parser.std)

    elif parser.func == 'comp_avalanche_time':
        comp_avalanche_time(iteration = parser.itr,repetition =parser.rep,std= parser.std)

    elif parser.func == 'comp_mutation_dist' :
        comp_mutation_dist(iteration = parser.itr,repetition =parser.rep,std= parser.std)
    elif parser.func == 'is_free_variation' :
        is_free_variation(iterations = parser.itr)

    elif parser.func == 'comp_cluster_sizes':
        comp_cluster_sizes(iterations = parser.itr)

    elif parser.func == 'comp_moving_vs_stationary':
        comp_moving_vs_stationary(iteration = parser.itr,repetition =parser.rep)

    elif parser.func == 'comp_diff_dim':
        comp_diff_dim(iterations = parser.itr)

    elif parser.func == 'get_fitness_dist':
        get_fitness_dist(iterations = parser.itr)