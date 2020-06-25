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
import networkx as nx
if path.isdir("src"):
    from src.plotting_functions import plot_setting
else:
    from plotting_functions import plot_setting
from random import random, gauss, expovariate, shuffle
import time
from itertools import count
import statistics
from collections import defaultdict
from scipy.spatial import distance
from copy import deepcopy
import math
import numpy as np
from pylab import arange
from scipy.ndimage import measurements
import matplotlib.animation as animation

# Automatically setting the local path to this repo for easy file writing and saving
dir_path = path.dirname(path.realpath(__file__))

class Lattice():
    def __init__(self, size=(10, 10),
                 rand_dist=('uniform',),
                 torus_mode=True,
                 neighbourhood='vonNeumann',
                 distance='euclidean',
                 free_percent=0.1,
                 mutate_chance=0.5,
                 iterations=2000,
                 age_fraction=1 / 5):
        """
        Creates the Graph
        :param size: if type is a 2d graph size needs to be tuple, if type= grid_graph size is a list []
        In the model lattice itself the Bak-Snappen model is simulated collecting all the attributes such as the
        avalanche time, distance between mutations. The neighbourhood can be chosen between vonNeumann and Moore.
        The distance can be chosen between networkx and euclidean.
        """

        # Characteristics of the Network
        self.size = size
        self.lattice = nx.grid_graph(list(size), periodic=torus_mode)
        self.random_dist, *self.random_dist_specification = rand_dist
        self.neighbourhood = neighbourhood
        self.distance_btw_neighbours = distance  # either "networkx" or "euclidean"
        self.free_percent = free_percent
        self.mutate_chance = mutate_chance
        self.iterations = iterations
        self.age_fraction = age_fraction

        # initialising counters and single variable
        self.time_step = 0
        self.old_min_value = -1
        self.avalanche_timer = 1
        # get an numpy array representation of the lattice
        self.array = np.zeros((size[0], size[1]))

        # collectors
        self.min_value_list = []
        self.threshold_list = defaultdict(list)
        self.average_fit_list = []
        self.average_age_list = []
        self.avalanche_time_list = defaultdict(list)
        self.cluster_size = defaultdict(dict)
        self.distance_btw_mutation_list = []
        self.neighbours_list = []
        self.latest_mutation_pos_list = []

        # Collects the nodes that have a fitness/age/free
        self.fitness_dict = {}
        self.age_dict = {}
        self.free_dict = {}
        self.colour_dict = {}

        # check for error
        self.check_error()

    def fitness_init(self, node):
        """
        Initialize the fitness values to the graph using 3 different random distribution
        uniform, exponential and gaussian
        """
        self.lattice.nodes[node]['fitness'] = random()

    def age_init(self, node):
        """
        Assigns an age to the nodes which is defined by the last time it underwent mutation
        """
        self.lattice.nodes[node]['age'] = 0

    def get_nodes_w_age(self):
        """
        Gets all the nodes that have attribute age
        """
        self.age_dict = nx.get_node_attributes(self.lattice, 'age')

    def get_nodes_w_fitness(self):
        """
        Gets all the nodes that have attribute fitness
        """
        self.fitness_dict = nx.get_node_attributes(self.lattice, 'fitness')

    def get_nodes_w_is_free(self):
        """
        Gets the attribute free of all the nodes since all of them
        have this attribute
        """
        self.free_dict = nx.get_node_attributes(self.lattice, 'is_free')

    def get_nodes_w_colour(self):
        """
        Gets the normalized age in order to get the groups
        """
        self.colour_dict = nx.get_node_attributes(self.lattice, 'colour')

    def free_init(self):
        """
        Assigns if the node is free or not and add the fitness and age to it
        It is calling function age_init and fitness_init
        """
        for node in self.lattice.nodes:
            if random() > self.free_percent:
                self.lattice.nodes[node]['is_free'] = False
                self.fitness_init(node)
                self.age_init(node)
            else:
                self.lattice.nodes[node]['is_free'] = True

    def update_age(self):
        """
        Update the age at each time step apart from the positions in the list pos
        :param pos: Pos is a list of cells that have been update so the respective nodes are set to zero
        """
        for node in self.age_dict.keys():
            if node in self.latest_mutation_pos_list:
                self.lattice.nodes[node]['age'] = 0
            elif node not in self.latest_mutation_pos_list:
                self.lattice.nodes[node]['age'] += 1

    def get_min(self):
        """
        Get the minimum fitness value and its position, add it to collector and
        add it to threshold list if it is a new max min fitness
        """
        self.min_pos, self.min_value = min(self.fitness_dict.items(), key=lambda x: x[1])

        # Add min value to collector
        self.min_value_list.append(self.min_value)

        # check whatever the list is empty if not append only new maximum threshold
        if not self.threshold_list:
            self.threshold_list['threshold'].append(self.min_value)
            self.threshold_list['time_step'].append(self.time_step)
        elif self.min_value > max(self.threshold_list['threshold']):
            self.threshold_list['threshold'].append(self.min_value)
            self.threshold_list['time_step'].append(self.time_step)

    def get_average(self):
        """
        Get the average fitness/age value
        """
        # compute the mean
        self.average_fit = statistics.mean([self.fitness_dict[key] for key in self.fitness_dict])
        self.average_age = statistics.mean([self.age_dict[key] for key in self.age_dict])

        # Add average fitness at each time step to the collector
        self.average_fit_list.append(self.average_fit)
        self.average_age_list.append(self.average_age)

    def get_neighbours(self, chosen_node):
        """
        Get the neighbours of the given node and return self.neighbours which is a list of tuples
        """
        if self.neighbourhood == 'vonNeumann':
            neighbours_list = list(self.lattice.neighbors(chosen_node))
            return neighbours_list

        elif self.neighbourhood == 'Moore':
            neighbours_list = list()
            # Calculate the neighbours for this object
            for x1 in range(-1, 2):
                for y1 in range(-1, 2):
                    # Do not loop over yourself
                    if (x1, y1) != (0, 0):
                        x2 = (chosen_node[0] + x1)
                        y2 = (chosen_node[1] + y1) % (self.size[0] - 1)
                        if y2 >= 0 and y2 <= self.size[0] - 1:
                            if x2 >= 0 and x2 <= self.size[0] - 1:
                                neighbours_list.append((x2, y2))
                            else:
                                if x2 == -1:
                                    x2 = self.size[0] - 1
                                if x2 == self.size[0]:
                                    x2 = 0
                                neighbours_list.append((x2, y2))

            return neighbours_list

    def mutation(self):
        """
        Mutates the position with the lowest fitness and its neighbours
        """
        # Make sure the latest mutation list is empty every run
        # and add the ones with the minimum fitness if it is not migrating otherwise
        self.latest_mutation_pos_list = []
        self.latest_mutation_pos_list.append(self.min_pos)
        # Mutate the one with lowest fitness
        if self.random_dist == 'uniform':
            self.lattice.nodes[self.min_pos]['fitness'] = random()
            # Mutate the neighbours
            for node in self.neighbours_list:
                if node in self.fitness_dict.keys():
                    # check if the node has the attribute fitness
                    self.lattice.nodes[node]['fitness'] = random()
                    self.latest_mutation_pos_list.append(node)

        elif self.random_dist == 'exponential':
            self.lattice.nodes[self.min_pos]['fitness'] = expovariate(self.random_dist_specification[0])
            # Mutate the neighbours
            for node in self.neighbours_list:
                # check if the node has the fitness attribute
                if node in self.fitness_dict.keys():
                    self.lattice.nodes[node]['fitness'] = expovariate(self.random_dist_specification[0])
                    self.latest_mutation_pos_list.append(node)

        elif self.random_dist == 'gauss':
            gaussian_random = gauss(self.lattice.nodes[self.min_pos]['fitness'], self.random_dist_specification[1])
            if gaussian_random >= 0 and gaussian_random <= 1:
                self.lattice.nodes[self.min_pos]['fitness'] = gaussian_random
            # Mutate the neighbours
            for node in self.neighbours_list:
                # check if the node has the attribute fitness
                if node in self.fitness_dict.keys():
                    gaussian_random = gauss(self.lattice.nodes[node]['fitness'], self.random_dist_specification[1])
                    if gaussian_random >= 0 and gaussian_random <= 1:
                        self.lattice.nodes[node]['fitness'] = gaussian_random
                    self.latest_mutation_pos_list.append(node)

    def move(self):
        """
        The node with lowest fitness get moved to free place when the mean fitness of its neighbors is above the total average while lower than new neighbors
        """
        free_nodes = list(dict(filter(lambda elem: elem[1], self.free_dict.items())).keys())
        shuffle(free_nodes)
        neighbours_fitness_list = []

        for neighbour in self.neighbours_list:
            if not self.lattice.nodes(data=True)[neighbour]['is_free']:
                neighbours_fitness_list.append(self.lattice.nodes(data=True)[neighbour]['fitness'])

        if neighbours_fitness_list and statistics.mean(neighbours_fitness_list) < self.average_fit:
            for free_node in free_nodes:
                neighbours_fitness_list = []
                for neighbour in self.get_neighbours(free_node):
                    if not self.lattice.nodes(data=True)[neighbour]['is_free']:
                        neighbours_fitness_list.append(self.lattice.nodes(data=True)[neighbour]['fitness'])

                if neighbours_fitness_list and statistics.mean(neighbours_fitness_list) > self.average_fit:
                    self.lattice.nodes[free_node].update(self.lattice.nodes[self.min_pos])

                    del self.lattice.nodes[self.min_pos]['fitness']
                    del self.lattice.nodes[self.min_pos]['age']
                    self.lattice.nodes[self.min_pos]['is_free'] = True
                    self.min_pos = free_node

                    self.get_nodes_w_fitness()
                    self.get_nodes_w_age()
                    self.get_nodes_w_is_free()
                    break
        self.mutation()

    def move_or_mutate(self):
        """
        Choosing between mutating the node with the lowest fitness or moving the node to free spaces
        """
        if random() < self.mutate_chance:
            self.mutation()
        else:
            self.move()

    def get_avalanche_time(self):
        """
        Counts the number of iteration need until the threshold is rised to a new max
        It appends to the list avalanche time list, the time step where a avalanche stops with the respective steps of
        the avalanche
        """
        # check if all the new mutation are above the max min fintess so if they rised the threshold
        if not all(self.lattice.nodes[node]['fitness'] > max(self.threshold_list['threshold']) for node in
                   self.latest_mutation_pos_list):
            self.avalanche_timer += 1
        elif all(self.lattice.nodes[node]['fitness'] > max(self.threshold_list['threshold']) for node in
                 self.latest_mutation_pos_list):
            self.avalanche_time_list['avalanche_time'].append(self.avalanche_timer)
            self.avalanche_time_list['time_step'].append(self.time_step)
            # rest avalanche_timer
            self.avalanche_timer = 1

    def get_dist_btw_mutation(self):
        """
        Get the distance between minimum fintess nodes and the next minimum fitness nodes
        Using the build in function the distances is compute by how many nodes to go  before reaching the target
        So if  source is 0,0 and target is 1,1 , the lenght will be 2 and not as in euclidean space 1.4
        So both networkx and euclidean method are implemented
        """
        if self.old_min_value != -1:
            if self.distance_btw_neighbours == 'networkx':
                dist = nx.shortest_path_length(self.lattice, self.min_pos, self.old_min_value)
                self.distance_btw_mutation_list.append(dist)
            elif self.distance_btw_neighbours == 'euclidean':
                dist = distance.euclidean(self.min_pos, self.old_min_value)
                self.distance_btw_mutation_list.append(dist)

        elif self.old_min_value == -1:
            # Assing the first value to self.old_min_value
            self.distance_btw_mutation_list.append(0)

        # Save current min_pos to compare to the next new min_pos
        self.old_min_value = deepcopy(self.min_pos)

    def run(self, things_to_collect):
        """
        Run the Bak-Sneppen model using the different rules
        Important to note that the order needs to be in this respective way
        For example avalanche time can work before mutation
        """
        # initialize the nodes with fitness and their age
        self.free_init()

        # get all the nodes that have explicitly attribute age/fitness
        self.get_nodes_w_fitness()
        self.get_nodes_w_age()
        self.get_nodes_w_is_free()

        for i in range(self.iterations):
            # nice if animated
            # plt.figure()
            # self.plot()

            # get the nodes with the minimum value
            self.get_min()

            # gets the average fitness/age of all the nodes at each time_step
            self.get_average()

            # get the neighbours
            self.neighbours_list = self.get_neighbours(self.min_pos)

            # assign new random number to the lowest fitness and its neighbours
            if "mutation" in things_to_collect or "all" in things_to_collect:
                self.mutation()

            if "moving" in things_to_collect or "all" in things_to_collect:
                self.move_or_mutate()

            # check if new mutation rise the threshold
            if "avalanche_time" in things_to_collect or "all" in things_to_collect:
                self.get_avalanche_time()

            # set the age of the nodes accordingly and needs to be placed after self.get_avalanche_time()
            if "update_age" in things_to_collect or "all" in things_to_collect:
                self.update_age()

            # update the dicts after mutation
            self.get_nodes_w_fitness()
            self.get_nodes_w_age()
            self.get_nodes_w_is_free()

            # get the distance between mutations
            if "get_dist_btw_mutation" in things_to_collect or "all" in things_to_collect:
                self.get_dist_btw_mutation()

            # get clusters
            if "get_cluster" in things_to_collect or "all" in things_to_collect:
                self.get_clusters()

            self.time_step += 1

    def plot(self, label='fitness'):
        """
        Visualise the graph and plot labels it labels = True
        """
        if label == 'fitness':
            values = set(self.fitness_dict.values())
            mapping = dict(zip(sorted(values), count()))
            nodes_fitness = self.fitness_dict.keys()

            # Get the color code and normalise it
            colors = [mapping[self.lattice.nodes[n]['fitness']] for n in nodes_fitness]
            colors = [color / (self.size[0] * self.size[1]) for color in colors]

            pos = dict((n, n) for n in self.lattice.nodes())
            nc = nx.draw_networkx_nodes(self.lattice, pos, nodelist=nodes_fitness, node_color=colors,
                                        with_labels=False, node_size=200, node_shape='h', cmap=plt.cm.jet)
            # plt.figure()
            # plt.colorbar(nc)
            # plt.axis('off')

        elif label == 'age':
            values = set(self.age_dict.values())
            mapping = dict(zip(sorted(values), count()))
            nodes_ages = self.age_dict.keys()

            # Get the color code and normalise it
            colors = [mapping[self.lattice.nodes[n]['age']] for n in nodes_ages]

            pos = dict((n, n) for n in self.lattice.nodes())
            nc = nx.draw_networkx_nodes(self.lattice, pos, nodelist=nodes_ages, node_color=colors,
                                        with_labels=False, node_size=200, node_shape='h', cmap=plt.cm.gnuplot)
            # plt.figure()
            # plt.colorbar(nc)
            # plt.axis('off')

        return nc

    def check_error(self):
        # Get possible errors
        if self.random_dist == ' exponential' and len(self.random_dist_specification) != 1:
            raise Exception("If distribution is set to exponential, 1 argument must be added !! "
                            "rand_dist=('exponential',1)")
        elif self.random_dist == 'gauss' and len(self.random_dist_specification) != 2:
            raise Exception("If distribution is set to gauss, 2 argument must be added (mean,std) !! "
                            "rand_dist=('gauss',mean,std)")

    def cluster_colour(self, node):
        """
        Introduces "colour" to the nodes for the cluster analysis.
        We define ranges of ages for the clustering.
        The variable "age_fraction" defines the size of each age range
        as a fraction.
        """

        # define max & min age
        max_age = max(self.age_dict.values())
        min_age = min(self.age_dict.values())
        # each cluster has an age range given by:
        cluster_age_range = (max_age - min_age) * self.age_fraction

        # make sure that in the beginning the number of groups is min 1
        if cluster_age_range == 0:
            cluster_age_range = 1

        self.lattice.nodes[node]['colour'] = math.ceil(self.lattice.nodes[node]['age'] / cluster_age_range)

    def draw_array(self, group_nr):
        """
        Redraw the array with the new normalized ages
        """
        for coord in self.colour_dict.keys():
            if self.colour_dict[coord] == group_nr:
                self.array[coord[0], coord[1]] = 1

    def reset_array(self):
        """
        Resets the array to zero since the function can't deal with it otherwise
        """
        self.array = np.zeros((self.size[0], self.size[1]))

    def get_clusters(self):
        """
        Get clusters and respective sizes
        """
        # first regroup similar ages to same value
        for node in self.age_dict.keys():
            self.cluster_colour(node)

        # collect in a dict all the normalised ages
        self.get_nodes_w_colour()

        # get the different groups all the different groups
        groups = np.unique(list(self.colour_dict.values()))

        # combine all the cluster sizes per iteration
        area_dist_per_itr = []

        # draw the array
        for group in groups:
            self.draw_array(group_nr=group)
            lw, num = measurements.label(self.array)
            area = measurements.sum(self.array, lw, index=arange(lw.max() + 1))
            # make sure to not include a zero in the array
            area = area[area != 0]
            #area_dist_per_itr = area_dist_per_itr + (list(area))
            # make sure to reset the array
            self.reset_array()

        # append it to the collector dict where each key corresponds to the time_step
            self.cluster_size[self.time_step][group] = area


if __name__ == "__main__":

    plot = True
    iterations = 2000
    t0 = time.time()
    # if rand_dist take 1 arg, rand_dist=('uniform',) !! Comma needed here
    lattice = Lattice(size=(20, 20), torus_mode=True, rand_dist=('uniform',),
                      free_percent=0, iterations=iterations, age_fraction=1 / 10)
    print(nx.info(lattice.lattice, n=None))
    lattice.run(["all"])
    t1 = time.time()

    print("The average fitness is {}".format(lattice.average_fit_list[-1]))
    print("TOTAL TIME NEEDED {}".format(t1 - t0))

    if plot:
        # make sure the default parameters are the same
        plot_setting()
        # The plotting need to be fixed
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
        ax1.plot(lattice.min_value_list, label='min_value')
        ax1.plot(lattice.average_fit_list, label='average_fitness')
        ax1.legend()
        ax1.set_title('Average fitness and maximum minimum fit')

        n, bins = np.histogram(lattice.avalanche_time_list['avalanche_time'], density=True)
        ax2.set_title('Avalanche size')
        ax2.set_yscale('log')
        ax2.set_xscale('log')
        ax2.plot(bins[:-1], n)

        n, bins = np.histogram(lattice.distance_btw_mutation_list, density=True)
        ax3.plot(bins[:-1], n, 100)
        ax3.set_title('The Distribution of the distances between mutations')
        ax3.set_yscale('log')
        ax3.set_xscale('log')

        plt.figure()
        plt.title('Random original fitness distribution on a 2D lattice')
        fitness = lattice.plot(label='fitness')
        plt.colorbar(fitness)
        plt.tight_layout()
        plt.savefig(path.join(dir_path, 'figures/lattice_itr=1.png'), dpi=300)
        
        plt.figure()
        plt.title('Age Distribution after {} iterations'.format(iterations))
        age = lattice.plot(label='age')
        plt.colorbar(age)
        plt.tight_layout()
        plt.savefig(path.join(dir_path, 'figures/lattice-age_itr={}.png'.format(iterations)), dpi=300)
        #

        number_of_frames = 2000

        def update_hist(num):
            plt.cla()
            plt.hist(np.concatenate([lattice.cluster_size[num][x] for x in lattice.cluster_size[num]]),bins=50)


        fig = plt.figure()
        hist = plt.hist(np.concatenate([lattice.cluster_size[0][x] for x in lattice.cluster_size[0]]),bins=50)

        animation = animation.FuncAnimation(fig, update_hist, number_of_frames,interval=20)
        plt.show()