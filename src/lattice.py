#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Tuesday June 16 2020
This code was implemented by
Louis Weyland, Hildebert Mouilé, Philippe Nicolau & Binjie Zhou.
"""


import matplotlib.pyplot as plt
import matplotlib.animation
import networkx as nx
from plotting_functions import plot_setting
from random import random, gauss, expovariate
import time
from itertools import count
import statistics
from collections import defaultdict
from scipy.spatial import distance
from copy import deepcopy



class Lattice():

    def __init__(self,size=(10,10,2),
                 rand_dist=('uniform',),
                 torus_mode=True,
                 neighbourhood='vonNeumann',
                 distance='euclidean'):
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
        self.random_dist,*self.random_dist_specification = rand_dist
        self.neighbourhood = neighbourhood
        self.distance_btw_neighbours = distance   # either "networkx" or "euclidean"

        # initialising counters and single variable
        self.time_step = 0
        self.old_min_value = -1
        self.avalanche_timer = 1

        # collectors
        self.min_value_list = []
        self.threshold_list = []
        self.average_fit_list = []
        self.average_age_list = []
        self.avalanche_time_list = defaultdict(list)
        self.distance_btw_mutation_list = []


        # check for error
        self.check_error()


    def fitness_init(self):
        """
        Initialize the fitness values to the graph using 3 different random distribution
        uniform, exponential and gaussian
        """
        if self.random_dist == 'uniform':
            for node in self.lattice.nodes:
                self.lattice.nodes[node]['fitness'] = random()
        elif self.random_dist == 'exponential':
            for node in self.lattice.nodes:
                self.lattice.nodes[node]['fitness'] = expovariate(self.random_dist_specification[0])
        elif self.random_dist == 'gauss':
            for node in self.lattice.nodes:
                self.lattice.nodes[node]['fitness'] = gauss(self.random_dist_specification[0],
                                                                  self.random_dist_specification[1])

    def age_init(self):
        """
        Assigns an age to the nodes which is defined by the last time it underwent mutation
        """
        for node in self.lattice.nodes:
            self.lattice.nodes[node]['age'] = 0

    def update_age(self):
        """
        Update the age at each time step appart from the positions in the list pos
        :param pos: Pos is a list of cells that have been update so the respective nodes are set to zero
        """
        for node in self.lattice.nodes:
            if node in self.latest_mutation_pos_list:
                self.lattice.nodes[node]['age'] = 0
            elif node not in self.latest_mutation_pos_list:
                self.lattice.nodes[node]['age'] += 1



    def get_min(self):
        """
        Get the minimum fitness value and its position, add it to collector and
        add it to threshold list if it is a new max min fitness
        """
        min_dict= nx.get_node_attributes(self.lattice, 'fitness')
        self.min_pos,self.min_value = min(min_dict.items(), key=lambda x: x[1])

        # Add min value to collector
        self.min_value_list.append(self.min_value)

        # check whatever the list is empty if not append only new maximum threshold
        if not self.threshold_list:
            self.threshold_list.append(self.min_value)
        elif self.min_value > max(self.threshold_list):
            self.threshold_list.append(self.min_value)

    def get_avergae(self):
        """
        Get the average fitness/age value
        """
        fitness_dict = nx.get_node_attributes(self.lattice, 'fitness')
        age_dict = nx.get_node_attributes(self.lattice,'age')
        # compute the mean
        self.average_fit = statistics.mean([fitness_dict[key] for key in fitness_dict])
        self.average_age = statistics.mean([age_dict[key] for key in age_dict])

        # Add average fitness at each time step to the collector
        self.average_fit_list.append(self.average_fit)
        self.average_age_list.append(self.average_age)

    def get_neighbours(self):
        """
        Get the neighbours of the lowest fitness and return self.neighbours which is a list of tuples
        """
        if self.neighbourhood == 'vonNeumann':
            self.neighbours = list(self.lattice.neighbors(self.min_pos))

        elif self.neighbourhood == 'Moore':

            raise Exception('Need to implement the Moor neighborhood !')

    def mutation(self):
        """
        Mutates the position with the lowest fitness and its neighbours
        and check the avalanche time is the max of the lowest fitness was increased
        """
        # Mutate the one with lowest fitness
        if self.random_dist == 'uniform':
            self.lattice.nodes[self.min_pos]['fitness'] = random()
            # Mutate the neighbours
            for node in self.neighbours:
                self.lattice.nodes[node]['fitness'] = random()
        elif self.random_dist == 'exponential':
            self.lattice.nodes[self.min_pos]['fitness'] = expovariate(self.random_dist_specification[0])
            # Mutate the neighbours
            for node in self.neighbours:
                self.lattice.nodes[node]['fitness'] = expovariate(self.random_dist_specification[0])
        elif self.random_dist == 'gauss':
            self.lattice.nodes[self.min_pos]['fitness'] = gauss(self.random_dist_specification[0],
                                                                self.random_dist_specification[1])
            # Mutate the neighbours
            for node in self.neighbours:
                self.lattice.nodes[node]['fitness'] = gauss(self.random_dist_specification[0],
                                                            self.random_dist_specification[1])

        # check if new mutation rised the threshold
        self.get_avalanche_time()

        # set the age of the nodes accordingly and needs to be placed after self.get_avalanche_time()
        self.update_age()

    def get_avalanche_time(self):
        """
        Counts the number of iteration need until the threshold is rised to a new max
        It appends to the list avalanche time list, the time step where a avalanche stops with the respective steps of
        the avalanche
        """
        # combine the node with the lowest fitness and its neighbours
        self.latest_mutation_pos_list = []
        self.latest_mutation_pos_list = deepcopy(self.neighbours)
        self.latest_mutation_pos_list.append(self.min_pos)

        # check if all the new mutation are above the max min fintess so if they rised the threshold
        if not all(self.lattice.nodes[node]['fitness'] > max(self.threshold_list) for node in  self.latest_mutation_pos_list):
            self.avalanche_timer += 1
        elif all(self.lattice.nodes[node]['fitness'] > max(self.threshold_list) for node in  self.latest_mutation_pos_list):
            self.avalanche_time_list['avalanche time'].append(self.avalanche_timer)
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

        elif self.old_min_value == -1 :
            # Assing the first value to self.old_min_value
            self.distance_btw_mutation_list.append(0)

        # Save current min_pos to compare to the next new min_pos
        self.old_min_value = deepcopy(self.min_pos)

    def run(self,iteration):
        """
        Run the Bak-Sneppen model using the different rules
        """
        # initialize the nodes with fitness and their age
        self.fitness_init()
        self.age_init()

        for i in range(iteration):

            # nice if animated
            #plt.figure()
            #self.plot()

            # get the nodes with the minimum vale
            self.get_min()

            self.get_avergae()

            # get the neighbours
            self.get_neighbours()
            # assign new random number to the lowest fitness and its neighbours
            self.mutation()
            # get the distance between mutations
            self.get_dist_btw_mutation()

            self.time_step += 1

        # The plotting need to be fixed
        plt.figure()
        plt.plot(self.min_value_list)
        plt.plot(self.average_fit_list)
        plt.plot(self.avalanche_time_list['time_step'],self.avalanche_time_list['avalanche time'])
        print("The average fitness is {}".format(self.average_fit_list[-1]))
        plt.show()

    def plot(self):
        """
        Visualise the graph and plot labels it labels = True
        """
        values = set(nx.get_node_attributes(self.lattice, 'fitness').values())
        mapping = dict(zip(sorted(values), count()))
        nodes = self.lattice.nodes()

        # Get the color code and normalise it
        colors = [mapping[self.lattice.nodes[n]['fitness']] for n in nodes]
        colors = [color/(self.size[0]*self.size[1]) for color in colors]

        pos = dict( (n, n) for n in self.lattice.nodes() )
        nc = nx.draw_networkx_nodes(self.lattice, pos, nodelist=nodes, node_color=colors,
                                    with_labels=False, node_size=100, cmap=plt.cm.jet)
        plt.colorbar(nc)
        plt.axis('off')
        plt.show()
        return nc

    def check_error(self):
        # Get possible errors
        if self.random_dist == ' exponential' and len(self.random_dist_specification) != 1:
            raise Exception("If distribution is set to exponential, 1 argument must be added !! "
                            "rand_dist=('exponential',1)")
        elif self.random_dist == 'gauss' and len(self.random_dist_specification) != 2:
            raise Exception("If distribution is set to gauss, 2 argument must be added (mean,std) !! "
                            "rand_dist=('gauss',mean,std)")


if __name__ == "__main__":
    t0 = time.time()
    # if rand_dist take 1 arg, rand_dist=('uniform',) !! Comma needed here
    lattice = Lattice(size=(20,20),torus_mode=True,rand_dist=('uniform',))
    print(nx.info(lattice.lattice, n=None))
    lattice.run(iteration=3000)
    t1 = time.time()

    print("TOTAL TIME NEEDED {}".format(t1-t0))

    plt.figure()
    plt.hist(lattice.distance_btw_mutation_list)


