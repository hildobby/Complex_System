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
                 distance='euclidean',
                 free_percent=0.1):
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
        self.free_percent = free_percent


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


    def fitness_init(self,node):
        """
        Initialize the fitness values to the graph using 3 different random distribution
        uniform, exponential and gaussian
        """
        self.lattice.nodes[node]['fitness'] = random()


    def age_init(self,node):
        """
        Assigns an age to the nodes which is defined by the last time it underwent mutation
        """
        self.lattice.nodes[node]['age'] = 0


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
        Update the age at each time step appart from the positions in the list pos
        :param pos: Pos is a list of cells that have been update so the respective nodes are set to zero
        """
        for node in nx.get_node_attributes(lattice.lattice,'age').keys():
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

    def get_average(self):
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
            if self.min_pos[0]==0:
                self.neighbours.append((self.size[0]-1,self.min_pos[1]))
            elif self.min_pos[0]==self.size[0]-1:
                #print((0,self.min_pos[1]))
                self.neighbours.append((0,self.min_pos[1]))
           # print(self.neighbours)
            
        elif self.neighbordhood ==  'Moore':
            #print(self.min_pos)
            self.neighbours = list()
            # Calculate the neighbours for this object
            for x1 in range(-1,2):
                for y1 in range(-1,2):
                    # Do not loop over yourself
                    if (x1,y1)!=(0,0):
                        x2 = (self.min_pos[0]+x1)
                        y2 = (self.min_pos[1]+y1) % (self.size[0]-1)
                        if y2>=0 and y2<=self.size[0]-1:
                            if x2>=0 and x2<=self.size[0]-1:
                                self.neighbours.append((x2,y2))
                            else:
                                if x2 == -1:
                                    x2=self.size[0]-1
                                if x2 == self.size[0]:
                                    x2=0
                                self.neighbours.append((x2,y2))

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
        Important to note that the order needs to be in this respective way
        For example avalanche time can work before mutation
        """
        # initialize the nodes with fitness and their age
        self.free_init()

        for i in range(iteration):

            # nice if animated
            #plt.figure()
            #self.plot()

            # get the nodes with the minimum vale
            self.get_min()

            # gets the average fitness/age of all the nodes at each time_step
            self.get_average()

            # get the neighbours
            self.get_neighbours()
            # assign new random number to the lowest fitness and its neighbours
            self.mutation()

            # check if new mutation rised the threshold
            self.get_avalanche_time()

            # set the age of the nodes accordingly and needs to be placed after self.get_avalanche_time()
            self.update_age()
            # get the distance between mutations
            self.get_dist_btw_mutation()

            self.time_step += 1


    def plot(self,label= 'fitness'):
        """
        Visualise the graph and plot labels it labels = True
        """
        if label == 'fitness':
            values = set(nx.get_node_attributes(self.lattice, 'fitness').values())
            mapping = dict(zip(sorted(values), count()))
            nodes_fitness = nx.get_node_attributes(lattice.lattice,'fitness').keys()

            # Get the color code and normalise it
            colors = [mapping[self.lattice.nodes[n]['fitness']] for n in nodes_fitness]
            colors = [color/(self.size[0]*self.size[1]) for color in colors]

            pos = dict( (n, n) for n in self.lattice.nodes() )
            nc = nx.draw_networkx_nodes(self.lattice, pos, nodelist=nodes_fitness, node_color=colors,
                                        with_labels=False, node_size=200,node_shape = 'h', cmap=plt.cm.jet)
            #plt.figure()
            #plt.colorbar(nc)
            #plt.axis('off')

        elif label == 'age':
            values = set(nx.get_node_attributes(self.lattice, 'age').values())
            mapping = dict(zip(sorted(values), count()))
            nodes_ages = nx.get_node_attributes(lattice.lattice,'age').keys()

            # Get the color code and normalise it
            colors = [mapping[self.lattice.nodes[n]['age']] for n in nodes_ages]

            pos = dict((n, n) for n in self.lattice.nodes())
            nc = nx.draw_networkx_nodes(self.lattice, pos, nodelist=nodes_ages, node_color=colors,
                                        with_labels=False, node_size=200,node_shape = 'h', cmap=plt.cm.gnuplot)
            #plt.figure()
            #plt.colorbar(nc)
            #plt.axis('off')

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

    plot=True
    iterations = 2000
    t0 = time.time()
    # if rand_dist take 1 arg, rand_dist=('uniform',) !! Comma needed here
    lattice = Lattice(size=(20,20),torus_mode=True,rand_dist=('uniform',),free_percent=0.1)
    print(nx.info(lattice.lattice, n=None))
    lattice.run(iteration=iterations)
    t1 = time.time()
    print("The average fitness is {}".format(lattice.average_fit_list[-1]))
    print("TOTAL TIME NEEDED {}".format(t1-t0))

    if plot:
        # make sure the default parameters are the same
        plot_setting()
        # The plotting need to be fixed
        fig, (ax1, ax2, ax3) = plt.subplots(1,3)
        ax1.plot(lattice.min_value_list,label='min_value')
        ax1.plot(lattice.average_fit_list,label='average_fitness')
        ax1.legend()
        ax1.set_title('Average fitness and maximum minimum fit')

        ax2.plot(lattice.avalanche_time_list['time_step'], lattice.avalanche_time_list['avalanche time'])
        ax2.set_title('Avalanche size over time')

        ax3.hist(lattice.distance_btw_mutation_list)
        ax3.set_title('The Distribution of the distances between mutations')


        plt.figure()
        plt.title('Fitness Distribution')
        fitness = lattice.plot(label='fitness')
        plt.colorbar(fitness)

        plt.figure()
        plt.title('Age Distribution after {} itr'.format(iterations))
        age = lattice.plot(label='age')
        plt.colorbar(age)


        plt.show()


