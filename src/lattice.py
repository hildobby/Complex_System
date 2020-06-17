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
from random import random, gauss
import time
from itertools import count



class Lattice():

    def __init__(self,size=(10,10,1),rand_dist=random,torus_mode=True):
        """
        Creates the Graph
        :param size: if type is a 2d graph size needs to be tuple, if type= grid_graph size is a list []
        """
        self.type = type
        self.size = size
        self.lattice = nx.grid_graph(list(size), periodic=torus_mode)
        self.random_dist = rand_dist
        self.time_step = 0
        self.threshold = []


    def random_init(self):
        """
        Initialize the fitness values to the graph
        """
        for node in self.lattice.nodes:
            self.lattice.nodes[node]['fitness'] = self.random_dist()


    def get_min(self):
        """
        Get the minimum fitness value and its position
        """
        min_dict= nx.get_node_attributes(self.lattice, 'fitness')
        self.min_pos,self.min_value = min(min_dict.items(), key=lambda x: x[1])
        #print("The position with lowest fintess of {} is {}".format(self.min_value,self.min_pos))


    def get_neighbours(self):
        """
        Get the neighbours of the lowest fitness
        """
        self.neighbours = list(self.lattice.neighbors(self.min_pos))

        # Add if not Neumann get neighbours


    def mutation(self):
        """
        Mutates the position with the lowest fitness and its neighbours
        """
        # Mutate the one with lowest fitness
        self.lattice.nodes[self.min_pos]['fitness'] = self.random_dist()
        # Mutate the neighbours
        for node in (self.neighbours):
            self.lattice.nodes[node]['fitness'] = self.random_dist()


    def run(self,iteration):
        """
        Run the Bak-Sneppen model using the different rules
        """
        # initialize the nodes with random values
        self.random_init()

        for i in range(iteration):

            # nice if animated
            self.plot()

            # get the nodes with the minimum vale
            self.get_min()

            self.threshold.append(self.min_value)
            # get the neigbours
            self.get_neighbours()
            # assign new random number to the lowest fitness and its neighbours
            self.mutation()

            self.time_step += 1

        plt.plot(self.threshold)
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


if __name__ == "__main__":

    lattice = Lattice(size=(40,40),torus_mode=False)
    t0 = time.time()
    lattice.run(iteration=8)
    t1 = time.time()
    print("Total time needed is {}".format(t1 - t0))


