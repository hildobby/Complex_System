#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Tuesday June 16 2020
This code was implemented by
Louis Weyland, Hildebert Mouilé, Philippe Nicolau & Binjie Zhou.
"""


import matplotlib.pyplot as plt
import networkx as nx
from src.plotting_functions import plot_setting
from random import random, gauss
import time
from numba import jit


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


    def random_init(self):
        """
        Initialize the fitness values to the graph
        """
        for node in self.lattice.nodes:
            self.lattice.nodes[node]['fitness'] = random()

    def plot(self,show_labels=False):
        """
        Visualise the graph and plot labels it labels = True
        """

        if show_labels == True:
            labels = nx.get_node_attributes(self.lattice, 'fitness')
            nx.draw(self.lattice, labels=labels, node_size=20)
        else:
            nx.draw(self.lattice, node_size=20)

        plot_setting()
        plt.show()
        plt.close()


if __name__ == "__main__":

    lattice = Lattice(size=(10,15,1),torus_mode=False)

    t0 = time.time()
    lattice.random_init()

    t1 = time.time()
    print("Total time needed is {}".format(t1 - t0))
    lattice.plot()

