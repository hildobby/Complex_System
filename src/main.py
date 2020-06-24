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
    from plotting_functions import plot_setting
    from lattice.py import Lattice
else:
    from src.plotting_functions import plot_setting
    from src.lattice import Lattice

import time
import numpy as np

# Get a comparison between the different random distribution

iterations = 2000
t0 = time.time()
# if rand_dist take 1 arg, rand_dist=('uniform',) !! Comma needed here



repetition = 10

uniform_list = []
gaussian_list = []
exponential_list = []

for i in range(repetition):
    print("Repetiton {}".format(i))

    uniform = Lattice(size=(20, 20), torus_mode=True, rand_dist=('uniform',), free_percent=0, iterations=iterations,
                      age_fraction=1 / 10)
    gaussian = Lattice(size=(20, 20), torus_mode=True, rand_dist=('gaussian', 0.5, 0.2), free_percent=0,
                      iterations=iterations, age_fraction=1 / 10)
    exponential = Lattice(size=(20, 20), torus_mode=True, rand_dist=('exponential', 1), free_percent=0,
                          iterations=iterations, age_fraction=1 / 10)


    uniform.run(["mutation"])
    gaussian.run("mutation")
    exponential.run("mutation")

    uniform_list.append(uniform.average_fit_list)
    gaussian_list.append(gaussian.average_fit_list)
    exponential_list.append(exponential.average_fit_list)



t1 = time.time()
print("TOTAL TIME NEEDED {}".format(t1 - t0))