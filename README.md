Expanding and exploring the Bak-Sneppen model written in Python 3+
=========================================

2 dimensions and 3 dimensions Bak-Snappen model, evaluation of the resulting model behaviour and effect of various parameters.

By: [Hildebert Moulié](https://github.com/hildobby), [Philippe Nicolau](https://github.com/PNicolau96), [Louis Weyland](https://github.com/LouisWW) & [Binjie Zhou](https://github.com/binjiezhou).

<p float="left" class="center">
  <img src="https://github.com/hildobby/Complex_System/blob/master/src/figures/lattice_itr%3D1.png" width="350" />
  <img src="https://github.com/hildobby/Complex_System/blob/master/src/figures/lattice-age_itr%3D2000.png" width="350" /> 
</p>

## Project Plan

1. Bring the 1 dimension Bak-Snappen model to 2 dimensions and 3 dimensions
2. Evaluate different results from the model:
    * Avalanche time and compare it between several models
    * Mutation distance
    * Cluster sizes
    * Evolutionary time
3. Test the data with the [powerlaw](https://pypi.org/project/powerlaw/) package

## Overview of some of the results

### Fitness over time
![Fitness over time](https://github.com/hildobby/Complex_System/blob/master/src/figures/average_fitness_s%3D(20%2C%2020)_itr%3D2000_rep%3D10_std%3D0.3.png)

### Impact of percentage of empty node
![Impact of percentage of empty node](https://github.com/hildobby/Complex_System/blob/master/presentation_content/images/emptynode1.png)


### Cluster size distribution for different grid sizes
![Cluster size distribution for different grid sizes](https://github.com/hildobby/Complex_System/blob/master/src/figures/)

## Working with the repository

Gathering the required packages
```
pip3 install -r requirements.txt
```

Running the model
```
python3 main.py
```

Generating the presentation from the `results_presentation.ipynb` Jupyter notebook
```
jupyter nbconvert results_presentation.ipynb --to slides --post serve
```

## References

* Alstott, J., Bullmore, E., & Plenz, D. (2014). powerlaw: a Python package for analysis of heavy-tailed distributions. PloS one, 9(1), e85777.
* Bak, P., & Sneppen, K. (1993). Punctuated equilibrium and criticality in a simple model of evolution. Physical review letters, 71(24), 4083.
* Fraiman, D. (2018). Bak-Sneppen model: Local equilibrium and critical value. Physical Review E, 97(4), 042123.
* Paczuski, M., Maslov, S., & Bak, P. (1996). Avalanche dynamics in evolution, growth, and depinning models. Physical Review E, 53(1), 414.
