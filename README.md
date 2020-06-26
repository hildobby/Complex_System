# Expanding and exploring the Bak-Sneppen model written in Python 3+
=========================================
2 dimensions and 3 dimensions Bak-Snappen model, evaluation of the resulting model behaviour and effect of various parameters.

By: [Hildebert Moulié](https://github.com/hildobby), [Philippe Nicolau](https://github.com/PNicolau96), [Louis Weyland](https://github.com/LouisWW) & [Binjie Zhou](https://github.com/binjiezhou).

![Initial state of the lattice fitnessses](https://github.com/hildobby/Complex_System/src/figures/lattice_itr=1.png) ![Age of the nodes on the lattice after 2000 iterations](https://github.com/hildobby/Complex_System/src/figures/lattice-age_itr=2000.png)

## Setup required packages

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

An overview of the results can also be found in [the Jupyter notebook file](https://github.com/hildobby/Complex_System/blob/master/results_overview.ipynb) named `results_overview.ipynb`.

## Project Plan

1. Bring the 1 dimension Bak-Snappen model to 2 dimensions and 3 dimensions
2. Evaluate different results from the model:
    * Avalanche time and compare it between several models
    * Mutation distance (Check if power law still applies here)
    * Cluster sizes
    * Evolutionary time
3. Test the data with the [powerlaw](https://pypi.org/project/powerlaw/) package



## Overview of some of the results




## References

* Bak, P., & Sneppen, K. (1993). Punctuated equilibrium and criticality in a simple model of evolution. Physical review letters, 71(24), 4083.
* Fraiman, D. (2018). Bak-Sneppen model: Local equilibrium and critical value. Physical Review E, 97(4), 042123.
* Paczuski, M., Maslov, S., & Bak, P. (1996). Avalanche dynamics in evolution, growth, and depinning models. Physical Review E, 53(1), 414.
* Alstott, J., Bullmore, E., & Plenz, D. (2014). powerlaw: a Python package for analysis of heavy-tailed distributions. PloS one, 9(1), e85777.
