# Complex System Simulation: Expanding and exploring the Bak-Snappen model

2 Dimensions Bak-Snappen model, evaluation of the resulting model behaviour and effect of various parameters.

By: [Hildebert Moulié](https://github.com/hildobby), [Philippe Nicolau](https://github.com/PNicolau96), [Louis Weyland](https://github.com/LouisWW) & [Binjie Zhou](https://github.com/binjiezhou).

## Setup required packages

Gathering the required packages
```
pip3 install -r requirements.txt
```

Running the model
```
python3 main.py
```

Looking at the presentation
```
jupyter nbconvert results_presentation.ipynb --to slides --post serve
```


An overview of the results can also be found in [the Jupyter notebook file](https://github.com/hildobby/Complex_System/blob/master/results_overview.ipynb) named `results_overview.ipynb`.

## Project Plan

1. Bring the 1 dimension Bak-Snappen model to 2 dimensions
2. Evaluate different results from the model:
    * Avalanche time and compare it between several models
    * Mutation distance (Check if power still applies here)
    * Cluster sizes
    * Evolutionary time
    * Different types of neighbourood (Moore and Von Neumann)
3. Test with real provided data
4. If time remains, repeat with 3 dimensions

[Here](https://docs.google.com/document/d/1rTodhozVX6pXBGTlviCNNwwcX3oG-x80GXP4DIutzHM/edit?usp=sharing) you will find the link to the Google Document of original project plan.

## An overview of some of the results




## Useful resources:

* [Concise explanation of a simple 1 dimension Bak-Sneppen model](http://csmgeo.csm.jmu.edu/geollab/fichter/gs102/2004handouts/Bak-Sneppenbrief.PDF)
* [Repo with a 2 dimensions Bak-Sneppen model for inspiration](https://github.com/voschezang/Spatial-Bak-Sneppen)



## References

* Bak, P., & Sneppen, K. (1993). Punctuated equilibrium and criticality in a simple model of evolution. Physical review letters, 71(24), 4083.
* Fraiman, D. (2018). Bak-Sneppen model: Local equilibrium and critical value. Physical Review E, 97(4), 042123.
* Paczuski, M., Maslov, S., & Bak, P. (1996). Avalanche dynamics in evolution, growth, and depinning models. Physical Review E, 53(1), 414.
