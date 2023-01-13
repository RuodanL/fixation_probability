# Fixation probability on hypergraphs with high symmetry

Here we provide the code that accompanies the following paper:

R. Liu, N. Masuda.
Fixation dynamics on hypergraphs.
Preprint on arXiv: [I will update here when I get the arXiv identifier number](pending URL)

## Jupyter notebook to create figures

```
figure_X.ipynb
```
contains the code to produce fig X in the manuscript, where $X = 2, 4, 5$, and E. $6$.
The folder [notebooks/data/](https://github.com/RuodanL/fixation_probability/tree/main/notebooks/data) contains the numerical results based on which Figures 2, 4, 5, and E.6 are produced.

## Python code for running evolutionary dynamics on empirical hypergraphs

The following Python codes are also included:

- bd-on-hypergraph.py -- to simulate the birth-death process on empirical networks and randomized networks.
- bd-on-one-mode-projection.py -- to simulate the birth-death process on one-mode projection of the empirical network.
